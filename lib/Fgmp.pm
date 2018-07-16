package Fgmp;
use strict; 
use IO::All; 
use Bio::SeqIO;
use Carp;
use Data::Dumper;
use feature 'say';  
use Getopt::Long;
use IPC::Open3 qw ( open3 );
use IPC::Run qw ( run start pump timer timeout );
use List::Util qw(max min);
use Scalar::Util qw(looks_like_number);

our ($VERSION, $DEBUG, $CALLER);
$VERSION = '1.0';

sub load_paths {
	my ($fmpdir,$wrkdir,$tmpdir) = ("","","");
	my $paths = io(@_);
	   $paths->autoclose(0);
	   while (my $lp = $paths->getline || $paths->getline){
	   chomp $lp;
		next if $lp =~ m/^#/;
		   ($fmpdir) = $lp if ( $lp =~ m/FGMP=/);
		    $fmpdir =~s/FGMP=//;
		   ($wrkdir) = $lp if ( $lp =~ m/WRKDIR=/);
		    $wrkdir =~s/WRKDIR=//;
		   ($tmpdir) = $lp if ( $lp =~ m/TMP=/);
		    $tmpdir =~s/TMP=//;
	}
	return($fmpdir,$wrkdir,$tmpdir);
}

sub count_num_of_seqs {
	my ($file) = @_; 
	my $num = `grep -c \">\" $file`;
	chomp($num);
	return($num);
}

sub how_many_locus {
	my ($in) = @_; 
	my $num1 = `grep -c LOCUS $in`;
	chomp($num1);
	return($num1)
}

sub execute {
        my ($cmmd) = @_;
        #warn "###\t$cmmd\n";
        system($cmmd)==0 || croak "cannot execute\t$cmmd:$!\n";
}

sub clean_files {
	my ($files,$directory,$tag) = @_; 
	
	my @listOfFiles = @$files;
	my $dir = $$directory;
	my $Tag = $$tag;

	if ($Tag eq 'FALSE'){
		foreach my $file (@listOfFiles){
			if (-e $file){
				execute("unlink $file");
			}
		}
	} else {
		execute("mkdir -p $dir");
			foreach my $file (@listOfFiles){
				if (-e $file){
					execute("mv $dir");
				}
		}
	}
}

sub split_and_run_sixpack {
	my ($multifasta) = @_; 

	#execute("csplit -s -f $multifasta.chunk_6p -z $multifasta \'/^>/\' \'{*}\'");
	execute("sed \'s/>/>multi_/g\' $multifasta > $multifasta.tmp");
	execute("sed \'s/\|/_/g\' $multifasta.tmp > $multifasta.tmp2");
	execute("mv $multifasta.tmp2 $multifasta.tmp");
	execute("seqretsplit -sequence $multifasta.tmp -outseq multi");
	execute("rm -rf $multifasta.sixpack_tmp") if (-d "$multifasta.sixpack_tmp");
	execute("mkdir $multifasta.sixpack_tmp"); 
	execute("mv multi*.fasta $multifasta.sixpack_tmp");
	#execute("mv $multifasta.chunk_6p\* $multifasta.sixpack_tmp");
	my $io = io("$multifasta.sixpack_tmp");
	my @contents = $io->all; 
	
	my $chk = "";
	foreach $chk (@contents) {
		execute("sixpack -sequence $chk -outfile $chk.sixpack -outseq $chk.orfs -orfminsize 50 -verbose false");
	}

	# concatenate
	execute("cat $multifasta.sixpack_tmp/\*.orfs > $multifasta.orfs");
	execute("rm -rf $multifasta.sixpack_tmp $multifasta.tmp");
}


sub multithread_exonerate {
	my ($candidateFasta, $cpuAvail, $proteins, $srcdir,$outdir) = @_;

	my $elementTmp = ""; 
	my @allChunks = ();
	my @runs_fas = ();
	my @run_exo = (); 

	my @cmd = ();

	my @seqs = extractSeqname($candidateFasta);
	my $numOfseqs = scalar(@seqs);
	my $numOfChunks = $cpuAvail;
	my $numOfseqPerChunks = sprintf("%.0f",($numOfseqs / $numOfChunks)); 
	
	my $i = 0;
	while ((@seqs) > 1){
		my @chunk = splice(@seqs,1,$numOfseqPerChunks);
		my $chunk = join("\t", @chunk);
		push(@allChunks,$chunk);		
		}

	# catch the remaining in case of impair number
		my $impair = join("",@seqs);

	# now adding the first element of the list
	unless (!defined($impair)){	
		# replacing the first element
		my $newFirstElement = "$allChunks[0]\t$impair";
		$allChunks[0] = $newFirstElement;
	}
	
	# now create chunk files
	my $chk = ""; 
	my $count = 0;

	my @chunksForFastaExtr = (); 
	 
	# each element of @allChunks is >sca1\t$scaf2 etc...
	foreach $chk (@allChunks){
		my @dataC = split /\t/, $chk;
		
		my $chktmp = "";
		my $chkind = ""; 
		foreach $chkind (@dataC){
			$chkind =~s/>//;
			$chktmp .="$chkind\n";
		}
		io("$candidateFasta.chunk$count.tmp")->write($chktmp);
		push(@chunksForFastaExtr,"$candidateFasta.chunk$count.tmp"); 
		$count++;
	}
	
	# now extract fasta seq for these chunk 
	# and launch exonerate
	my $ext = ""; 
	foreach $ext (@chunksForFastaExtr){
		# You need to update the path later
		push(@runs_fas,"perl $srcdir/retrieveFasta.pl $outdir/$ext $outdir/$candidateFasta > $outdir/$ext.fas");
		push(@run_exo,"exonerate --model protein2genome --percent 5 -q $proteins --showtargetgff Y -t $outdir/$ext.fas --showvulgar F --showalignment T --ryo \'%qi,%ql,%qab,%qae,%ti,%tl,%tab,%tae,%et,%ei,%es,%em,%r,%pi,%ps,%C\' > $outdir/$ext.p2g");

	}
	return($numOfseqs,$numOfChunks,$numOfseqPerChunks,\@runs_fas,\@run_exo);
}
sub multithread_augustus {
        my ($candidateFasta, $cpuAvail, $proteins, $srcdir, $augdir,$outdir,$speciestag) = @_;

        my $elementTmp = "";
        my @allChunks = ();
        my @runs_fas = ();
        my @run_aug = ();
	my @run_gff2aa = ();
	my @toconcat = ();
	
        my @cmd = ();

        my @seqs = extractSeqname($candidateFasta);
        my $numOfseqs = scalar(@seqs);
        my $numOfChunks = $cpuAvail;
        my $numOfseqPerChunks = sprintf("%.0f",($numOfseqs / $numOfChunks));

        my $i = 0;
        while ((@seqs) > 1){
                my @chunk = splice(@seqs,1,$numOfseqPerChunks);
                my $chunk = join("\t", @chunk);
                push(@allChunks,$chunk);
                }

        # catch the remaining in case of impair numer
                my $impair = join("",@seqs);

        # now adding the first element of the list
        unless (!defined($impair)){
                # replacing the first element
                my $newFirstElement = "$allChunks[0]\t$impair";
                $allChunks[0] = $newFirstElement;
        }

        # now create chunk files
        my $chk = "";
        my $count = 0;

        my @chunksForFastaExtr = ();

        # each element of @allChunks is >sca1\t$scaf2 etc...
        foreach $chk (@allChunks){
                my @dataC = split /\t/, $chk;

                my $chktmp = "";
                my $chkind = "";
                foreach $chkind (@dataC){
                        $chkind =~s/>//;
                        $chktmp .="$chkind\n";
                }
                io("$candidateFasta.chunk$count.tmp")->write($chktmp);
                push(@chunksForFastaExtr,"$candidateFasta.chunk$count.tmp");
                $count++;
        }

        # now extract fasta seq for these chunks
        my $ext = "";
        foreach $ext (@chunksForFastaExtr){
		#push(@run_aug,"$augdir/bin/augustus --species=$speciestag --AUGUSTUS_CONFIG_PATH=$augdir/config $outdir/$ext.fas > $outdir/$ext.gff");
		execute("perl $srcdir/retrieveFasta.pl $ext $candidateFasta > $outdir/$ext.fas");   		
		push(@run_aug,"augustus --species=$speciestag --AUGUSTUS_CONFIG_PATH=$augdir/config $outdir/$ext.fas > $outdir/$ext.gff");

		#push(@run_gff2aa,"$augdir/scripts/getAnnoFasta.pl $outdir/$ext.gff");
		push(@run_gff2aa,"getAnnoFasta.pl $outdir/$ext.gff");
		push(@toconcat,"$outdir/$ext.aa");
	 }
	return($numOfseqs,$numOfChunks,$numOfseqPerChunks,\@run_aug,\@run_gff2aa,\@toconcat);
}

sub extractSeqname{
	my ($file) = @_; 
	
	my @list = (); 
	my $in = io("$file"); 
	   $in->autoclose(0); 
       while(my $line = $in->getline || $in->getline){
	chomp $line; 
		if ($line =~ m/^>/){
			# becarefull, you need to propely format the contig name
			push(@list,$line);
		}
	} 
	return (@list);
}

sub execute_and_returnWhendone {
	my (@jobs) = @_; 

	sub launch {
		open (local *CHILD_STDIN, '<', '/dev/null') || croak $!;
		return open3('<&CHILD_STDIN','>&STDOUT','>&STDERR',@_);
	}	

	my %children = ();
	for my $cmd (@jobs){
		#warn "Command: $cmd started at ".localtime."\n";
		my $pid = launch($cmd);
		$children{$pid} = $cmd;
	}
	
	while (%children){
		my $pid = wait();
		die $! if $pid < 1;
		my $cmd = delete($children{$pid});
		#warn "Command: $cmd ended at ".localtime."with \$? = $?.\n";
	}
	my $st = keys (%children);
	return($st);
}

sub extractCandidateRegion {
	my ($tblastn) = @_; 

	my %h1 = ();
	my @locations = (); 
	
	my $tb = io("$tblastn");
	   $tb->autoclose(0); 
	   while( my $tbline = $tb->getline || $tb->getline){
	   chomp $tbline; 
		next if $tbline =~ m/^#/;
		my @datatb = split /\t/, $tbline; 
		my ($target,$sstart,$send) = ($datatb[0],$datatb[1],$datatb[2]); 
		
		if ($sstart > $send) { # reverse
			my $rsstart = $send;
			my $rsend   = $sstart;
			push(@locations,$rsstart,$rsend);
			@{$h1{$target}} = @locations;
		} else {
			push(@locations,$sstart,$send); 
			@{$h1{$target}} = @locations;
		}
		
		# See it before
		if (exists ($h1{$target})){
			# add new HSPs
			my @range = @{$h1{$target}};
			# update
			push(@range, $sstart,$send);
			@{$h1{$target}} = @range
		} else {
			@{$h1{$target}} = @locations;
		}
		@locations = (); 
	}

	# compute the borders
	my %boundaries = ();
	my @bounds = (); 

	my $el = ""; 
	foreach $el (keys %h1){
		my $startBoundaries = "";
		
		# meaning that 5000 (5kb) extends the end of the contig
		if (min(@{$h1{$el}}) < 5000000){
			$startBoundaries = 0;
		} else {
			$startBoundaries = (min(@{$h1{$el}}) - 5000000);
		}
		
		# you might need the size of the contig
		my $endBoundaries = (max(@{$h1{$el}}) + 5000000);

		# populating the hash
		push(@bounds,$startBoundaries,$endBoundaries);		 
		@{$boundaries{$el}} = @bounds;

		# cleaning @bounds 
		@bounds = ();
	} 
	return (%boundaries);
}

sub exportCandidateRegions {
	my ($info,$fasta) = @_; 

	my %info = %$info;

	my $buff = ""; 
	my $seqio = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
		while(my $seq = $seqio->next_seq) {
  		my $name = $seq->id;
		my $sequence = $seq->seq;
	
		if (defined ($info{$name})){
			my @coord = @{$info{$name}};
			# If 5kb is longer the contig, substr simply goes reverse
			my $chunk = substr($sequence,$coord[0],$coord[1]);
			$buff .=">$name\n";
			$buff .="$chunk\n";
		}
	}	
	io("$fasta.candidates.fa")->write($buff);
}

sub run_on_single_node {
	my ($genomeFas,$protFas,$direct) = @_; 
	execute("exonerate --model protein2genome --percent 5 -q $protFas --showtargetgff Y -t $direct/$genomeFas --showvulgar F --showalignment T --ryo \'%qi,%ql,%qab,%qae,%ti,%tl,%tab,%tae,%et,%ei,%es,%em,%r,%pi,%ps,%C\' > $direct/$genomeFas.p2g"); 
}
sub search_in_reads {
	my ($reads,$protein,$fgmpdir,$threads) = @_;

	my @markers = load($protein);	
	
	my $ids = `grep \'^>\' $reads`;
	   $ids =~s/\n//g;
	   
	my @ids = split />/,$ids;
	my $num_seqs = scalar @ids;

	my %samples = reservoir_sampling(@ids);

	my %markersFound = ();
	my $it = "";
	my $trials = 0;

	foreach $it ( keys %samples) {
        	my @sample = @{$samples{$it}};
        	generateFasta(\@sample,$it,$reads,\@ids,$fgmpdir);  # generate for blast $it.fa
        	
		my %makers = runBlastx("$reads.sampled.$it.fa",$protein,$threads);
        	my %new = compare(\%makers, \%markersFound);
		
		my $new = scalar (keys %new);
        	my $previous = scalar (keys %markersFound);

        	if ($new < ($previous * 0.1)){
                	$trials++;
                	warn"...only $new markers found - threshold not satisfied - attempt # $trials / 20\n";
        	} else {
                	warn"...new markers found:\t$new (previous $previous)";
      		  }
        	# update %markerFound
        	# becareful will erase previous data but
        	my $t = "";
		foreach $t ( keys %makers ) {
                	$markersFound{$t} = $makers{$t}
        		}
        	last if ($trials == 20);
	}
	my $founds = scalar (keys %markersFound);
	return($founds);
}

sub compare {
        my ($h1,$h2) = @_;

        my %h1 = %$h1;
        my %h2 = %$h2;
        my %new = ();
        my ($m1,$m2) = ("","");
        foreach $m1 ( keys %h1 ) {
                unless ( $h2{$m1} ) {
                        $new{$m1} = $h1{$m1};
                }
        }
        return(%new);
}

sub runBlastx {
        my ($query,$target,$threads) = @_;
        execute("blastx -db $target -query $query -evalue 0.01 -num_threads $threads -outfmt 6 -max_target_seqs 1 -out $query.blastx");
	my %founds = extractMarkers("$query.blastx");
        return(%founds);
}
sub extractMarkers {
        my %h = ();
        my $blastx = io(@_);
           $blastx->autoclose(0);
           while ( my $bl = $blastx->getline || $blastx->getline ) {
           chomp $bl;
                my @dat = split /\t/,$bl;
                my ($q,$s,$eval) = ($dat[0],$dat[1],$dat[10]);

                if ($h{$s}) {
                        my @tmp = @{$h{$s}};
                        unless ($eval > $tmp[1]){
                                @{$h{$s}} = ($q,$eval);
                        }
                } else {
                        @{$h{$s}} = ($q,$eval);
                }
        }
        return(%h);
}
sub generateFasta {
        my ($list,$iteration,$readsFile,$ids,$dir) = @_;

        my @ids = @$ids;

        my $buffer = "";
        my $id = "";
        foreach $id ( @$list ) {
		my @clean = split /\s+/,$ids[$id];
                #$buffer .="$ids[$id]\n";
		$buffer .="$clean[0]\n";
        }
        io("$iteration.tmp")->write($buffer);
        execute("perl $dir/src/retrieveFasta.pl $iteration.tmp $readsFile > $readsFile.sampled.$iteration.fa");
        #execute("rm $iteration.tmp");
}
sub reservoir_sampling {
        my (@list) = @_;
        
	my %samples = ();

        my $num_of_samples = '1000';
        my $sample_size = '10000';
        my $counter = 1;
	my $num_seqs = scalar (@list);

        while ( $counter < ($num_of_samples + 1)){
         my @sampled = (1 .. $sample_size);
         for my $j ( $sample_size + 1 .. $num_seqs) {
                 $sampled[ rand(@sampled) ] = $j if rand() <= ($sample_size / $j)
                }

        my $tmp = "";
        my @buf = ();
        foreach $tmp ( @sampled ) {
                push(@buf, $tmp);
        }

        @{$samples{"s.$counter"}} = @buf;
        $counter++;
        }
        return(%samples);
}

sub load {
        my @array = ();
        my $f = io(@_);
          $f->autoclose(0);
          while ( my $fl = $f->getline || $f->getline){
          chomp $fl;
                if ( $fl =~ m/^>/){
                        $fl =~s/>//;
                        push(@array, $fl);
                }
        }
        return(@array);
}
sub check_multicopies {
	my ($hmmreportFile,$multicopiesFile) = @_; 
	
	my %copiesWithLowerThanExpected = ();
	my %multicopies = loadcsv($multicopiesFile);

	my $buffHmms = "";
	my $m = "";
	foreach $m ( keys %multicopies ) {
		my $n = how_many_copies($m,$hmmreportFile);
		if ( $n < $multicopies{$m} ){
			$buffHmms .="$m\tdetected:\t$n,\texpected:\t$multicopies{$m}\tLOW\n";
		} elsif ( $n > $multicopies{$m} ) {
			#$buffHmms .="$m\tfound:\t$n,expected:\t$multicopies{$m}\tHIGH\n";
		} else {
			next;
		}		
	}
	io("$hmmreportFile.multicopies_check.csv")->write($buffHmms);	
}

sub how_many_copies {
	my ($hmmtofind,$hmmfile,$report) = @_;
	my @arrhmms = ();
	my $hm1 = io($hmmfile);
	   $hm1->autoclose(0);
	   while ( my $hml1 = $hm1->getline || $hm1->getline ) {
	   chomp $hml1; 
		next if $hml1 =~ m/^#/;
        	next if $hml1 =~ m/^$hmmtofind/; # those are proteins

		if ( $hml1 =~ m/$hmmtofind/ ) {
			my @dathm1 = split /\s+/, $hml1;
			my ($prot,$hm,$eval) = ($dathm1[0],$dathm1[3],$dathm1[6]);
		               if ( $eval < 1e-20 ) {
                	                push(@arrhmms,$prot);
                                }
		} else {
			next;
		}
	}
	my %uniqProtein = ();
	my $uni = "";
	foreach $uni ( @arrhmms ) {
		$uniqProtein{$uni}++;
	}
	return(scalar(keys %uniqProtein ));
}
sub loadcsv {
	my %multi = ();
	my $csv = io(@_);
	   $csv->autoclose(0);
	   while ( my $csvl = $csv->getline || $csv->getline ) {
	   chomp $csvl; 
		my @csvdat = split /,/, $csvl;
	        $multi{$csvdat[0]} = $csvdat[1];
	} 
	return(%multi);
}
sub retrieveFasta {
	my ($list,$proteomdb,$extractedSeqs) = @_;

	my %id_hash = ();

	my $ids = io($list);
	   $ids->autoclose(0);
	   while(my $idsl= $ids->getline || $ids->getline){
	   chomp $idsl;
		$id_hash{$idsl} = 1;
	}

	my $buf = "";
	my $new=Bio::SeqIO->new(-file=> $proteomdb, -format=>"fasta");
	while (my $seq=$new->next_seq){
        	if (defined $id_hash{$seq->id}){
        	$buf .=">$seq->id\n$seq->seq\n";
        	}
	}
	io("$extractedSeqs")->write($buf);
}
sub filter_unfiltByScore {
	my ($hmmsearchOut,$cutoff,$consMarkers,$tag,$nhmmerOut,$fucesNames,$readsFile,$multicopiesSearch) = @_;
	
#	die "
#	# This scipt takes a hmmsearch output ( domtbl format) file, a profile thresholds file and
#	#  a file containing the name of 245 fgmp markers ( one per line)
#	# and give a text file containing the name of the maker and the best scoring
#	# protein ( one per line)
#	# filter_unfiltByScore.pl < hmmseach > < profiles_cutoff.tbl > < 593 FGMP> < TAG > < nhmmer output > < UCEs prefix >> file
#	#
#	#
#	# options
#	# -v              verbose
#	# --cutoff        alternative cutoff
#	# 
#	" unless (@ARGV > 2);
#
	 my $completeCutoff  = '0.7';
#
#	 GetOptions(
#             "cutoff=f" => \$completeCutoff
#          );

# ---------------------------------------------- #
#
# Check if the files exist
#
# ---------------------------------------------- #
if (!(-e "$hmmsearchOut")){
        say "File does not exist: $hmmsearchOut";
        exit(1);
}

if (!(-e "$cutoff")){
        say "File does not exist: $cutoff";
        exit(1);
}

if (!(-e"$consMarkers")){
        say "File does not exist: $consMarkers";
        exit(1);
}

if (!(-e "$nhmmerOut")){
        warn"File does not exist: $nhmmerOut\n";
        exit(1);
}

if (!("$fucesNames")){
        warn"Cannot find fUCEs names: $fucesNames\n";
        exit(1);
}

if (!("$multicopiesSearch")){
        warn"Cannot find multicopies genes search file: $multicopiesSearch\n";
        exit(1);
}
#---------------------------------------------- #
# parse nhmmer data
# --------------------------------------------- #
# return how many fUCEs are found in the genome
my ($fucesFound,$fucesMissing,$totalfuces) = parse_nhmerData($nhmmerOut,$fucesNames);
#
my %fucesFound = %$fucesFound;
my %fucesMissing = %$fucesMissing;
my $allfuces = $$totalfuces;

# ---------------------------------------------- #
#  load multicopies genes data
#  ---------------------------------------------- #

my @multicopiesData = loadMulticopiesData($multicopiesSearch);

# ---------------------------------------------- #
#   load search in reads data
# ---------------------------------------------- #

my $reads = loadreadsData($readsFile);

# ---------------------------------------------- #
# parse cutoff file and record model score cutoff
# and length
# ---------------------------------------------- #
# record hmm score and the length of the original
# untrimmed alignment
#

&report("...LOADING PROFILES DATA");

my ($hmmcutoff,$hmmAlgnlen,$hmmRatio) = parseCutoffData($cutoff);
my %hmmcutoff = %$hmmcutoff;
my %hmmAlgnlen = %$hmmAlgnlen;
my %hmmRatio = %$hmmRatio;

# load the name of markers
&report("...LOADING markers data");

my %markers = loadMarkerName("$consMarkers");


# load sequence length
#
&report("...LOADING SEQS LENGTH AND MODELS TO SEARCH");

my ($seqLen, $hmmName) = extractSeqLenAndmodelNames($hmmsearchOut);
my %seqLen = %$seqLen;
my %hmmName = %$hmmName;

# ---------------------------------------------- #
# retrieve the best protein for each model
# based on score
# ---------------------------------------------- #
#
&report("...ANALYZING HMMSEARCH OUTPUT");

# containers
my %seq2scoreG = ();
my %seq2model = ();
my %significantmodel = ();

my $model = "";
foreach $model (keys %hmmName){
	
	# how many proteins map against this model
	# key = seq, value = score
	
	my %seq2scorel = extractSeqForThismodel($model,$hmmsearchOut);
	
	# now pick for this model the best candidate
	my $cand = "";
        my $scoretmp = 10;
        my $bestCand = "";
        my $empty = "";

        foreach $cand (keys %seq2scorel){
                unless ($seq2scorel{$cand} eq 'NOT FOUND') { # skip those models without hits

		 # check that the best seq pass the filtering cuoff
		 if ($seq2scorel{$cand} >= ($hmmcutoff{$model} * 0.3)){
			if ($markers{$model}){
				$significantmodel{$model} = $model;
			} else {
				#
			}
		# now iterate to pick the best seq based on score
		 unless ($seq2scorel{$cand} < $scoretmp){
			$scoretmp = $seq2scorel{$cand};
                                        $bestCand = $cand;
                                        }
                        } else {
				# do something the seqs that failed to pass the score cutoff
			}
                }
        }
        $seq2model{$bestCand} = $model if ($bestCand ne $empty);
        $seq2scoreG{$bestCand} = $scoretmp if ($scoretmp ne $empty); # record the score of those seq
	}
	my $fgmpMarkers = keys %significantmodel;

	# catch missing markers
	my %missingFGMPmarkers = ();

	my $m = "";
	foreach $m (keys %markers){
        	unless ($significantmodel{$m}){
                	$missingFGMPmarkers{$m} = $m;
        	}
	}


# ---------------------------------------------- #
# EXTRACT ALN LENGTH FOR BEST CANDIDATES
# ---------------------------------------------- #
&report("...EXTRACTING ALN LENGTH FOR THE BEST CANDIDATE");
my %seq2alnLen = ();
my $el = "";
foreach $el (keys %seq2model){
        my $alnLength = extractLenAln($el,$seq2model{$el},$hmmsearchOut);
        $seq2alnLen{$el} = $alnLength;
}

# ---------------------------------------------- #
# TRACK ABBERANT PROTEINS
# ---------------------------------------------- #
&report("...TRACKING ABBERANT PROTEINS");

# key = seq xxx = value = AB|NORM
my ($abberantProteins,$countAbb) = isAbberrant(\%seqLen,\%hmmAlgnlen,\%seq2model);
my %abberantProteins = %$abberantProteins;

# ---------------------------------------------- #
# check that protein cover at least 70% of the
# alignment (i.e. the original untrimmed aligment
# used to build the profile)
# ---------------------------------------------- #
#
&report("COMPLETENESS ESTIMATION");
my %goodPreds = ();
my $complete = 0;
my $partial = 0;

# counters
my ($v) = "";
my $count = 1;

my $total = scalar (keys %seq2model);

	# at this stage of the model present in this hash have a sequence that has passed the score filter
	# only the best sequence per model should arrive here

	my $logfull  = "# SEQ\tMODEL\tSTATUS\tCOMPLETNESS\tCONS\n";
	   $logfull .= "# ------------------------------------------------------------ #\n";

	foreach $v (keys %seq2model){
        	&report("----> JOBS DONE $count / $total");
        	my $threshold = ($hmmAlgnlen{$seq2model{$v}} * $completeCutoff);
	 # check that the length of the sequence is at least 70% that of the untrimmed alignment used to build the profil
	 if ($seq2alnLen{$v} >= $threshold){
                $complete++;
	 	# check if that is one of the completeness markers
	 	if ($significantmodel{$seq2model{$v}}){
                        $logfull .="$v\t$seq2model{$v}\tCOMPLETE\tCOMPL_MARK\t$abberantProteins{$v}\t$hmmRatio{$seq2model{$v}}\n";
                } else {
                         $logfull .="$v\t$seq2model{$v}\tCOMPLETE\tXXXX\t$abberantProteins{$v}\t$hmmRatio{$seq2model{$v}}\n";
                }

        } elsif ($seq2alnLen{$v} <  $threshold){
                if ($significantmodel{$seq2model{$v}}){
                        $logfull .="$v\t$seq2model{$v}\tPARTIAL\tCOMPL_MARK\t$abberantProteins{$v}\t$hmmRatio{$seq2model{$v}}\n";
                    } else {
                        $logfull .="$v\t$seq2model{$v}\tPARTIAL\tXXXX\t$abberantProteins{$v}\t$hmmRatio{$seq2model{$v}}\n";
                }
                $partial++;
        } else {
		# 
	}
        $count++;
}
&report("DONE!");
&report("PRINTING FULL REPORT");
io("$hmmsearchOut.full_report")->write($logfull);

&report("PRINTING SUMMARY REPORT");
my $totalUnfilteredPreds = keys %seqLen;
my $totalSeq  = keys %seq2model;
my $totalFgmp = scalar (keys %markers);
my $completness = sprintf("%.1f",(($fgmpMarkers / $totalFgmp) * 100));
my $missingFGMPmarkers = scalar %missingFGMPmarkers;
my $ucesFound = scalar(keys %fucesFound);
my $ucesFoundPercent =   sprintf("%.1f", (($ucesFound / $allfuces) * 100));
my $markersInReads = "";
if ($reads eq 'NA'){
	$markersInReads = 'NA';
} else {
	$markersInReads = sprintf("%.1f",(($reads / 593) * 100));
}
# printing report
        my $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tFungal Highly conserved elements (FHCEs) report\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tNUM OF UCEs found:\t$ucesFound out of $allfuces\n";
           $logsum .= "|\tUCEs completion estimation:\t($ucesFoundPercent\%)\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\t### MISSING fUCEs\n\n";

           foreach my $misfucs (keys %fucesMissing){
                        $logsum .="|\t$misfucs\n";
           }

           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tSUMMARY STATISTICS\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tTOTAL NUM OF PREDS ANALYZED:\t$totalUnfilteredPreds\n";
           $logsum .= "|\tNo. BEST FILTERED PREDS:\t$totalSeq\n";
           $logsum .= "|\tNUM OF MARKERS COMPLETE:\t$complete\n";
           $logsum .= "|\tNUM OF MARKERS PARTIAL:\t$partial\n";
           $logsum .= "|\tNUM OF ABBERANT PROTEINS:\t$countAbb\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tSTATUS\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tNUM OF FGMP COMPLETENESS MARKERS:\t$fgmpMarkers\n";
           $logsum .= "|\tEstimate completeness:\t$completness (%)\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\t### MISSING MARKERS\n\n";

           foreach my $mis (keys %missingFGMPmarkers){
                        $logsum .="|\t$mis\n";
           }

           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tSUMMARY OF READS ANALYSIS\n";
           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tNumber of markers detected in reads:\t$reads\n";
           $logsum .= "|\tEstimate completeness:\t$markersInReads (%)\n";

           $logsum .= "| ---------------------------------------------------------- |\n";
           $logsum .= "|\tSUMMARY OF MULTICOPY GENES ANALYSIS\n";
           $logsum .= "| ---------------------------------------------------------- |\n";

           $logsum .= "|\t### MULTICOPIES GENES WITH LOW COPIES NUMBER\n";

           if ( @multicopiesData == 0 ) {
                $logsum .= "|\t33 MULTICOPY GENES ARE IN MULTICOPIES\nno evicence of collapsed regions\n";
         } else {
           foreach my $mult ( @multicopiesData ) {
                        if ( $mult eq '' ){
                                $logsum .= "|\t33 MULTICOPY GENES ARE IN MULTICOPIES\n\tNo evicence of collapsed regions\n";
                        } else {
                                $logsum .="|\t$mult\n";
                        }
                }
        }
 $logsum .= "| ---------------------------------------------------------- |\n\n";
           $logsum .= "| These results are based on the set of genes selected by OHC & JES #\n\n";
           $logsum .= "| Key:\n";
           $logsum .= "| Proteins = 593 conserved fungal genes\n";
           $logsum .= "| DNA = $allfuces fungal highly conserved elements\n";
           $logsum .= "| \%completeness = percent of 593 FCGs in the dataset\n";
           $logsum .= "|---------------------------------------------------------- |\n";


         io("$hmmsearchOut.summary_report")->write($logsum);
}

sub extractLenAln {
        my ($seqname, $modelname,$filehmmseach) = @_;

        my $len = 0;
        my $hse = io("$filehmmseach");
           $hse->autoclose(0);
           while(my $hseline = $hse->getline || $hse->getline){
           chomp $hseline;
                next if $hseline =~ m/^#/;
                if (($hseline =~ m/$seqname/) && ($hseline =~ m/$modelname/)){
                        my @dat = split /\s+/, $hseline;
                        my ($envStart,$envEnd) = ($dat[19],$dat[20]);
                        my $lenTmp = ($envEnd - $envStart) + 1;
                        $len += $lenTmp;
                }
        }
        return $len;
}
sub isAbberrant{
        my ($seqlen,$hmmlen,$seq2mod) = @_;

        my %seqslen = %$seqlen;
        my %hmmlen = %$hmmlen;
        my %seq2mod = %$seq2mod;

        my $abcount = 0;
        my %abb = ();

	# the names of seqs
	my @seqs = keys %seq2mod;
	
	# find
	# my $id = "";
	my $empty = "";
	foreach my $id (keys %seq2mod){
                if ($seqslen{$id} > ($hmmlen{$seq2mod{$id}} * 2)){
                                $abb{$id} = 'ABERRANT';
                                $abcount++;
                        } else {
				$abb{$id} = 'NORMAL';
                        }
        }
        return(\%abb,$abcount);
}
sub extractSeqLenAndmodelNames {
        my ($inHmmseach) = @_;

        my %h4 = ();
        my %h5 = ();

        my $inH = io("$inHmmseach");
           $inH->autoclose(0);
           while (my $inL = $inH->getline || $inH->getline){
           chomp $inL;
                        next if $inL =~ m/^#/;
                        my @dataL = split /\s+/, $inL;
                        my ($id,$len,$modelName) = ($dataL[0],$dataL[2],$dataL[3]);
                        $h4{$id} = $len;
                        $h5{$modelName} = 'x';
           }
        return(\%h4,\%h5);
}
sub loadMarkerName {
        my ($data) = @_;
        my %h1 = ();
        my $in = io("$data");
           $in->autoclose(0);
           while ( my $ln = $in->getline || $in->getline){
           chomp $ln;
                my @dat1 = split /_/, $ln;
                $h1{$dat1[0]} = $dat1[0];
            }
        return (%h1);
}

sub extractSeqForThismodel {
        my ($mod,$hmmscanFile) = @_;

        my %h = ();
        my $file = io("$hmmscanFile");
        $file->autoclose(0);
        while (my $f = $file->getline || $file->getline){
        chomp $f;

        next if $f =~ m/^#/;
        if ($f =~ m/$mod/){
                my @inline = split /\s+/, $f;
                my ($seqname,$model,$score,$envStarfileenvEnd) = ($inline[0],$inline[3],$inline[7],$inline[19],$inline[20]);
		if ($model eq $mod){
			# the score of a given seq against a particular model does not change
			$h{$seqname} = $score;
                        }
                else {
                        next;
                }
            }
        }
        my @found = keys %h;
	if (@found == 0){ # check if empty hash
 	$h{$mod} = 'NOT FOUND';
                }
        return(%h);
}

sub parseCutoffData{
        my ($data) = @_;

        my (%h1,%h2,%h3) = ();

        my $in = io("$data");
        $in->autoclose(0);
        while (my $line = $in->getline || $in->getline){
        chomp $line;
                my @dat = split /\t/, $line;
                $h1{$dat[0]} = $dat[2];  # hmm score cutoff
                $h2{$dat[0]} = $dat[3];  # untrim aln length used for hmm construction
                $h3{$dat[0]} = $dat[4];  # cons ratio
        }
        return (\%h1,\%h2,\%h3);
}
sub report {
        my ($log) = @_;
        warn"--- LOG:\t$log\n";
}

sub parse_nhmerData {
	my ($report,$names) = @_;

        my (%found,%missing) = ((),());

        # read all names
        my %all = ();
        my $n = io($names);
           $n->autoclose(0);
           while ( my $nl = $n->getline || $n->getline ) {
           chomp $nl;
                $all{$nl} = $nl;
        }

        my %detected = ();
        my $nh = io($report);
           $nh->autoclose(0);
           while ( my $nhl = $nh->getline || $nh->getline ) {
           chomp $nhl;
                next if $nhl =~ m/^#/;
                #my ($g) = $nhl =~/\s+(FG\d+\.1)\s+/;
                my @nd = split /\s+/, $nhl;
                my $g = $nd[2];
                $detected{$g} = 'x';
        }

	# soustract
	my $i = "";
	
	foreach $i ( keys %all ) {
		if ($detected{$i}){
			$found{$i} = $i;
		} else {
			$missing{$i} = $i;
		}
	}
	my $totalfucs = scalar (keys %all);
	return(\%found,\%missing, \$totalfucs);
}
sub loadMulticopiesData {
        my @mData = "";
	my $mucop = io(@_);
	$mucop->autoclose(0);
	while ( my $mucopl = $mucop->getline || $mucop->getline ) {
	chomp $mucopl;
		push(@mData,$mucopl);
	}
	return(@mData);
}
sub loadreadsData {
        my $rf = io(@_);
           $rf->autoclose(0);
            while (my $rfl = $rf->getline || $rf->getline ){
            chomp $rfl;
                return($rfl);
        }
}
1;

