package Fgmp;
use strict; 
use IO::All; 
use Carp;
use feature 'say';  
use IPC::Open3 qw ( open3 );
use IPC::Run qw ( run start pump timer timeout );
use Data::Dumper;
use List::Util qw(max min);
use Bio::SeqIO;


sub load_paths {
	my ($fmpdir,$wrkdir,$tmpdir) = ("","","");
	
	my $paths = io(@_);
	   $paths->autoclose(0);
	   while (my $lp = $paths->getline || $paths->getline){
	   chomp $lp;
		next if $lp =~ m/^#/;
		   ($fmpdir) = $lp =~/FGMP=(.*)/;
		   ($wrkdir) = $lp =~/WRKDIR=(.*)/;
		   ($tmpdir) = $lp =/TMP=(.*)/;
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
        warn "###\t$cmmd\n";
        system($cmmd)==0 || croak "cannot execute\t$cmmd:$!\n";

}
sub clean_files {
	my (@files, $tag,$directory) = @_; 

	execute("mkdir -p $directory/$tag-temporyfiles");

		foreach my $file (@files){
			if ( -e $file){
				execute("mv $file $directory/$tag-temporyfiles");
			}
		# compress
		#execute("tar -cvzf $tag-temporyfiles.tar.gz $directory/$tag-temporyfiles");
		#execute("rm -rf $directory/$tag-temporyfiles"); 
	}
}
sub split_and_run_sixpack {
	my ($multifasta) = @_; 

	execute("csplit -s -f $multifasta.chunk_6p -z $multifasta \'/^>/\' \'{*}\'");
	execute("rm -rf $multifasta.sixpack_tmp") if (-d "$multifasta.sixpack_tmp");
	execute("mkdir $multifasta.sixpack_tmp"); 
	execute("mv $multifasta.chunk_6p\* $multifasta.sixpack_tmp");
	my $io = io("$multifasta.sixpack_tmp");
	my @contents = $io->all; 
	
	my $chk = "";
	foreach $chk ( @contents) {
		execute("sixpack -sequence $chk -outfile $chk.sixpack -outseq $chk.orfs -orfminsize 50");
	}

	# concatenate
	execute("cat $multifasta.sixpack_tmp/\*.orfs > $multifasta.orfs");
}

sub transeq {
	 my ($multifasta1) = @_;
	
	execute("transeq -clean -frame 6 -trim -sequence $multifasta1 -outseq $multifasta1.translated"); 


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
	# don't worry I catch it here 
	# the problem with this is that it will print 10 chunk of 8 seqs, but miss one seq

	# hey this might bug if there more cpus than candidate, then the ratio will be less than 1	
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
        # don't worry I catch it here
        # the problem with this is that it will print 10 chunk of 8 seqs, but miss one seq

        # hey this might bug if there more cpus than candidate, then the ratio will be less than 1
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

        # now extract fasta seq for these chunk
        # and launch exonerate
        my $ext = "";
        foreach $ext (@chunksForFastaExtr){
                # You need to update the path later
            #    push(@runs_fas,"perl $srcdir/retrieveFasta.pl $outdir/$ext $outdir/$candidateFasta > $outdir/$ext.fas");
		push(@run_aug,"$augdir/bin/augustus --species=$speciestag --AUGUSTUS_CONFIG_PATH=$augdir/config $outdir/$ext.fas > $outdir/$ext.gff");
       		push(@run_gff2aa,"$augdir/scripts/getAnnoFasta.pl $outdir/$ext.gff");
		push(@toconcat,"$outdir/$ext.aa");
	 }
        #return($numOfseqs,$numOfChunks,$numOfseqPerChunks,\@runs_fas,\@run_aug,\@run_gff2aa);
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
		warn "Command: $cmd started at ".localtime."\n";
		my $pid = launch($cmd);
		$children{$pid} = $cmd;
	}
	
	while (%children){
		my $pid = wait();
		die $! if $pid < 1;
		my $cmd = delete($children{$pid});
		warn "Command: $cmd ended at ".localtime."with \$? = $?.\n";
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
                	warn"...only $new markers found - thresold not satisfied - attempt # $trials / 20\n";
        	} else {
                	warn"...new markers found:\t$new ( previous $previous)";
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
                $buffer .="$ids[$id]\n";
        }
        io("$iteration.tmp")->write($buffer);
        execute("perl $dir/src/retrieveFasta.pl $iteration.tmp $readsFile > $readsFile.sampled.$iteration.fa");
        execute("rm $iteration.tmp");
}
sub reservoir_sampling {
        my (@list) = @_;
        
	my %samples = ();

        my $num_of_samples = '1000';
        my $sample_size = '1000';
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
1;


