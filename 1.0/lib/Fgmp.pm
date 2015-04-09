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
	my ($fileP) = @_; 

	my @settings = ();
	my $paths = io("$fileP");
	   $paths->autoclose(0);
	   while (my $lp = $paths->getline || $paths->getline){
	   chomp $lp;
		next if $lp =~ m/^#/;
		my @datapaths = split/=/, $lp;
		push(@settings, $datapaths[1]);
	}
	return(@settings[0..4]);
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
		execute("tar -cvzf $tag-temporyfiles.tar.gz $directory/$tag-temporyfiles");
		execute("rm -rf $directory/$tag-temporyfiles"); 
	}
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
		my @datatb = split /\t/, $tbline; 
		my ($target,$sstart,$send) = ($datatb[1],$datatb[8],$datatb[9]); 
		
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
		if (min(@{$h1{$el}}) < 5000){
			$startBoundaries = 0;
		} else {
			$startBoundaries = (min(@{$h1{$el}}) - 5000);
		}
		
		# you might need the size of the contig
		my $endBoundaries = (max(@{$h1{$el}}) + 5000);

		# populating the hash
		push(@bounds,$startBoundaries,$endBoundaries);		 
		@{$boundaries{$el}} = @bounds;

		# clean @bounds 
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
1;


