#!/usr/bin/env perl

##########################################################################
#                                                                        #
#                               FGMP                                     #
#                                                                        #
##########################################################################
#                                                                        #
#                Fungal Gene Mapping Project	                         #
#                                                                        #
#              Copyright (C) 20013-2016      Ousmane H. Cisse            #
#                                            Jason E. Stajich            #
#                                                                        #
# This program is free software; you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation; either version 2 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# This program is distributed in the hope that it will be useful, but    #
# WITHOUT ANY WARRANTY; without even the implied warranty of             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      #
# General Public License for more details.                               #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with this program; if not, write to the Free Software            #
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                        #
##########################################################################

use warnings;
use strict;
use feature 'say'; 
use Data::Dumper; 
use Carp; 
use IO::All; 
use Getopt::Long qw(:config no_ignore_case bundling);
use Fgmp; 
use IPC::Cmd qw[can_run run run_forked];
use Config;

# ----------------------------------- #
# 	VARIABLES
# ----------------------------------- #
# main variables

my $SOFTWARE = "fgmp";
my $VERSION = "1.0.0";
my @clean = ();

# loading paths
my ($FGMP,$WRKDIR,$TMPDIR) = ("","","");
if (-e 'fgmp.config'){
	  ($FGMP,$WRKDIR,$TMPDIR) = Fgmp::load_paths('fgmp.config');
} else { 
	croak "$0 requires config file: fgmp.config\n";
}

# ----------------------------------- #
# 	READING OPTIONS
# ----------------------------------- #
# print help and exit if no arguments is given *
&show_help unless @ARGV; 

# getOptions variables
my($genome,$protein,$output,$blastdb,$hmm_profiles,$hmm_prefix,$cutoff_file,$mark_file,$fuces_hmm,$fuces_prefix,$tag,$reads,$multicopies,$verbose_str,$verbose_flg,$quiet_flg,$temp_flg,$help_flg,$threads,$augTraingCutoff,$nsamples,$nsampleSize) = (undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,0,0,0,0,4,50,1000,10000);

# Reading options
&which_Options(); 

if (!($genome) && ($reads)){
	my $makersInReads = Fgmp::search_in_reads($reads,$protein,$FGMP,$threads,$nsamples,$nsampleSize);
	$makersInReads = 'NA' unless (defined($makersInReads));
	my $buf .="# no. of markers detected:\t$makersInReads\tof\t593 markers\n";
	io("$reads.SEARCH_IN_READS.report")->write($buf);
	&die("#\tFGMP will only report search in reads (see $reads.SEARCH_IN_READS.report) -- only reads provided, no genome!");
}

# Reading options
&which_Options();

# ----------------------------------- #
# 	CHECKS
# ----------------------------------- #

# check platform
my $platform = $Config{osname};
#say "TEST\t$Config{archname}\n";

# check if the softs are installed -
foreach my $soft ( qw ( makeblastdb tblastn exonerate hmmsearch sixpack csplit )){
	my $full_path = can_run($soft) || croak "$soft is not installed\n";	
}

# check that the no. of cpus requested are available
if ($platform =~ m/Linux/i){
	my $nb_cpus_on_system = `nproc`;
	chomp $nb_cpus_on_system;
	croak "ERROR:\tNB OF CPUS REQUESTED: $threads is superior to nb. of cpus available on this system ($nb_cpus_on_system)" unless ($threads <= $nb_cpus_on_system);
} elsif ($platform =~ m/Darwin/){
	my $nb_cpus_on_system = `sysctl -n hw.ncpu1`;
	chomp $nb_cpus_on_system;
	croak "ERROR:\tNB OF CPUS REQUESTED: $threads is superior to nb. of cpus available on this system ($nb_cpus_on_system)" unless ($threads <= $nb_cpus_on_system);
}

# ----------------------------------- #
# FIND CANDIDATE REGIONS IN THE GENOME
# ----------------------------------- #
&report("INFO\tfinding candidate regions - TBLASTN");
unless (-e "$genome.candidates.fa"){
	&run_find_candidate_regions($genome,$protein,$threads);
}

# Implements sixpack to translate from candidate regions, because exonerate misses some regions
if (defined ($threads) && ($threads >= 2)){
	Fgmp::split_and_run_sixpack("$genome.candidates.fa");
}

# ----------------------------------- #
# REFINING MAPPING ON SELECTED REGIONS
# ----------------------------------- #
&report("INFO\tREFINING MAPPING ON SELECTED REGIONS"); 

# initial mapping
unless (-e "$genome.candidates.fa.p2g"){

	# chunk the candate regions and runn exonerate in parallel # for speed
	if (defined ($threads) && ($threads >= 2)){
 		my ($nb_seqs,$nb_chunk,$nb_seq_per_chunk,$fastaJobs,$exonerateJobs) = Fgmp::multithread_exonerate("$genome.candidates.fa","$threads","$protein","$FGMP/src",$TMPDIR);
	
		&report("CMD:LAUNCHING MULTI-THREAD EXONERATE\n\tNB OF CPUs: \t$threads\n\tNB SEQS TO ANALYZE: $nb_seqs\n\tNB OF CHUNKs: $nb_chunk\n\tAVE NB OF SEQS PER CHUNKS: $nb_seq_per_chunk");
		
		# In case someone wants to export commands and launch on a cluster for example
		#io("$genome.run_exonerate_on_chunk.sh")->write($chunkTodo);
	
		# run exonerate jobs on specified nodes and wait until done;
		my $status_fas = Fgmp::execute_and_returnWhendone(@$fastaJobs);
		
		# wait for fasta files to be created, $status_fas should be 0 if all runs complete successfully
		if ($status_fas == '0'){
			my $status_exo	= Fgmp::execute_and_returnWhendone(@$exonerateJobs);	
			
			# now concatenate the chunkfile into one p2g file
			if ( $status_exo == '0'){
				Fgmp::execute("cat $TMPDIR/$genome*.chunk*.p2g > $TMPDIR/$genome.candidates.fa.p2g"); 
				#Fgmp::execute("rm $TMPDIR/$genome*.chunk*");
			} else {
				# do something cannot concatenate 
			}
		} else {
			# do something : something went wrong with the fasta files
		}
	} else {
		Fgmp::run_on_single_node("$genome.candidates.fa",$protein,$WRKDIR);
	}
}

# recovering exonerate translated matches
unless (-e "$WRKDIR/$genome.candidates.fa.p2g.aa"){
	Fgmp::execute("cat $WRKDIR/$genome.candidates.fa.p2g | grep -v '^#' | grep -v \'exonerate:protein2genome:local\' > $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
	Fgmp::execute("perl -pE 's/XXX/$tag/' $FGMP/src/recoverCDS.sh > $FGMP/src/recoverCDS.$tag.sh");  
	Fgmp::execute("bash $FGMP/src/recoverCDS.$tag.sh $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
}

# check if the exonerate file is empty, because even empty exonerate generates a minimal outp
my $fileSize = -s ("$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa");
if (($fileSize > '373') && (-e "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa")){ # the size of null report from exonerate
	Fgmp::execute("perl $FGMP/src/exonerate2proteins.pl $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa > $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
} else {
	report("MSG\tEXONERATE failed");
}
push (@clean, "$WRKDIR/$genome.candidates.fa", "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g", "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa", "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
 
# ----------------------------------- #
# Ab initio predictions
# ----------------------------------- #
&report("INFO\tAB INITIO PREDS");

#Fgmp::execute("$FGMP/utils/augustus-3.0.3/scripts/exonerate2hints.pl --in=$WRKDIR/$genome.candidates.fa.p2g --out=$WRKDIR/$genome.trainingSet");
Fgmp::execute("exonerate2hints.pl --in=$WRKDIR/$genome.candidates.fa.p2g --out=$WRKDIR/$genome.trainingSet");
#Fgmp::execute("$FGMP/utils/augustus-3.0.3/scripts/gff2gbSmallDNA.pl $WRKDIR/$genome.trainingSet $WRKDIR/$genome 100 $WRKDIR/$genome.trainingSet.gb");
Fgmp::execute("gff2gbSmallDNA.pl $WRKDIR/$genome.trainingSet $WRKDIR/$genome 100 $WRKDIR/$genome.trainingSet.gb");

my $numOfGenesIngb = Fgmp::how_many_locus("$genome.trainingSet.gb");

if ($numOfGenesIngb >= $augTraingCutoff){
	#Fgmp::execute("$FGMP/utils/augustus-3.0.3/scripts/randomSplit.pl $genome.trainingSet.gb $augTraingCutoff");
	Fgmp::execute("randomSplit.pl $genome.trainingSet.gb $augTraingCutoff");
	
	# check if this dir already exists then erase if there
	
	Fgmp::execute("rm -rf $FGMP/augustus_tmp") if (-e "$FGMP/augustus_tmp");
	Fgmp::execute("mkdir -p $FGMP/augustus_tmp");
	#Fgmp::execute("rm -rf $FGMP/utils/augustus-3.0.3/config/species/$genome") if (-e "$FGMP/utils/augustus-3.0.3/config/species/$genome"); 
	#Fgmp::execute("perl $FGMP/utils/augustus-3.0.3/scripts/new_species.pl --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config > /dev/null 2>&1");
	Fgmp::execute("cp -R $FGMP/utils/config $FGMP/augustus_tmp/");

	#Fgmp::execute("rm $FGMP/augustus_tmp/config/species/$genome/$genome\_parameters.cfg") if (-e "$FGMP/augustus_tmp/config/species/$genome/$genome\_parameters.cfg");
	#Fgmp::execute("rm -rf $FGMP/augustus_tmp/config/species/$genome") if (-e "$FGMP/augustus_tmp/config/species/$genome");
	Fgmp::execute("new_species.pl --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config > /dev/null 2>&1");
	# training
 	#Fgmp::execute("$FGMP/utils/augustus-3.0.3/src/etraining --species=$genome $genome.trainingSet.gb.train --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config > /dev/null 2>&1");
 	Fgmp::execute("etraining --species=$genome $genome.trainingSet.gb.train --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config > /dev/null 2>&1");
	
	#Fgmp::execute("$FGMP/utils/augustus-3.0.3/bin/augustus --species=$genome $genome.trainingSet.gb.test --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config | tee firsttest.$genome");
	#Fgmp::execute("augustus --species=$genome $genome.trainingSet.gb.test --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config | tee firsttest.$genome"); # causes augustus to run forever, bug in augustus?

	 if (defined ($threads) && ($threads >= 2)){
		#my ($nb_seqsAug,$nb_chunkAug,$nb_seq_per_chunkAug,$augustusJobs,$gff2aaJobs,$concatElements) = Fgmp::multithread_augustus("$genome.candidates.fa","$threads","$protein","$FGMP/src","$FGMP/utils/augustus-3.0.3",$WRKDIR,$genome);
		my ($nb_seqsAug,$nb_chunkAug,$nb_seq_per_chunkAug,$augustusJobs,$gff2aaJobs,$concatElements) = Fgmp::multithread_augustus("$genome.candidates.fa","$threads","$protein","$FGMP/src","$FGMP/augustus_tmp",$WRKDIR,$genome);
		&report("CMD:LAUNCHING MULTI-THREAD AUGUSTURS\n\tNB OF CPUs: \t$threads\n\tNB SEQS TO ANALYZE: $nb_seqsAug\n\tNB OF CHUNKs: $nb_chunkAug\n\tAVE NB OF SEQS PER CHUNKS: $nb_seq_per_chunkAug");

 	 	my $status_aug = Fgmp::execute_and_returnWhendone(@$augustusJobs);
        
        	# wait for fasta files to be created , $status_fas should be 0 if all runs complete
                if ($status_aug == '0'){
				my $statuts_gff2aa = Fgmp::execute_and_returnWhendone(@$gff2aaJobs);
				if ( $statuts_gff2aa == '0'){
					
					# clean before
					Fgmp::execute("echo \"\" > $WRKDIR/$genome.candidates.fa.aa");
					my $element_aa = "";
						foreach $element_aa (@$concatElements){
                                			Fgmp::execute("cat $element_aa >> $WRKDIR/$genome.candidates.fa.aa") if (-e "$element_aa");
							push (@clean,$element_aa);	
					}
				} else {
                                # do something gff to amina acids conversation has failed
                        }
                } else {
                        # do something : something went wrong with the fasta files
                }
	 } else {
	 	#Fgmp::execute("$FGMP/utils/augustus-3.0.3/bin/augustus --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config $WRKDIR/$genome.candidates.fa > $WRKDIR/$genome.candidates.fa.gff");
	 	Fgmp::execute("augustus --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config $WRKDIR/$genome.candidates.fa > $WRKDIR/$genome.candidates.fa.gff");

		#Fgmp::execute("perl $FGMP/utils/augustus-3.0.3/scripts/getAnnoFasta.pl $WRKDIR/$genome.candidates.fa.gff");
		Fgmp::execute("getAnnoFasta.pl $WRKDIR/$genome.candidates.fa.gff");
	}
	
} else {
	warn"insufficient number of genes for augustus training :-(\n";
	warn"SKIPPING AUGUSTUS PREDICTIONS\n";	
}

	# need script to evaluate the quality of the preds

push (@clean, "$genome.trainingSet", "$genome.trainingSet.gb", "tee firsttest.$genome", "$genome.trainingSet.gb.train", "$genome.trainingSet.gb.test", "$genome.gff");

# ----------------------------------- #
# MERGING PREDS
# ----------------------------------- #
&report("INFO\tMERGING TRANSLATED CDS AND AUGUSTUS PREDICTIONS");
	
	Fgmp::execute("cat $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins > $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
	Fgmp::execute("cat $WRKDIR/$genome.candidates.fa.aa >> $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.aa");
	
 	Fgmp::execute("perl $FGMP/src/rename.pl $WRKDIR/$genome.unfiltered > $WRKDIR/$genome.unfiltered.renamed");	
	my $countUn = Fgmp::count_num_of_seqs("$WRKDIR/$genome.unfiltered.renamed");	
	warn"No. of unfiltered predictions\t$countUn\n";

#}
push(@clean, "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins", "$WRKDIR/$genome.aa","$WRKDIR/$genome.unfiltered");
 
# ----------------------------------- #
# FILTERING
# ----------------------------------- #
&report("INFO\tFILTERING PREDICTIONS");
if (-s "$genome.preds.filtered"){
	my $countFil = Fgmp::count_num_of_seqs("$genome.unfiltered.renamed");
	warn "A filtered file already exists $genome.preds.filtered and contains $countFil\n";
	} else {
		# protein makers
 		Fgmp::execute("hmmsearch --cpu $threads --domtblout $WRKDIR/$genome.unfiltered.renamed.hmmsearch $hmm_profiles $WRKDIR/$genome.unfiltered.renamed > $WRKDIR/$genome.unfiltered.renamed.hmmsearch.log");

		# fUCEs
  		Fgmp::execute("nhmmer -E 1e-15 --noali  --cpu $threads --dfamtblout $WRKDIR/$genome.nhmmer.out $fuces_hmm $WRKDIR/$genome > /dev/null 2>&1");
		
		# search in reads
		my $makersFoundInReads = "";
		if (defined($reads)){
			$makersFoundInReads = Fgmp::search_in_reads($reads,$protein,$FGMP,$threads);
			Fgmp::execute("echo NA > makersFoundInReads.txt");
		} else {
			$makersFoundInReads = 'NA';
			Fgmp::execute("echo NA > makersFoundInReads.txt");
		}	

		# search multicopies
		Fgmp::check_multicopies("$WRKDIR/$genome.unfiltered.renamed.hmmsearch",$multicopies);
		
		# filtering	
		#Fgmp::execute("perl $FGMP/src/filter_unfiltByScore.pl $WRKDIR/$genome.unfiltered.renamed.hmmsearch $cutoff_file $mark_file $tag $WRKDIR/$genome.nhmmer.out $fuces_prefix makersFoundInReads.txt $WRKDIR/$genome.unfiltered.renamed.hmmsearch.multicopies_check.csv --cutoff 0.7"); 
		Fgmp::filter_unfiltByScore("$WRKDIR/$genome.unfiltered.renamed.hmmsearch",$cutoff_file,$mark_file,$tag,"$WRKDIR/$genome.nhmmer.out",$fuces_prefix,"makersFoundInReads.txt","$WRKDIR/$genome.unfiltered.renamed.hmmsearch.multicopies_check.csv");
	# extract fasta
	Fgmp::execute("grep \'\^Seq\' $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report | cut -f 1 > $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp");
	#Fgmp::execute("$FGMP/src/retrieveFasta.pl $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp $WRKDIR/$genome.unfiltered.renamed > $WRKDIR/$genome.bestPreds.fas");
	Fgmp::retrieveFasta("$WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp","$WRKDIR/$genome.unfiltered.renamed","$WRKDIR/$genome.bestPreds.fas");
}
	
push(@clean,"$genome.unfiltered.renamed","$genome.unfiltered.renamed.hmmsearch","$genome.unfiltered.renamed.hmmsearch.log","$genome.unfiltered.renamed.hmmsearch.full_report.tmp","$genome.db.nhr","$genome.db.nin","$genome.db.nsq","$genome.db.nin,makersFoundInReads.txt");

# ----------------------------------- #
# CLEANING
# ----------------------------------- #
Fgmp::clean_files(\@clean,\$TMPDIR,\$temp_flg);

# ----------------------------------- #
#       SUBS
# ----------------------------------- #
sub run_find_candidate_regions {
	my ($genome_file, $prot,$cps) = @_; 
	
	# TODO : check that fasta header are properly formatted

	# run BLAST
   	Fgmp::execute("makeblastdb -in $genome_file -dbtype nucl -out $genome_file.db  > /dev/null 2>&1");
    	Fgmp::execute("tblastn -db $genome_file.db -query $prot -word_size 5 -max_target_seqs 5 -evalue 0.01 -seg yes -num_threads $cps -outfmt  \"7 sseqid sstart send sframe bitscore qseqid\" > $genome_file.tblastn");
	#Fgmp::execute("grep -v \'#\' $genome_file.tblastnOut | cut -f 1 | sort -u > $genome_file.candidates");
	my (%adjusted) = Fgmp::extractCandidateRegion("$genome_file.tblastn");	
	Fgmp::exportCandidateRegions(\%adjusted,$genome_file);

}

sub report {
	warn "###\t@_\n";
}

sub which_Options {
	GetOptions(
		"g|genome=s"		=> \$genome, 	# genome fasta file
		"p|protein=s"		=> \$protein,	# proteins
		"o|output=s"		=> \$output,	# outputfile name
		"d|blastdb=s"		=> \$blastdb,
		"hmm_profiles=s" 	=> \$hmm_profiles, # hmm profiles for filtering predictions
		"hmm_prefix=s"		=> \$hmm_prefix,
		"c|cutoff_file=s"	=> \$cutoff_file,
		"m|mark_file=s"		=> \$mark_file,
		"fuces_hmm=s"		=> \$fuces_hmm, # fUCEs hmm
		"fuces_prefix=s"	=> \$fuces_prefix, #fUCEs names
		"r|reads=s"		=> \$reads, 
		"mu|multicopies=s"	=> \$multicopies,
		"t|tag=s"		=> \$tag,
		"v|verbose"		=> \$verbose_flg,	# verbose
		"q|quiet"		=> \$quiet_flg,	# quiet mode
		"tmp"			=> \$temp_flg,	# keep tempory files
		"h|help|?"		=> \$help_flg,	# print help
		"T|threads=i"		=> \$threads, 	# number of threads 
		"A|augTraingCutoff=i"	=> \$augTraingCutoff,
		"nsamples=s"		=> \$nsamples,
		"nsampleSize=s"		=> \$nsampleSize,

	) || &show_help(); 
	&show_help if $help_flg;
	
	
	# check if files exists
	&die("FATAL ERROR!!! --genome or --blastdb not specified")
		unless (defined($genome) || defined($blastdb) || defined($reads));

	# need to check here that the $genome is in the correct fasta format
	# &die("FATAL ERROR!!! --verbose and --quiet are mutually exclusive")
	#	if ($verb_flg && $quiet_flg);
	

	# default settings
	$protein        = "$FGMP/data/593_cleanMarkers.fa" if (!(defined($protein)));
	$output 	= "output" if (!(defined($output))); 
	$hmm_profiles	= "$FGMP/data/593_cleanMarkers.hmm" if (!(defined($hmm_profiles)));
	$hmm_prefix	= "OMA" if (!(defined($hmm_prefix)));
	$cutoff_file	= "$FGMP/data/profiles_cutOff.tbl" if (!(defined($cutoff_file)));
	$mark_file	= "$FGMP/data/593_cleanMarkers.txt" if (!(defined($mark_file)));
	$fuces_hmm	= "$FGMP/data/all.hmm" if (!(defined($fuces_hmm)));
	$fuces_prefix   = "$FGMP/data/all.txt" if (!(defined($fuces_prefix))); 
	$multicopies	= "$FGMP/data/multicopies_lowerbound.csv" if (!(defined($multicopies)));
	$tag		= "OMA" if (!(defined($tag)));
	$verbose_str	= " -v " if $verbose_flg;
	$threads	= 4 if (!(defined($threads)));
	$augTraingCutoff = 50 if (!(defined($augTraingCutoff)));
	$nsamples	= 1000 if (!(defined($nsamples)));
	$nsampleSize	= 10000 if (!(defined($nsampleSize)));
 	$temp_flg	= 'FALSE';	
}
sub die() {
	warn"@_\n"; # need to added cleaning stuff here
	say "\n#####\nEnding FGMP";
	exit(1);

}
sub show_help(){
	open(HELP, "| cat");
	print HELP <<"+++EndOfHelp+++";

			$SOFTWARE

SOFTWARE:

		$SOFTWARE - $VERSION

USAGE

	$SOFTWARE [options] -g < genome_fasta_file >

DESCRIPTION


REQUIRES
	$SOFTWARE requieres the installations of the following softwares
	
	- hmmer (HMMER 3.0) 
	- NCBI blast+
	- Exonerate
	- BioPerl xxx
	- IO::All
	- Emboss sixpack & csplit


ENVIRONMENT VARIABLES
	You can specific the path where the $SOFTWARE can find the default files
	with the shell variable \"$SOFTWARE\".


	o Using a Bourne-SHell
		export FGMP="path"
		export FGMPTMP="path"
		export PERL5LIB="\$PERL5LIB:\$FGMP/lib" 

	Check exonrate cmmd

COMMAND-LINE OPTIONS

		Available options and a short description are listed here; 

	-g, --genome		genome in fasta format
	
	-p, --protein		protein seeds

	-o, --output		output file prefix

	-d, --blastdb		blast database for the genome sequence

	-c, --cutoff_file	profiles cutoff file

	-m, --mark_file		completeness markers

	-r, --reads		reads

	--fuces_hmm		fungal Ultra Conserved Elements (hmms)
	
	--fuces_prefix		 fungal Ultra Conserved Elements (names - one per line please!)

	--multicopies		default: multicopy genes from 1FKG data
	
	-t, --tag		tag to use OMA for fgmp, FUNY (Funybase) or CEG (cegma)

	-T, --threads		Specify the number of processor threads to use

	-v, --verbose		show progress

	-q, quiet		suppress show log

	-h, --help		show this help

	--tmp			keep temporary files

	-augTraingCutoff	specify the num of genes for augustus training

	--nsampleSize		Specify the # of samples for "SEARCH IN READS" module

	--nsampleSize		Specify the sample size for the "SEARCH IN READS" module

BUGS:
	Please report bugs to 'ousmanecis\@gmail.com'.

AUTHORS:

	$SOFTWARE has been developped by Ousmane H. Cisse and Jason E. Stajich.



GNU-GPL (C) 		date				 $SOFTWARE
+++EndOfHelp+++
	close (HELP);
	exit(1);
}
