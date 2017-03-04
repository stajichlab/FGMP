#!/usr/bin/env perl

##########################################################################
#                                                                        #
#                               FGMP                                     #
#                                                                        #
##########################################################################
#                                                                        #
#                Fungal Genome Mapping Project	                         #
#                                                                        #
#              Copyright (C) 20013-2017      Ousmane H. Cisse            #
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

use strict;
use warnings;
use feature 'say'; 
use Data::Dumper; 
use Carp; 
use IO::All; 
use Getopt::Long qw(:config no_ignore_case bundling);
use FGMP; 
use IPC::Cmd qw[can_run run run_forked];

# ----------------------------------- #
# 	VARIABLES
# ----------------------------------- #
# main variables

my $SOFTWARE = "fgmp";
my $VERSION = "1.0";
my @clean = ();

# loading paths
my ($FGMP,$WRKDIR,$TMP) = ("","","");
if (-e 'fgmp.config'){
	  ($FGMP,$WRKDIR,$TMP) = Fgmp::load_paths('fgmp.config');
} else { 
	croak "$0 requires config file: fgmp.config\n";
}

# ----------------------------------- #
# 	READING OPTIONS
# ----------------------------------- #
# print help and exit if no arguments is given *
&show_help unless @ARGV; 

# getOptions variables

my ($genome,$protein,$output,$blastdb,$hmm_profiles,
$hmm_prefix,$cutoff_file,$mark_file,$fuces_hmm,
$fuces_prefix,$tag,$reads,$verbose_str,$verbose_flg,
$quiet_flg,$temp_flg,$help_flg,$threads,
$augTraingCutoff);

my ($threads,$augTrainingCutff) =  (4,50);

# Reading options
&which_Options(); 

if ( ! $genome && $reads){
    my $makersInReads = FGMP::search_in_reads($reads,$protein,
					      $FGMP,$threads);
    my $buf .= "# no. of reads detected:\t$makersInReads\n";
    io("$reads.SEARCH_IN_READS.report")->write($buf);
    &fgmp_die("#\tstop here -- only reads provided, no genome!");
}

# Reading options
&which_Options();


# ----------------------------------- #
# 	CHECKS
# ----------------------------------- #
# check if the softs are installed - TODO: catch nicely the die message
# AUDIT This should be replaced with path of these tool set in config file
foreach my $soft ( qw ( makeblastdb tblastn exonerate hmmsearch sixpack csplit )){
    my $full_path = can_run($soft) || croak "$soft is not installed\n";	
}

# check that the no. of cpus requested are available
# AUDIT this is fragile code- not all systems have nproc - consider module or better
# AUDIT bullet-proofing of this code
my $nb_cpus_on_system = `nproc`;
chomp $nb_cpus_on_system;

croak "ERROR:\tNUM OF CPUS REQUESTED: More threads ($threads) than available CPUS ($nb_cpus_on_system)" 
    if ($threads > $nb_cpus_on_system);

# ----------------------------------- #
# FIND CANDIDATE REGIONS IN THE GENOME
# ----------------------------------- #
&report("INFO\tfinding candidate regions - TBLASTN");
unless (-e "$genome.candidates.fa"){
	&run_find_candidate_regions("$WRKDIR/$genome",$protein,$threads);
}

# Implements sixpack to translate from candidate regions, because exonerate miss some regions
if (defined ($threads) && ($threads >= 2)){
	FGMP::split_and_run_sixpack("$genome.candidates.fa");
}

# ----------------------------------- #
# REFINING MAPPING ON SELECTED REGIONS
# ----------------------------------- #
&report("INFO\tREFINING MAPPING ON SELECTED REGIONS"); 

# initial mapping
unless (-e "$genome.candidates.fa.p2g") {

    # chunk the candate regions and runn exonerate in parallel # for speed
    if (defined ($threads) && ($threads >= 2)){
	my ($nb_seqs,$nb_chunk,$nb_seq_per_chunk,$fastaJobs,
	    $exonerateJobs) = FGMP::multithread_exonerate("$genome.candidates.fa","$threads","$protein","$FGMP/src",$WRKDIR);

	&report("CMD:LAUNCHING MULTI-THREAD EXONERATE\n\tNB OF CPUs: \t$threads\n\tNB SEQS TO ANALYZE: $nb_seqs\n\tNB OF CHUNKs: $nb_chunk\n\tAVE NB OF SEQS PER CHUNKS: $nb_seq_per_chunk");

	# In case someone want to export commands and launch on a cluster for example
	#io("$genome.run_exonerate_on_chunk.sh")->write($chunkTodo);;

	# run exonerate jobs on specified nodes and wait until done;
	my $status_fas = FGMP::execute_and_returnWhendone(@$fastaJobs);

	# wait for fasta files to be created, $status_fas should be 0 if all runs complete successfully
	if ($status_fas == '0'){
	    my $status_exo = FGMP::execute_and_returnWhendone(@$exonerateJobs);	
	    
	    # now concatenate the chunkfile into one p2g file
	    if ( $status_exo == '0'){
		FGMP::execute("cat $WRKDIR/$genome*.chunk*.p2g > $WRKDIR/$genome.candidates.fa.p2g"); 
		FGMP::execute("rm $genome*.chunk*");
	    } else {
				# do something cannot conca 
	    }
	} else {
	    # do something : something went wrong with the fasta files
	}
    } else {
	FGMP::run_on_single_node("$genome.candidates.fa",$protein,$WRKDIR);
    }
}

# recovering exonerate translated matches
unless (-e "$WRKDIR/$genome.candidates.fa.p2g.aa") {
    FGMP::execute("cat $WRKDIR/$genome.candidates.fa.p2g | grep -v '^#' | grep -v \'exonerate:protein2genome:local\' > $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
#	Fgmp::execute("$FGMP/src/recoverCDS.sh $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
    FGMP::execute("perl -pE 's/XXX/$tag/' $FGMP/src/recoverCDS.sh > $FGMP/src/recoverCDS.$tag.sh");  
    FGMP::execute("bash $FGMP/src/recoverCDS.$tag.sh $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
}

# check if the exonerate file is empty, because even empty exonerate generates a minimal output
# AUDIT make this a stored variable ?
my $fileSize = -s ("$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa");
if ( ($fileSize > '373') # the size of null report from exonerate
     && (-e "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa") ) { 
   #AUDIT can we do this another way? make it a routine?
    FGMP::execute("perl $FGMP/src/exonerate2proteins.pl $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa > $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
} else {
    report("MSG\tEXONERATE failed");
}
# AUDIT - consider using File::Spec->cafile() instead for making path
push (@clean, "$WRKDIR/$genome.candidates.fa", 
      "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g", 
      "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa", 
      "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
 
# ----------------------------------- #
## abinitio predictions
# ----------------------------------- #
&report("INFO\tAB INITIO PREDS");

#AUDIT consider parameterizing the exe path as variable instead of putting
#      it in here. also makes new version install/upgrade easier
FGMP::execute("perl $FGMP/utils/augustus-3.0.3/scripts/exonerate2hints.pl --in=$WRKDIR/$genome.candidates.fa.p2g --out=$WRKDIR/$genome.trainingSet");
FGMP::execute("perl $FGMP/utils/augustus-3.0.3/scripts/gff2gbSmallDNA.pl $WRKDIR/$genome.trainingSet $WRKDIR/$genome 100 $WRKDIR/$genome.trainingSet.gb");

#AUDIT rewrite for better english and readability - "count_genbank_loci":
my $numOfGenesIngb = FGMP::how_many_locus("$genome.trainingSet.gb");

if ($numOfGenesIngb >= $augTraingCutoff) {

#AUDIT consider parameterizing the exe path as variable instead of putting
#      it in here
    FGMP::execute("perl $FGMP/utils/augustus-3.0.3/scripts/randomSplit.pl $genome.trainingSet.gb $augTraingCutoff");

#AUDIT use perl only option instead of relying on system rm call?
    # check if this dir already exists then erase if there
    FGMP::execute("rm -rf $FGMP/utils/augustus-3.0.3/config/species/$genome") if (-e "$FGMP/utils/augustus-3.0.3/config/species/$genome"); 
#AUDIT parameterize the path with a variable?
    FGMP::execute("perl $FGMP/utils/augustus-3.0.3/scripts/new_species.pl --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config > /dev/null 2>&1");

    # training
#AUDIT parameterize the path with a variable?
    FGMP::execute("$FGMP/utils/augustus-3.0.3/src/etraining --species=$genome $genome.trainingSet.gb.train --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config > /dev/null 2>&1");
#AUDIT parameterize the path with a variable? also makes version upgrade easier
    FGMP::execute("$FGMP/utils/augustus-3.0.3/bin/augustus --species=$genome $genome.trainingSet.gb.test --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config | tee firsttest.$genome");


#AUDIT parameterize so version of augustus can be changed    
    if (defined ($threads) && ($threads >= 2)){
	my ($nb_seqsAug,$nb_chunkAug,$nb_seq_per_chunkAug,
	    $augustusJobs,$gff2aaJobs,
	    $concatElements) = FGMP::multithread_augustus
		("$genome.candidates.fa","$threads","$protein",
		 "$FGMP/src","$FGMP/utils/augustus-3.0.3",$WRKDIR,$genome);
	
	&report("CMD:LAUNCHING MULTI-THREAD AUGUSTURS\n\tNB OF CPUs: \t$threads\n\tNB SEQS TO ANALYZE: $nb_seqsAug\n\tNB OF CHUNKs: $nb_chunkAug\n\tAVE NB OF SEQS PER CHUNKS: $nb_seq_per_chunkAug");
	
	my $status_aug = FGMP::execute_and_returnWhendone(@$augustusJobs);

	# wait for fasta files to be created , $status_fas should be 0 if all runs complete
	
	if ($status_aug == 0 ) {
	    my $statuts_gff2aa = FGMP::execute_and_returnWhendone(@$gff2aaJobs);
	    if ( $statuts_gff2aa == 0){

		# clean before
		# AUDIT - this could be redone without redirect/system call
		FGMP::execute("echo \"\" > $WRKDIR/$genome.candidates.fa.aa");
		my $element_aa = "";
		#AUDIT this should be just an open file statement?
		foreach $element_aa (@$concatElements){
		    FGMP::execute("cat $element_aa >> $WRKDIR/$genome.candidates.fa.aa") if (-e "$element_aa");
		    push (@clean,$element_aa);	
		}
	    } else {
		# do soementmething gff to amino acids conversation has failed
		#AUDIT warnings needed
	    }
	} else {
	    # do something : something went wrong with the fasta files
	    #AUDIT warnings/die needed?
	}
    } else {
	# AUDIT parameterize the path to augustus so version upgrade is simpler
	# How are configs working if already been created?
	FGMP::execute("$FGMP/utils/augustus-3.0.3/bin/augustus --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config $WRKDIR/$genome.candidates.fa > $WRKDIR/$genome.candidates.fa.gff");

	# AUDIT parameterize
	FGMP::execute("perl $FGMP/utils/augustus-3.0.3/scripts/getAnnoFasta.pl $WRKDIR/$genome.candidates.fa.gff");
    }
} else {
    warn"insufficient number of genes for augustus training :-(\n";
    warn"SKIPPING AUGUSTUS PREDICTIONS\n";	
}

# need script to evaluate the quality of the preds

push (@clean, "$genome.trainingSet", 
      "$genome.trainingSet.gb", 
      "tee firsttest.$genome", #AUDIT I don't understand what this is doing here
      "$genome.trainingSet.gb.train", 
      "$genome.trainingSet.gb.test", "$genome.gff");

# ----------------------------------- #
# MERGING PREDS
# ----------------------------------- #
&report("INFO\tMERGING TRANSLATED CDS AND AUGUSTUS PREDICTIONS");
#AUDIT This can be rewritten with an open and print statements within perl
#AUDIT rather than a system call
FGMP::execute("cat $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins > $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
#AUDIT This can be rewritten with an open and print statements within perl
#AUDIT rather than a system call

FGMP::execute("cat $WRKDIR/$genome.candidates.fa.aa >> $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.aa");
#FGMP::execute("cat $WRKDIR/$genome.candidates.fa.orfs >> $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.orfs");
	
# test transeq
#AUDIT This can be rewritten with an open and print statements within perl
#AUDIT rather than a system call
FGMP::execute("cat $WRKDIR/$genome.candidates.fa.translated >> $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.translated");
	
#AUDIT This can be rewritten as a routine within the module rather
#AUDIT than an external system call of pelr
FGMP::execute("perl $FGMP/src/rename.pl $WRKDIR/$genome.unfiltered > $WRKDIR/$genome.unfiltered.renamed");

my $countUn = FGMP::count_num_of_seqs("$WRKDIR/$genome.unfiltered.renamed");	
warn"No. of unfiltered predictions\t$countUn\n";


push(@clean, "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins", 
     "$WRKDIR/$genome.aa","$WRKDIR/$genome.unfiltered");
 
# ----------------------------------- #
# FILTERING
# ----------------------------------- #
&report("INFO\tFILTERING PREDICTIONS");
if (-s "$genome.preds.filtered") {
    my $countFil = FGMP::count_num_of_seqs("$genome.unfiltered.renamed");
    warn "A filtered file already exists $genome.preds.filtered and contains $countFil\n";
} else {
    # protein makers
    FGMP::execute("hmmsearch --cpu $threads --domtblout $WRKDIR/$genome.unfiltered.renamed.hmmsearch $hmm_profiles $WRKDIR/$genome.unfiltered.renamed > $WRKDIR/$genome.unfiltered.renamed.hmmsearch.log");
    
    # fUCEs
    FGMP::execute("nhmmer -E 1e-5 --noali  --cpu $threads --dfamtblout $WRKDIR/$genome.nhmmer.out $fuces_hmm $WRKDIR/$genome > /dev/null 2>&1");
    
    # search in reads
    my $makersFoundInReads = "";
    if (defined($reads)){
	$makersFoundInReads = FGMP::search_in_reads($reads,$protein,$FGMP,$threads);
    }	
    
    # need to clean before

    
    # filtering	
# AUDIT this needs to be better parameterize and/or a subroutine?
    FGMP::execute("perl $FGMP/src/filter_unfiltByScore.pl $WRKDIR/$genome.unfiltered.renamed.hmmsearch $cutoff_file $mark_file $tag $WRKDIR/$genome.nhmmer.out $fuces_prefix $makersFoundInReads --cutoff 0.7"); 
	
    # extract fasta
    #AUDIT should this be within perl code no need to create a file this way?
    #AUDIT could combine this with step 2 below in one routine
    FGMP::execute("grep \'\^Seq\' $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report | cut -f 1 > $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp");

    #AUDIT combine with above
    FGMP::execute("$FGMP/src/retrieveFasta.pl $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp $WRKDIR/$genome.unfiltered.renamed > $WRKDIR/$genome.bestPreds.fas");
}

push(@clean,"$genome.unfiltered.renamed",
     "$genome.unfiltered.renamed.hmmsearch",
     "$genome.unfiltered.renamed.hmmsearch.log",
     "$genome.unfiltered.renamed.hmmsearch.full_report.tmp",
     "$genome.db.nhr","$genome.db.nin");

# ----------------------------------- #
# CLEANING
# ----------------------------------- #

FGMP::clean_files(@clean,$TMP,$temp_flg);

#FGMP::clean_files(@clean,"$genome",$WRKDIR);

# ----------------------------------- #
#       SUBS
# ----------------------------------- #
sub run_find_candidate_regions {
    my ($genome_file, $prot,$cpus) = @_; 

    #AUDIT - I would make variables here that define the db, output paths
    #AUDIT - and pass these around to the cmdline arguments below for simplicity
    # TODO : check that fasta header are properly formatted

    # run BLAST
    #AUDIT BLAST exes need to be set in a global config file so that exact
    #AUDIT path can be set.
    FGMP::execute("makeblastdb -in $genome_file -dbtype nucl -out $genome_file.db  > /dev/null 2>&1");
    #AUDIT BLAST exes need to be set in a global config file so that exact
    #AUDIT path can be set.

    FGMP::execute("tblastn -db $genome_file.db -query $prot -word_size 6 -max_target_seqs 5 -evalue 0.01 -seg yes -num_threads $cpus -outfmt  \"7 sseqid sstart send sframe bitscore qseqid\" > $genome_file.tblastn");

    #FGMP::execute("grep -v \'#\' $genome_file.tblastnOut | cut -f 1 | sort -u > $genome_file.candidates");

    #AUDIT use hash reference here as the returned value so that can 
    #AUDIT can be single mem ref passed around?
    #AUDIT - will you cache tblastn results to speed up on subsequent re-runs? 
    #AUDIT - if not - run it and have output parsed on the fly to avoid creating the file?
    my (%adjusted) = FGMP::extractCandidateRegion("$genome_file.tblastn");	
    
    FGMP::exportCandidateRegions(\%adjusted,$genome_file);
}

sub report {
    for my $msg ( @_ ) {
	warn "###\t$mesg\n";
    }
}


sub which_Options {
    #AUDIT - scope is a little confusing here - these are global variables
    #AUDIT - set in this function, not totally best practices.
    GetOptions(
	"g|genome=s"		=> \$genome, # genome fasta file
	"p|protein=s"		=> \$protein,	     # proteins
	"o|output=s"		=> \$output,	     # outputfile name
	"d|blastdb=s"		=> \$blastdb,
	"hmm_profiles=s" 	=> \$hmm_profiles, # hmm profiles for filtering predictions
	"hmm_prefix=s"		=> \$hmm_prefix,
	"c|cutoff_file=s"	=> \$cutoff_file,
	"m|mark_file=s"		=> \$mark_file,
	"fuces_hmm=s"		=> \$fuces_hmm,		   # fUCEs hmm
	"fuces_prefix=s"	=> \$fuces_prefix, #fUCEs names
	"r|reads=s"		=> \$reads, 
	"t|tag=s"		=> \$tag,
	"v|verbose"		=> \$verbose_flg,	# verbose
	"q|quiet"		=> \$quiet_flg,		# quiet mode
	"tmp"			=> \$temp_flg, # keep tempory files
	"h|help|?"		=> \$help_flg, # print help
	"T|threads=i"		=> \$threads,  # number of threads 
	"A|augTraingCutoff=i"	=> \$augTraingCutoff,
	) || &show_help(); 
    &show_help if $help_flg;

    # check if files exists
    &fgmp_die("FATAL ERROR!!! --genome or --blastdb not specified")
	unless (defined($genome) || defined($blastdb) || defined($reads));

    # need to check here that the $genome is in the correct fasta format
    #&fgmp_die("FATAL ERROR!!! --verbose and --quiet are mutually exclusive")
    #	if ($verb_flg && $quiet_flg);


    # default settings
    #AUDIT this needs some top level parameters
    $protein        = "$FGMP/data/593_cleanMarkers.fa" if (!(defined($protein)));
    $output 	= "output" if (!(defined($output))); 
    #AUDIT this needs some top level parameters
    $hmm_profiles	= "$FGMP/data/593_cleanMarkers.hmm" if (!(defined($hmm_profiles)));

    #AUDIT this needs some top level parameters setting for all the hmms
    #AUDIT and marker files? This seems buried, can we make it part of the
    #AUDIT config file

    $hmm_prefix	     = "OMA" unless defined $hmm_prefix;
    $cutoff_file     = "$FGMP/data/profiles_cutOff.tbl" unless defined $cutoff_file;
    $mark_file	     = "$FGMP/data/593_cleanMarkers.txt" unless defined $mark_file;
    $fuces_hmm	     = "$FGMP/data/172_fUCEs.hmm" unless defined $fuces_hmm;
    $fuces_prefix    = "$FGMP/data/172_fUCEs.txt" unless defined $fuces_prefix;
    $tag	     = "OMA" unless defined($tag);
    $verbose_str     = " -v " if $verbose_flg; #AUDIT can this be a flag? unclear why passing this way
    $threads         = 4 unless defined $threads;
    $augTraingCutoff = 50 unless defined $augTraingCutoff;
}

sub fgmp_die {
	warn "@_\n"; # need to added cleaning stuff here
	say "\n#####\nEnding FGMP";
	exit(1);
}

#AUDIT Rewrite this as POD and then show_help is going to be
#AUDIT simply a function that calls out `perldoc $0`

sub show_help() {
	open(HELP, "| cat");
	print HELP <<"+++EndOfHelp+++";

			$SOFTWARE

SOFTWARE:

		$SOFTWARE - $VERSION

		xxxxxxxxxxxxxxx

USAGE

	$SOFTWARE [options] -g < genome_fasta_file

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
	
	--fuces_prefix		fungal Ultra Conserved Elements (names - one per line please!)

	-t, --tag		tag to use OMA for fgmp, FUNY (Funybase) or CEG (cegma)

	-T, --threads		Specify the number of processor threads to use

	-v, --verbose		show progress

	-q, quiet		suppress show log

	-h, --help		show this help

	--tmp			keep temporary files

	-augTraingCutoff	specify the num of genes for augustus training

BUGS:
	Please report bugs on GitHub https://github.com/stajichlab/FGMP/issues

AUTHORS:

	$SOFTWARE has been developped by Ousmane H. Cisse and Jason E. Stajich.



GNU-GPL (C) 		date				 $SOFTWARE
+++EndOfHelp+++
	close (HELP);
	exit(1);
}
