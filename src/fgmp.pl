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
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/lib";

use Getopt::Long qw(:config no_ignore_case bundling);
use FGMP;
use IPC::Cmd qw[can_run run run_forked];
use Config;
use Cwd;
use Time::Piece;

# ----------------------------------- #
# 	VARIABLES
# ----------------------------------- #
# main variables

my $SOFTWARE = "fgmp";
my $VERSION  = "1.0.2";
my @clean    = ();

my $FGMP = cwd();

# ----------------------------------- #
# 	READING OPTIONS
# ----------------------------------- #
# print help and exit if no arguments is given
&show_help unless @ARGV;

# getOptions variables
my (
        $genome,       $protein,      $output,      $blastdb,
        $hmm_profiles, $hmm_prefix,   $cutoff_file, $mark_file,
        $fuces_hmm,    $fuces_prefix, $tag,         $reads,
        $multicopies,  $verbose_str,  $verbose_flg, $quiet_flg,
        $temp_flg,     $help_flg,     $threads,     $augTraingCutoff,
        $nsamples,     $nsampleSize
  )
  = (
        undef, undef, undef, undef, undef, undef, undef, undef,
        undef, undef, undef, undef, undef, undef, 0,     0,
        0,     0,     4,     50,    1000,  10000
  );

# Reading options
&which_Options();

if ( !($genome) && ($reads) ) {
        my $makersInReads =
          FGMP::search_in_reads( $reads, $protein, $FGMP, $threads, $nsamples,
                $nsampleSize );
        $makersInReads = 'NA' unless ( defined($makersInReads) );
        my $buf .=
          "# no. of markers detected:\t$makersInReads\tof\t593 markers\n";
        io("$reads.SEARCH_IN_READS.report")->write($buf);
        &die(
"#\tFGMP will only report search in reads (see $reads.SEARCH_IN_READS.report) -- only reads provided, no genome!"
        );
}

# Reading options
&which_Options();

# ----------------------------------- #
# 	CHECKS
# ----------------------------------- #
# Retrieve platform info
my $platform = $Config{osname};

# Check if the software are installed
foreach my $soft (qw ( makeblastdb tblastn exonerate hmmsearch sixpack csplit ))
{
        my $full_path = can_run($soft) || croak "$soft is not installed\n";
}

# Check that the no. of cpus requested are available
if ( $platform =~ m/Linux/i ) {
        my $nb_cpus_on_system = `nproc`;
        chomp $nb_cpus_on_system;
        croak
"ERROR:\tNB OF CPUS REQUESTED: $threads is superior to nb. of cpus available on this system ($nb_cpus_on_system)"
          unless ( $threads <= $nb_cpus_on_system );
}
elsif ( $platform =~ m/Darwin/ ) {
        my $nb_cpus_on_system = `sysctl -n hw.ncpu1`;
        chomp $nb_cpus_on_system;
        croak
"ERROR:\tNB OF CPUS REQUESTED: $threads is superior to nb. of cpus available on this system ($nb_cpus_on_system)"
          unless ( $threads <= $nb_cpus_on_system );
}

# Decompress hmms if needed
if ( -e "$FGMP/data/593_cleanMarkers.hmm.gz" ) {
        FGMP::execute("gunzip $FGMP/data/593_cleanMarkers.hmm.gz");
        &report(
"found compressed file: $FGMP/593_cleanMarkers.hmm.gz - decompressing ..."
        );
}
else {
        &report(
"file: $FGMP/data/593_cleanMarkers.hmm.gz is already decompressed - good!"
        );
}

# ----------------------------------- #
# FIND CANDIDATE REGIONS IN THE GENOME
# ----------------------------------- #
&report("Finding candidate regions");
unless ( -e "$genome.candidates.fa" ) {
        &run_find_candidate_regions( $genome, $protein, $threads );
}

# run sixpack to translate from candidate regions, because exonerate misses some regions
if ( defined($threads) && ( $threads >= 2 ) ) {
        FGMP::split_and_run_sixpack("$genome.candidates.fa");
}

# ----------------------------------- #
# REFINING MAPPING ON SELECTED REGIONS
# ----------------------------------- #
&report("Refining mapping on candidate regions");

# initial mapping
unless ( -e "$genome.candidates.fa.p2g" ) {

        # Identify candidate regions and run exonerate in parallel
        if ( defined($threads) && ( $threads >= 2 ) ) {
                my (
                        $nb_seqs,   $nb_chunk, $nb_seq_per_chunk,
                        $fastaJobs, $exonerateJobs
                  )
                  = FGMP::multithread_exonerate( "$genome.candidates.fa",
                        "$threads", "$protein", "$FGMP/src", $FGMP );
                &report("Mapping proteins to assembly chunks using Exonerate");

                # To implement, in case someone want to export commands
                #io("$genome.run_exonerate_on_chunk.sh")->write($chunkTodo);

                # run exonerate jobs on specified nodes and wait until done;
                my $status_fas = FGMP::execute_and_returnWhendone(@$fastaJobs);

# wait for fasta files to be created, $status_fas should be 0 if all runs complete successfully
                if ( $status_fas == '0' ) {
                        my $status_exo =
                          FGMP::execute_and_returnWhendone(@$exonerateJobs);

                        # now concatenate the chunkfile into one p2g file
                        if ( $status_exo == '0' ) {
                                FGMP::execute(
"cat $FGMP/$genome*.chunk*.p2g > $FGMP/$genome.candidates.fa.p2g"
                                );

                                #FGMP::execute("rm $FGMP/$genome*.chunk*");
                        }
                        else {
                                &report(
"Potential error:\tFailed to concatenate p2g files"
                                );
                        }
                }
                else {
                        &report(
"Potential error: something went wrong with the fasta files"
                        );
                }
        }
        else {
                FGMP::run_on_single_node( "$genome.candidates.fa", $protein,
                        $FGMP );
        }
}

# recovering exonerate translated matches
unless ( -e "$FGMP/$genome.candidates.fa.p2g.aa" ) {
        FGMP::execute(
"cat $FGMP/$genome.candidates.fa.p2g | grep -v '^#' | grep -v \'exonerate:protein2genome:local\' > $FGMP/$genome.candidates.fa.withoutGFF.p2g"
        );
        FGMP::execute(
"perl -pE 's/XXX/$tag/' $FGMP/src/recoverCDS.sh > $FGMP/src/recoverCDS.$tag.sh"
        );
        FGMP::execute(
"bash $FGMP/src/recoverCDS.$tag.sh $FGMP/$genome.candidates.fa.withoutGFF.p2g"
        );
}

# check if the exonerate file is empty, because even empty exonerate generates a minimal output
my $fileSize = -s ("$FGMP/$genome.candidates.fa.withoutGFF.p2g.aa");
if (       ( $fileSize > '373' )
        && ( -e "$FGMP/$genome.candidates.fa.withoutGFF.p2g.aa" ) )
{    # the size of null report from exonerate
        FGMP::execute(
"perl $FGMP/src/exonerate2proteins.pl $FGMP/$genome.candidates.fa.withoutGFF.p2g.aa > $FGMP/$genome.candidates.fa.withoutGFF.p2g.aa.proteins"
        );
}
else {
        report("Exonerate failed");
}
push( @clean,
        "$FGMP/$genome.candidates.fa",
        "$FGMP/$genome.candidates.fa.withoutGFF.p2g",
        "$FGMP/$genome.candidates.fa.withoutGFF.p2g.aa",
        "$FGMP/$genome.candidates.fa.withoutGFF.p2g.aa.proteins" );

# ----------------------------------- #
# Ab initio predictions
# ----------------------------------- #
&report("Training Augustus");

FGMP::execute(
"exonerate2hints.pl --in=$FGMP/$genome.candidates.fa.p2g --out=$FGMP/$genome.trainingSet"
);
FGMP::execute(
"gff2gbSmallDNA.pl $FGMP/$genome.trainingSet $FGMP/$genome 100 $FGMP/$genome.trainingSet.gb &> /dev/null"
);

my $numOfGenesIngb = FGMP::how_many_locus("$genome.trainingSet.gb");

if ( $numOfGenesIngb >= $augTraingCutoff ) {
        FGMP::execute("randomSplit.pl $genome.trainingSet.gb $augTraingCutoff");

        # check if this dir already exists then erase if there

        FGMP::execute("rm -rf $FGMP/augustus_tmp")
          if ( -e "$FGMP/augustus_tmp" );
        FGMP::execute("mkdir -p $FGMP/augustus_tmp");
        FGMP::execute("cp -R $FGMP/utils/config $FGMP/augustus_tmp/");

        FGMP::execute(
"new_species.pl --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config > /dev/null 2>&1"
        );

        # training
        FGMP::execute(
"etraining --species=$genome $genome.trainingSet.gb.train --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config > /dev/null 2>&1"
        );

        if ( defined($threads) && ( $threads >= 2 ) ) {
                my ( $nb_seqsAug, $nb_chunkAug, $nb_seq_per_chunkAug,
                        $augustusJobs, $gff2aaJobs, $concatElements )
                  = FGMP::multithread_augustus( "$genome.candidates.fa",
                        "$threads", "$protein", "$FGMP/src",
                        "$FGMP/augustus_tmp", $FGMP, $genome );
                &report("Running Augustus");
                my $status_aug =
                  FGMP::execute_and_returnWhendone(@$augustusJobs);

# wait for fasta files to be created , $status_fas should be 0 if all runs complete
                if ( $status_aug == '0' ) {
                        my $statuts_gff2aa =
                          FGMP::execute_and_returnWhendone(@$gff2aaJobs);
                        if ( $statuts_gff2aa == '0' ) {

                                # clean before
                                FGMP::execute(
"echo \"\" > $FGMP/$genome.candidates.fa.aa"
                                );
                                my $element_aa = "";
                                foreach $element_aa (@$concatElements) {
                                        FGMP::execute(
"cat $element_aa >> $FGMP/$genome.candidates.fa.aa"
                                        ) if ( -e "$element_aa" );
                                        push( @clean, $element_aa );
                                }
                        }
                        else {
                                &report(
"Potential error: the GFF convertion to translated peptides failed"
                                );
                        }
                }
                else {
                        &report(
"Potential error: possibly ill formated fasta files"
                        );
                }
        }
        else {
                FGMP::execute(
"augustus --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/augustus_tmp/config $FGMP/$genome.candidates.fa > $FGMP/$genome.candidates.fa.gff"
                );
                FGMP::execute(
                        "getAnnoFasta.pl $FGMP/$genome.candidates.fa.gff");
        }

}
else {
        &report(
"Insufficient number of genes for augustus training :-(\nSkip Augustus step"
        );
}

# need script to evaluate the quality of the preds

push( @clean,
        "$genome.trainingSet",         "$genome.trainingSet.gb",
        "tee firsttest.$genome",       "$genome.trainingSet.gb.train",
        "$genome.trainingSet.gb.test", "$genome.gff" );

# ----------------------------------- #
# MERGING PREDS
# ----------------------------------- #
&report("Merging translated CDS and augustus predictions");

FGMP::execute(
"cat $FGMP/$genome.candidates.fa.withoutGFF.p2g.aa.proteins > $FGMP/$genome.unfiltered"
) if ( -e "$FGMP/$genome.candidates.fa.withoutGFF.p2g.aa.proteins" );
FGMP::execute("cat $FGMP/$genome.candidates.fa.aa >> $FGMP/$genome.unfiltered")
  if ( -e "$FGMP/$genome.candidates.fa.aa" );
FGMP::execute(
"perl $FGMP/src/rename.pl $FGMP/$genome.unfiltered > $FGMP/$genome.unfiltered.renamed"
);
my $countUn = FGMP::count_num_of_seqs("$FGMP/$genome.unfiltered.renamed");
&report("$countUn predictions to compare to HMM models");

push( @clean,
        "$FGMP/$genome.candidates.fa.withoutGFF.p2g.aa.proteins",
        "$FGMP/$genome.aa",
        "$FGMP/$genome.unfiltered" );

# ----------------------------------- #
# FILTERING
# ----------------------------------- #
&report("Filtering gene models");
if ( -s "$genome.preds.filtered" ) {
        my $countFil = FGMP::count_num_of_seqs("$genome.unfiltered.renamed");
        warn
"A filtered file already exists $genome.preds.filtered and contains $countFil\n";
}
else {
        # Protein makers
        FGMP::execute(
"hmmsearch --cpu $threads --domtblout $FGMP/$genome.unfiltered.renamed.hmmsearch $hmm_profiles $FGMP/$genome.unfiltered.renamed > $FGMP/$genome.unfiltered.renamed.hmmsearch.log"
        );

        # fUCEs
        FGMP::execute(
"nhmmer -E 1e-15 --noali  --cpu $threads --dfamtblout $FGMP/$genome.nhmmer.out $fuces_hmm $FGMP/$genome > /dev/null 2>&1"
        );

        # search in reads
        my $makersFoundInReads = "";
        if ( defined($reads) ) {
                $makersFoundInReads =
                  FGMP::search_in_reads( $reads, $protein, $FGMP, $threads );
                FGMP::execute("echo NA > makersFoundInReads.txt");
        }
        else {
                $makersFoundInReads = 'NA';
                FGMP::execute("echo NA > makersFoundInReads.txt");
        }

        # search multicopies
        FGMP::check_multicopies( "$FGMP/$genome.unfiltered.renamed.hmmsearch",
                $multicopies );

        # filtering
        FGMP::filter_unfiltByScore(
                "$FGMP/$genome.unfiltered.renamed.hmmsearch",
                $cutoff_file,
                $mark_file,
                $tag,
                "$FGMP/$genome.nhmmer.out",
                $fuces_prefix,
                "makersFoundInReads.txt",
"$FGMP/$genome.unfiltered.renamed.hmmsearch.multicopies_check.csv"
        );

        # extract fasta
        FGMP::execute(
"grep \'\^Seq\' $FGMP/$genome.unfiltered.renamed.hmmsearch.full_report | cut -f 1 > $FGMP/$genome.unfiltered.renamed.hmmsearch.full_report.tmp"
        );
        FGMP::retrieveFasta(
                "$FGMP/$genome.unfiltered.renamed.hmmsearch.full_report.tmp",
                "$FGMP/$genome.unfiltered.renamed",
                "$FGMP/$genome.bestPreds.fas"
        );
}
push( @clean,
        "$genome.unfiltered.renamed",
        "$genome.unfiltered.renamed.hmmsearch",
        "$genome.unfiltered.renamed.hmmsearch.log",
        "$genome.unfiltered.renamed.hmmsearch.full_report.tmp",
        "$genome.db.nhr",
        "$genome.db.nin",
        "$genome.db.nsq",
        "$genome.db.nin,makersFoundInReads.txt" );

# ----------------------------------- #
# CLEANING
# ----------------------------------- #
FGMP::clean_files( \@clean, \$FGMP, \$temp_flg );

# ----------------------------------- #
#       SUBS
# ----------------------------------- #
sub run_find_candidate_regions {
        my ( $genome_file, $prot, $cps ) = @_;

        # To implement: check that fasta header are properly formatted

        # run BLAST
        FGMP::execute(
"makeblastdb -in $genome_file -dbtype nucl -out $genome_file.db  > /dev/null 2>&1"
        );
        FGMP::execute(
"tblastn -db $genome_file.db -query $prot -word_size 5 -max_target_seqs 5 -evalue 0.01 -seg yes -num_threads $cps -outfmt  \"7 sseqid sstart send sframe bitscore qseqid\" > $genome_file.tblastn"
        );
        my (%adjusted) = FGMP::extractCandidateRegion("$genome_file.tblastn");
        FGMP::exportCandidateRegions( \%adjusted, $genome_file );

}

sub report {
        my $now = localtime->strftime("%d-%m-%Y %I:%M:%S %p");
        warn "[ $now ]\t@_\n";
}

sub which_Options {
        GetOptions(
                "g|genome=s"  => \$genome,     # genome fasta file
                "p|protein=s" => \$protein,    # proteins
                "o|output=s"  => \$output,     # outputfile name
                "d|blastdb=s" => \$blastdb,
                "hmm_profiles=s" =>
                  \$hmm_profiles,    # hmm profiles for filtering predictions
                "hmm_prefix=s"        => \$hmm_prefix,
                "c|cutoff_file=s"     => \$cutoff_file,
                "m|mark_file=s"       => \$mark_file,
                "fuces_hmm=s"         => \$fuces_hmm,       # fUCEs hmm
                "fuces_prefix=s"      => \$fuces_prefix,    #fUCEs names
                "r|reads=s"           => \$reads,
                "mu|multicopies=s"    => \$multicopies,
                "t|tag=s"             => \$tag,
                "v|verbose"           => \$verbose_flg,     # verbose
                "q|quiet"             => \$quiet_flg,       # quiet mode
                "tmp"                 => \$temp_flg,        # keep tempory files
                "h|help|?"            => \$help_flg,        # print help
                "T|threads=i"         => \$threads,         # number of threads
                "A|augTraingCutoff=i" => \$augTraingCutoff,
                "nsamples=s"          => \$nsamples,
                "nsampleSize=s"       => \$nsampleSize,

        ) || &show_help();
        &show_help if $help_flg;

        # check if files exists
        &die("FATAL ERROR!!! --genome or --blastdb not specified")
          unless ( defined($genome) || defined($blastdb) || defined($reads) );

        # need to check here that the $genome is in the correct fasta format
        # &die("FATAL ERROR!!! --verbose and --quiet are mutually exclusive")
        #	if ($verb_flg && $quiet_flg);

        # default settings
        $protein = "$FGMP/data/593_cleanMarkers.fa"
          if ( !( defined($protein) ) );
        $output       = "output" if ( !( defined($output) ) );
        $hmm_profiles = "$FGMP/data/593_cleanMarkers.hmm"
          if ( !( defined($hmm_profiles) ) );
        $hmm_prefix  = "OMA" if ( !( defined($hmm_prefix) ) );
        $cutoff_file = "$FGMP/data/profiles_cutOff.tbl"
          if ( !( defined($cutoff_file) ) );
        $mark_file = "$FGMP/data/593_cleanMarkers.txt"
          if ( !( defined($mark_file) ) );
        $fuces_hmm    = "$FGMP/data/all.hmm" if ( !( defined($fuces_hmm) ) );
        $fuces_prefix = "$FGMP/data/all.txt" if ( !( defined($fuces_prefix) ) );
        $multicopies  = "$FGMP/data/multicopies_lowerbound.csv"
          if ( !( defined($multicopies) ) );
        $tag             = "OMA"  if ( !( defined($tag) ) );
        $verbose_str     = " -v " if $verbose_flg;
        $threads         = 4      if ( !( defined($threads) ) );
        $augTraingCutoff = 50     if ( !( defined($augTraingCutoff) ) );
        $nsamples        = 1000   if ( !( defined($nsamples) ) );
        $nsampleSize     = 10000  if ( !( defined($nsampleSize) ) );
        $temp_flg        = 'FALSE';
}

sub die() {
        warn "@_\n";    # need to added cleaning stuff here
        say "\n#####\nEnding FGMP";
        exit(1);

}

sub show_help() {
        open( HELP, "| cat" );
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
        close(HELP);
        exit(1);
}
