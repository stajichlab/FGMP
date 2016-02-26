#!/opt/perl/5.16.3/bin/perl -w

use strict;
use feature 'say'; 
use Data::Dumper; 
use Carp; 
use IO::All; 
use Getopt::Long qw(:config no_ignore_case bundling);
use Fgmp; 
use IPC::Cmd qw[can_run run run_forked];

# perl fgmp_draft.pl -g sample.dna -p sample.prot --augTraingCutoff 5
# perl fgmp_draft.pl -g sample.dna -p sample.prot -T 3

# ----------------------------------- #
# 	VARIABLES
# ----------------------------------- #
# main variables

my $SOFTWARE = "fgmp";
my $VERSION = "1.0";
my @clean = ();

#my $FGMP 	= defined($ENV{'FGMP'})	? $ENV{'FGMP'} : '.'; # to implement later

#my $BIN		= defined($ENV{'FGMP'}) ? "$FGMP/bin" : '.';
#my $SRC	="/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0/src";  # should be bin and I will do when ready to compile stuff

#my $FGMPdata	= defined($ENV{'FGMP'}) ? "$FGMP/data" : "../data";

#my $TMP		= defined($ENV{'FGMP'}) ? $ENV{'FGMPTMP'} : '/tmp';

#my $TMPROOT	= "fgmp_$$";
#my $FGMPTMP	= "$TMP/$TMPROOT";

# loading paths

my ($FGMP,$WRKDIR) = "";
if (-e 'fgmp.config'){
	  ($FGMP,$WRKDIR) = Fgmp::load_paths('fgmp.config');
} else { 
	croak "$0 requires config file: fgmp.config\n";
}

# ----------------------------------- #
# 	READING OPTIONS
# ----------------------------------- #
# print help and exit if no arguments is given *
&show_help unless @ARGV; 

# getOptions variables
my($genome,$protein,$output,$blastdb,$hmm_profiles,$hmm_prefix,$cutoff_file,$mark_file,$fuces_hmm,$fuces_prefix,$tag,$reads,$verbose_str,$verbose_flg,$quiet_flg,$temp_flg,$help_flg,$threads,$augTraingCutoff) = (undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,0,0,0,0,4,50);

# Reading options
&which_Options(); 

# ----------------------------------- #
# 	CHECKS
# ----------------------------------- #
# check if the softs are installed - TODO: catch nicely the die message
foreach my $soft ( qw ( makeblastdb tblastn exonerate hmmsearch sixpack csplit )){
	my $full_path = can_run($soft) || croak "$soft is not installed\n";	
}

# check that the no. of cpus requested are available
my $nb_cpus_on_system = `nproc`;
chomp $nb_cpus_on_system;
croak "ERROR:\tNB OF CPUS REQUESTED: $threads is superior to nb. of cpus available on this system ($nb_cpus_on_system)" unless ($threads <= $nb_cpus_on_system);

# ----------------------------------- #
# FIND CANDIDATE REGIONS IN THE GENOME
# ----------------------------------- #
&report("INFO\tfinding candidate regions - TBLASTN");
unless (-e "$genome.candidates.fa"){
	&run_find_candidate_regions("$WRKDIR/$genome",$protein,$threads);
}

# test to confirm that the problem come from the candidate section -  before I fix the probleme of "candidate regions
#Fgmp::execute("cp $genome $genome.candidates.fa");

# Implements sixpack to translate from candidate regions, because exonerate miss some regions
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
 		my ($nb_seqs,$nb_chunk,$nb_seq_per_chunk,$fastaJobs,$exonerateJobs) = Fgmp::multithread_exonerate("$genome.candidates.fa","$threads","$protein","$FGMP/src",$WRKDIR);
	
		&report("CMD:LAUNCHING MULTI-THREAD EXONERATE\n\tNB OF CPUs: \t$threads\n\tNB SEQS TO ANALYZE: $nb_seqs\n\tNB OF CHUNKs: $nb_chunk\n\tAVE NB OF SEQS PER CHUNKS: $nb_seq_per_chunk");
		
		# In case we want to write somewhere command file
		#io("$genome.run_exonerate_on_chunk.sh")->write($chunkTodo);;
	
		# run exonerate jobs on specified nodes and wait until done;
		my $status_fas = Fgmp::execute_and_returnWhendone(@$fastaJobs);
		
		# wait for fasta files to be created , $status_fas should be 0 if all runs complete
		if ($status_fas == '0'){
			my $status_exo	= Fgmp::execute_and_returnWhendone(@$exonerateJobs);	
			
			# now concatenate the chunkfile into one p2g file
			if ( $status_exo == '0'){
				Fgmp::execute("cat $WRKDIR/$genome*.chunk*.p2g > $WRKDIR/$genome.candidates.fa.p2g"); 
				#Fgmp::execute("rm $genome*.chunk*");
			} else {
				# do something cannot conca 
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
#	Fgmp::execute("$FGMP/src/recoverCDS.sh $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
	Fgmp::execute("perl -pE 's/XXX/$tag/' $FGMP/src/recoverCDS.sh > $FGMP/src/recoverCDS.$tag.sh");  
	Fgmp::execute("bash $FGMP/src/recoverCDS.$tag.sh $WRKDIR/$genome.candidates.fa.withoutGFF.p2g");
}

# check if the exonerate file is empty, because even empty exonerate generates a minimal outp
my $fileSize = -s ("$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa");
if (($fileSize > '373') && (-e "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa")){ # the size of null report from exonerate
	Fgmp::execute("$FGMP/src/exonerate2proteins.pl $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa > $WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
} else {
	report("MSG\tEXONERATE has failed");
}
push (@clean, "$WRKDIR/$genome.candidates.fa", "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g", "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa", "$WRKDIR/$genome.candidates.fa.withoutGFF.p2g.aa.proteins");
 
# ----------------------------------- #
## abinitio predictions
# ----------------------------------- #
&report("INFO\tAB INITIO PREDS");

Fgmp::execute("$FGMP/utils/augustus-3.0.3/scripts/exonerate2hints.pl --in=$WRKDIR/$genome.candidates.fa.p2g --out=$WRKDIR/$genome.trainingSet");
Fgmp::execute("$FGMP/utils/augustus-3.0.3/scripts/gff2gbSmallDNA.pl $WRKDIR/$genome.trainingSet $WRKDIR/$genome 100 $WRKDIR/$genome.trainingSet.gb");

my $numOfGenesIngb = Fgmp::how_many_locus("$genome.trainingSet.gb");

if ($numOfGenesIngb >= $augTraingCutoff){
		Fgmp::execute("$FGMP/utils/augustus-3.0.3/scripts/randomSplit.pl $genome.trainingSet.gb $augTraingCutoff");
	
	# check if this dir already exists then erase if there
	Fgmp::execute("rm -rf $FGMP/utils/augustus-3.0.3/config/species/$genome") if (-e "$FGMP/utils/augustus-3.0.3/config/species/$genome"); 
	Fgmp::execute("perl $FGMP/utils/augustus-3.0.3/scripts/new_species.pl --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config > /dev/null 2>&1");

	# training
 	Fgmp::execute("$FGMP/utils/augustus-3.0.3/src/etraining --species=$genome $genome.trainingSet.gb.train --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config > /dev/null 2>&1");
 	Fgmp::execute("$FGMP/utils/augustus-3.0.3/bin/augustus --species=$genome $genome.trainingSet.gb.test --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config | tee firsttest.$genome");

	 if (defined ($threads) && ($threads >= 2)){
		my ($nb_seqsAug,$nb_chunkAug,$nb_seq_per_chunkAug,$augustusJobs,$gff2aaJobs,$concatElements) = Fgmp::multithread_augustus("$genome.candidates.fa","$threads","$protein","$FGMP/src","$FGMP/utils/augustus-3.0.3",$WRKDIR,$genome);

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
	 	Fgmp::execute("$FGMP/utils/augustus-3.0.3/bin/augustus --species=$genome --AUGUSTUS_CONFIG_PATH=$FGMP/utils/augustus-3.0.3/config $WRKDIR/$genome.candidates.fa > $WRKDIR/$genome.candidates.fa.gff");
	 	Fgmp::execute("perl $FGMP/utils/augustus-3.0.3/scripts/getAnnoFasta.pl $WRKDIR/$genome.candidates.fa.gff");
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
	#Fgmp::execute("cat $WRKDIR/$genome.candidates.fa.orfs >> $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.orfs");
	
	# test transeq
	Fgmp::execute("cat $WRKDIR/$genome.candidates.fa.translated >> $WRKDIR/$genome.unfiltered") if (-e "$WRKDIR/$genome.candidates.fa.translated");
	
 	Fgmp::execute("$FGMP/src/rename.pl $WRKDIR/$genome.unfiltered > $WRKDIR/$genome.unfiltered.renamed");	
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
 		Fgmp::execute("nhmmer -E 1e-5 --noali  --cpu $threads --dfamtblout $WRKDIR/$genome.nhmmer.out $fuces_hmm $WRKDIR/$genome > /dev/null 2>&1");
		
		# search in reads
		my $makersFoundInReads = "";
		if (defined($reads)){
			$makersFoundInReads = Fgmp::search_in_reads($reads,$protein,$FGMP,$threads);
		}	

	# need to clean before

#	Fgmp::execute("perl -pi.old -E 's/.tmp.aln.trimed//' $genome.unfiltered.renamed.hmmsearch");
#	Fgmp::execute("perl -pi.old -E 's/.aln.trim//' $genome.unfiltered.renamed.hmmsearch");

	# filtering	
	Fgmp::execute("perl $FGMP/src/filter_unfiltByScore.pl $WRKDIR/$genome.unfiltered.renamed.hmmsearch $cutoff_file $mark_file $tag $WRKDIR/$genome.nhmmer.out $fuces_prefix $makersFoundInReads --cutoff 0.7"); 
	
	# extract fasta
	Fgmp::execute("grep \'\^Seq\' $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report | cut -f 1 > $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp");
	Fgmp::execute("$FGMP/src/retrieveFasta.pl $WRKDIR/$genome.unfiltered.renamed.hmmsearch.full_report.tmp $WRKDIR/$genome.unfiltered.renamed > $WRKDIR/$genome.bestPreds.fas");
}
push(@clean,"$genome.unfiltered.renamed","$genome.unfiltered.renamed.hmmsearch","$genome.unfiltered.renamed.hmmsearch.log","$genome.unfiltered.renamed.hmmsearch.full_report.tmp","$genome.db.nhr","$genome.db.nin");





# ----------------------------------- #
# CLEANING
# ----------------------------------- #
#Fgmp::clean_files(@clean,$TMP,$temp_flg);
#Fgmp::clean_files(@clean,"$genome",$WRKDIR);



# ----------------------------------- #
#       SUBS
# ----------------------------------- #
sub run_find_candidate_regions {
	my ($genome_file, $prot,$cps) = @_; 
	
	# TODO : check that fasta header are properly formatted

	# run BLAST+ (too slow)
#	Fgmp::execute("makeblastdb -in $genome_file -dbtype nucl -parse_seqids -out $genome_file.db  > /dev/null 2>&1"); # the -parse_seqids makes the program with names with space
   	 Fgmp::execute("makeblastdb -in $genome_file -dbtype nucl -out $genome_file.db  > /dev/null 2>&1");
	 
    	Fgmp::execute("tblastn -db $genome_file.db -query $prot -word_size 5 -max_target_seqs 5 -evalue 0.01 -seg yes -num_threads $cps -outfmt  \"7 sseqid sstart send sframe bitscore qseqid\" > $genome_file.tblastn");
	#Fgmp::execute("grep -v \'#\' $genome_file.tblastnOut | cut -f 1 | sort -u > $genome_file.candidates");
	my (%adjusted) = Fgmp::extractCandidateRegion("$genome_file.tblastn");	

	# with blastall 
	#Fgmp::execute("formatdb -i $genome_file -p F");
	#Fgmp::execute("blastall -p tblastn -d $genome_file -i $prot  -v 5 -b 5 -a $cps -e 0.01 -m 8 -o $genome_file.tblastn");
        #my (%adjusted) = Fgmp::extractCandidateRegion("$genome_file.tblastn"); 
	
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
		"t|tag=s"		=> \$tag,
		"v|verbose"		=> \$verbose_flg,	# verbose
		"q|quiet"		=> \$quiet_flg,	# quiet mode
		"tmp"			=> \$temp_flg,	# keep tempory files
		"h|help|?"		=> \$help_flg,	# print help
		"T|threads=i"		=> \$threads, 	# number of threads 
		"A|augTraingCutoff=i"	=> \$augTraingCutoff,
	) || &show_help(); 
	&show_help if $help_flg;
	
	
	# check if files exists
	&die("FATAL ERROR!!! --genome or --blastdb not specified")
		unless (defined($genome) || defined($blastdb));

	# need to check here that the $genome is in the correct fasta format
	#&die("FATAL ERROR!!! --verbose and --quiet are mutually exclusive")
	#	if ($verb_flg && $quiet_flg);
	

	# default settings
#	$protein	= "$FGMP/data/OneRep.fa" if (!(defined($protein)));
	$protein        = "$FGMP/data/593_cleanMarkers.fa" if (!(defined($protein)));
	$output 	= "output" if (!(defined($output))); 
	$hmm_profiles	= "$FGMP/data/593_cleanMarkers.hmm" if (!(defined($hmm_profiles)));
	$hmm_prefix	= "OMA" if (!(defined($hmm_prefix)));
	$cutoff_file	= "$FGMP/data/profiles_cutOff.tbl" if (!(defined($cutoff_file)));
	$mark_file	= "$FGMP/data/593_cleanMarkers.txt" if (!(defined($mark_file)));
	$fuces_hmm	= "$FGMP/data/172_fUCEs.hmm" if (!(defined($fuces_hmm)));
	$fuces_prefix   = "$FGMP/data/172_fUCEs.txt" if (!(defined($fuces_prefix))); 
	$tag		= "OMA" if (!(defined($tag)));
	$verbose_str	= " -v " if $verbose_flg;
	$threads	= 4 if (!(defined($threads)));
	$augTraingCutoff = 50 if (!(defined($augTraingCutoff)));
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

		xxxxxxxxxxxxxxx

USAGE

	$SOFTWARE [options] -g < genome_fasta_file>

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

	-t, --tag		tag to use OMA for fgmp, FUNY (Funybase) or CEG (cegma)

	-T, --threads		Specify the number of processor threads to use

	-v, --verbose		show progress

	-q, quiet		suppress show log

	-h, --help		show this help

	--tmp			keep temporary files

	-augTraingCutoff	specify the num of genes for augustus training

BUGS:
	Please report bugs to 'ousmane.cisse\@ucr.edu'.

AUTHORS:

	$SOFTWARE has been developped by Ousmane H. Cisse and Jason E. Stajich.



GNU-GPL (C) 		date				 $SOFTWARE
+++EndOfHelp+++
	close (HELP);
	exit(1);
}
