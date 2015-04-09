# FGMP v 1.0
Fungal Genomics Mapping Pipeline

Released: 4/8/25

# *** Contents *** #

+ 1. FGMP description ?
+ 2. Installing FGMP
+ 3. File Listing
+ 4. running FGMP
+ 5. Testing FGMP
+ 6. Authors and help
+ 7. Citing FGMP

# ---------------------------------------- #
# 1. FGMP description

FGMP (Fungal Gene Mapping Project) is a bioinformatic pipeline designed to 
provide in unbiased manner a set of high quality gene models from any fungal
genome assembly. The strategy is based on the screening of the genome using a 
set of highly diversified fungal proteins, that is espected to represent a realistic 
snapshot of fungal diversity). This approach is likely to capture homolog from any 
fungal genomes. FMP combines ab initio predictions and translated matches to produce 
its final predictions.
Using an empirical approach, we identified 245 protein markers that arose early in eukaryotes, 
are conserved in all fungi. This set of genes display  50% overlap
with CEGMA 248 eukaryotic markers.  

A local version of FGMP can be installed on UNIX platforms and it requires
pre-installation of Perl, NCBI BLAST, HMMER, exonerate and augustus. 

The pipeline uses information from the selected genes of 32 fungi by first using TBLASTN to 
identify candidate regions in a new genome. Gene structures are first delinated using exonerate and augustus,
and validated using HMMER. At the end FGMP produces a set of best predictions and an estimation of the genome 
completeness of the genome analyzed. 

FGMP source code and documentation are available under the GNU GENERAL PUBLIC LICENSE.

# ---------------------------------------- #
# 2. running FGMP 

The FGMP distribution includes several directories and files. Source code 
as well as documented are provided in the distribution. To decompress 
simply enter the following command in the shell: 

CMD:	tar -xvzf FGMP.tar.gz

The directory 'fgmp_v1.0' will be created in your current directory. 

FGMP requieres the pre-installation of the following software: 

- hmmer (HMMER v3.0)	http://hmmer.janelia.org/
- NCBI BLASTALL (tested using version 2.2.21)
- exonerate (tested using version 2.2.0)
- augustus (tested using version 3.0.3)


# ---------------------------------------- #
# 3. Files listing

The FGMP distribution includes the following files and directories:

- data/	Proteins, profiles and cutoff.
- lib/	Perl modules
- sample/	DNA and protein fasta file to test FGMP.
- sample_output/	Results provided by FGMP for the test files.
- src/	Source code of FGMP.
- GNULicence	This software is registed under the GNU licence.	
- READE.md	This file.
- fgmp.config	a file that needs to be set up with path directories
- utils/	contains augustus and some useful scripts (see below)

- countDomains.pl - a parser that count the number of domains in the output of HMMER
- simplifyNcbiHeader.pl - a parser to remove unecessary information from NCBI contigs header. 

# ---------------------------------------- #
# 4. running FGMP

**** IMPORTANT ****

IMPORTANT - the 'fgmp.config' file needs to be placed where are the fasta files

In the file 'fgmp.config' there several environmental variables
that need to be set up by the users.

It contains the following variables:

FGMP	- the path to the FGMP folder
WRKDIR	- working directory	(where are the fasta files)

For example: 
FGMP=/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0
WRKDIR=/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0/sample

FGMP uses a custom library Fgmp.pm. You need set the PERL5LIB environment variable 
to use or you can simply copy the modules to the Perl module directory that is 
available to your Perl installation.

Before running FGMP, type the following commands:
	
	export FGMP=/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0
	export PERL5LIB="$PERL5LIB:$FGMP/lib"

To run FGMP using the 7773 seed proteins and genome completeness evaluation using the
245 fungal markers run:
	fgmp.pl -g < genomic_fasta_file > 

You can specify the number of cpus to use using the -T option, which will be passes
to all subsequent softwares.

	-p, --protein		fasta file of the protein sequence
				(default: \$FGMP/data/OneRep.fa)
	-c, --cutoff_file	Profile cutoffs
				(default: \$FGMP/data/profiles_cutOff.tbl)
	--hmm_profiles		Directory that contains hmm files
				(default: \$FGMP/data/hmm_profiles)
	

Example:	fgmp.pl -g sample.dna -p sample.prot --hmm_profiles hmm_profiles \
		-c profiles_cutOff.tbl

If you have another set of proteins you want use instead, simply provide them. 


# ---------------------------------------- #
# 5. TESTING FGMP

launch the following command and compare the output with the sample files in 'sample_output'

	 perl ../src/fgmp.pl -g sample.dna -p sample.prot 2> log

# running on cluster (example - must be a better way)
	
	echo "export FGMP=/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0" > run_fgmp.sh	
	echo "export PERL5LIB="$PERL5LIB:$FGMP/lib"" >> run_fgmp.sh
	echo "cd ~/bigdata/Project3_cema/Version1/src/fgmp/1.0/sample" >> run_fgmp.sh
	echo "perl ../src/fgmp.pl -g sample.dna -p sample.prot -T9" >> run_fgmp.sh		
	qsub -l nodes=1:ppn=9 run_fgmp.sh

FGMP creates some intermediate files during the annotation.

# - final files ( to be kept)
- sample.dna.bestPreds.fas	: predicted best predictions (fasta format)
- sample.dna.unfiltered.renamed.hmmsearch.full_report	: detailed analysis of best predictions
- sample.dna.unfiltered.renamed.hmmsearch.summary_report : summary

# - intermediate files : (likely to be removed)
- sample.dna.tblastn : 	tblastn output
- sample.dna.candidates.fa	: Genomic regions extracted based on Tblastn matches coordinates (fasta format)
- sample.dna.candidates.fa.p2g	: Alignment of 7773 proteins to sample.dna.candidates.fa
- sample.dna.candidates.fa.p2g.aa : intermediary file contains exonerate alignment matches
- sample.dna.candidates.fa.p2g.aa.proteins	: translated CDS (amino acids)
- sample.dna.trainingSet	: augustus training set
- sample.dna.trainingSet.gb	: augustus training set ( enebank format)
- sample.dna.unfiltered : merged fasta files containing both exonerate CDS and augustus predictions
- sample.dna.unfiltered.renamed : renamed 'sample.dna.unfiltered' to avoid name clashes
- sample.dna.unfiltered.renamed.hmmsearch : Hmmsearch output

If you want to keep intermediary files use the --keep option ( working on that ...)

# ---------------------------------------- #
# 6. AUTHORS AND HELP

FGMP has been developped by Ousmane H. Cisse (ousmanecis@gmail.com) and Jason E. Stajich 
(jason.stajich@ucr.edu).
 
FGMP home page is at https://github.com/stajichlab/FGMP

# ---------------------------------------- #
# 7. Citing FGMP

in preparation.
