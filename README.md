# FGMP 
Fungal Genome Mapping Pipeline

## Contents

+ 1. FGMP description ?
+ 2. Installing FGMP
+ 3. File Listing
+ 4. running FGMP
+ 5. Testing FGMP
+ 6. Authors and help
+ 7. Citing FGMP

----------------------------------------
## Introduction

FGMP (Fungal Genome Mapping Project) is a bioinformatic pipeline designed to 
provide in an unbiased manner an estimation of genome completeness of a fungal
genome assembly. The strategy is based on the screening of the genome using a 
set of highly diversified fungal proteins. This approach is likely to capture homologs from any 
fungal genome. FGMP is based on 593 protein markers and 31 highly conserved fungal genomic segments .  

A local version of FGMP can be installed on UNIX platforms. The tool requires the
pre-installation of Perl, NCBI BLAST, HMMER, EXONERATE and AUGUSTUS. 

The pipeline uses information from the selected genes of 25 fungi by first using TBLASTn to 
identify candidate regions in a new genome. Gene structures are delineated using EXONERATE and AUGUSTUS
and validated using HMMER. At the end of the process FGMP produces a set of best predictions and an estimation of the genome 
completeness. 

FGMP uses NHMMER to screen the genome with a set of highly conserved DNA segments. When raw reads are provided in fasta format, 
FGMP searches protein markers directly in the unassembled reads. Reads are randomly sampled using reservoir sampling
algorithm and screened using BLASTx.

FGMP source code and documentation are available under the GNU GENERAL PUBLIC LICENSE.

## Installation
+ System requirements
	- Perl 5 (tested with the version 20)
	- BioPerl-1.6.924 http://bioperl.org
	- HMMER v3.0    http://hmmer.org/
	- NCBI BLASTALL (tested using version 2.2.31+) ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/
	- Exonerate (tested using version 2.2.0) https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
	- Augustus (tested using version 3.0.3) http://augustus.gobics.de/

These software can be installed via bioconda

```shell
conda install hmmer
conda install blast
conda install exonerate
conda install augustus
conda install emboss
conda install perl-ipc-run
conda install -c bioconda perl-bioperl
```

Note: Augustus version 3.3 might cause compiling issues on Mac OS (but not on linux). We recommend using the version 3.2.3 for Mac.

Install
```shell
echo "Installation via github"
git clone https://github.com/stajichlab/FGMP.git
```

## Files listing

+ Description
	The FGMP distribution includes the following files and directories:

	- data/				Proteins, profiles and cutoff
	- lib/				Perl modules
	- sample/			Sample dna
	- sample_output/		Sample output
	- src/				Source code of FGMP.	
	- utils/			FGMP scripts

+ running FGMP

To use FGMP with default settings run:
```shell
	fgmp.pl -g < genomic_fasta_file > fgmp_report.out
```

You can specify the number of cpus to use using the -T option, which will be passed
to all subsequent softwares.

	-p, --protein		fasta file of the protein sequence
				(default: $FGMP/data/593_cleanMarkers.fa)

	--fuces_hmm		Directory that contains hmm files
				(default: $FGMP/data/172_fUCEs.hmm)

	-c, --cutoff_file	Profile cutoffs
				(default: $FGMP/data/profiles_cutOff.tbl)
				
	--hmm_profiles		Directory that contains hmm files
				(default: $FGMP/data/593_cleanMarkers.hmm)

	
+ Testing

launch the following command and compare the output with the sample files in 'sample_output'
```shell
	./fgmp.pl -g sample_test.dna -T 3 --tag OMA
```

+ Output 	
FGMP will create some intermediate files during the annotation.

final files:
	- sample.dna.bestPreds.fas: predicted best predictions (fasta format)

	- sample.dna.unfiltered.renamed.hmmsearch.full_report: detailed analysis of best predictions

	- sample.dna.unfiltered.renamed.hmmsearch.summary_report: summary

intermediate files: 
	- sample.dna.tblastn: 	tblastn output

	- sample.dna.candidates.fa: Genomic regions extracted based on Tblastn matches coordinates (fasta format)

	- sample.dna.candidates.fa.p2g: Alignment of 593 proteins to candidates.fa

	- sample.dna.candidates.fa.p2g.aa: exonerate alignment matches

	- sample.dna.candidates.fa.p2g.aa.proteins: translated CDS (amino acids)

	- sample.dna.trainingSet: augustus training set

	- sample.dna.trainingSet.gb: augustus training set (genbank format)

	- sample.dna.unfiltered: unfiltered predicted peptides

	- sample.dna.unfiltered.renamed : renamed predicted peptides to avoid name conflits

	- sample.dna.unfiltered.renamed.hmmsearch : Hmmsearch output

+ Search in reads (experimental)
Reads can be searched using the following command:
```shell
	perl $FGMP/src/fgmp.pl -g Ncrassa_V2.fasta -T $CPU -r sd_merge.fq.fasta
```

FGMP main repository page is at https://github.com/stajichlab/FGMP

Development:  https://github.com/ocisse/FGMP
(original code and changes to be merged to this repository where needed)
 
## Citation
Cisse, O. H. and Stajich J.E. FGMP: assessing fungal genome completeness and gene content.
bioRxiv 049619; doi: https://doi.org/10.1101/049619 (2016).
