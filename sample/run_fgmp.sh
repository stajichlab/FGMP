#PBS -l nodes=1:ppn=4,mem=8gb -N FGMP -j oe -l walltime=200:00:00

module load ncbi-blast/2.2.31+
module load exonerate/2.2.0
module load hmmer/3.1b2
module load emboss/6.6.0
module load perl/5.20.2

CPU=6

if [ $PBS_NUM_PPN ]; then
 CPU=$PBS_NUM_PPN
fi

export FGMP= <FGMP DIRECTORY>
export PERL5LIB="$FGMP/lib"

# run
perl $FGMP/src/fgmp.pl -g sample.fasta --fuces_hmm all.hmm --fuces_prefix all.txt -T 6

# search in reads
perl $FGMP/src/fgmp.pl --nsamples 10 --nsampleSize 100 --reads SRR3110858.fastq
