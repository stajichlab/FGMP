#!/opt/linux/centos/7.x/x86_64/pkgs/perl/5.20.2/bin/perl -w

use Data::Dumper; 
use Carp; 
use feature 'say';
use Bio::SeqIO;

my ($list,$proteomdb) = @ARGV;
my %id_hash = ();

open (FILE,$list) || croak "cannot open $list:$!\n";
while(<FILE>){
chomp; 
	$id_hash{$_} = 1;
} 

$new=Bio::SeqIO->new(
			-file=> $proteomdb,
			-format=>"fasta");

while ($seq=$new->next_seq){
	if (defined $id_hash{$seq->id}){
		print ">", $seq->id,"\n",$seq->seq,"\n";
	}
}
