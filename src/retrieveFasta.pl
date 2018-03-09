#!/usr/bin/env perl

use strict;
use warnings;
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
<<<<<<< HEAD:src/retrieveFasta.pl
} 
=======

}
>>>>>>> stajichlab/master:1.0/src/retrieveFasta.pl

$new=Bio::SeqIO->new(
			-file=> $proteomdb,
			-format=>"fasta");

while ($seq=$new->next_seq){
	if (defined $id_hash{$seq->id}){
		print ">", $seq->id,"\n",$seq->seq,"\n";
	}
}
