#!/usr/bin/env perl

<<<<<<< HEAD
=======

>>>>>>> stajichlab/v2.0-hyphaltip_audit
use strict;
use warnings;
use Data::Dumper;
use Carp;
use feature 'say';
use Bio::SeqIO;

# need better later
my ($list,$proteomdb) = @ARGV;
my %id_hash = ();

open (FILE,$list) || croak "cannot open $list:$!\n";
while(<FILE>){
chomp;
	$id_hash{$_} = 1;
<<<<<<< HEAD
<<<<<<< HEAD:src/retrieveFasta.pl
} 
=======

}
>>>>>>> stajichlab/master:1.0/src/retrieveFasta.pl
=======

}
>>>>>>> stajichlab/v2.0-hyphaltip_audit

$new=Bio::SeqIO->new(
			-file=> $proteomdb,
			-format=>"fasta");

while ($seq=$new->next_seq){
	if (defined $id_hash{$seq->id}){
		print ">", $seq->id,"\n",$seq->seq,"\n";
	}
}
