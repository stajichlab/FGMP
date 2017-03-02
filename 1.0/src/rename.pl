#!/usr/bin/perl -w

use Carp;
use strict;
use warnings;
use feature 'say'; 
use Data::Dumper; 


my $fasta = shift; 
my $i = 0;

open F,$fasta || croak "cannot open $fasta:$!\n";
while(<F>){
chomp; 
	if (/^>/){
		$_ =~s/>//;
		say">Seq_$i\t$_";
		$i++
	} else {
		say;
	}
}
close F;

