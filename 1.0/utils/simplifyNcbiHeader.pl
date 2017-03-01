#!/usr/bin/env perl

use feature 'say'; 
use Data::Dumper; 
use Carp; 
use strict; 


my $in = shift; 

open I,$in || croak"$!\n";
while(<I>){
chomp; 
	if(/^>/){
		my @line = split /\|/;
		my $id = $line[3];
		say ">$id";

	}else {
		say;
	}
}	
close I;
