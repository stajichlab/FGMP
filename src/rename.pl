#!/opt/linux/centos/7.x/x86_64/pkgs/perl/5.20.2/bin/perl -w

use Carp;
use strict;
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

