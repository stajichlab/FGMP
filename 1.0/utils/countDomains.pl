#!/opt/perl/5.16.3/bin/perl -w

use strict; 
use Carp; 
use feature 'say'; 
use Data::Dumper;
use IO::All;


# record all Pfam domains

my $file = shift; 

my @domains = (); 

open F, $file || croak "cannot open $file:$!\n";
while(<F>){
chomp;
	push(@domains, $_);	
}

close F;

my %count = (); 
my $i = ''; 

foreach $i (@domains){
	$count{$i}++; 
}

# 
my $x = ''; 
foreach $x (keys %count){
	say "$x\t$count{$x}";
}

#say Dumper \%count; 
