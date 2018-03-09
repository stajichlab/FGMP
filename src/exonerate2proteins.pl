#!/usr/bin/env perl -w

use strict; 
use Carp; 
use feature 'say'; 
use IO::All; 

my %IUPCA = ("Ala"=>"A","Asx"=>"B","Cys"=>"C","Asp"=>"D","Glu"=>"E","Phe"=>"F","Gly"=>"G","His"=>"H","Ile"=>"I","Lys"=>"K","Leu"=>"L","Met"=>"M","Asn"=>"N","Pro"=>"P","Gln"=>"Q","Arg"=>"R","Ser"=>"S","Thr"=>"T","Val"=>"V", "Trp" => "W", "Unk" => "X", "Tyr"=>"Y","Glx" => "Z");
my %seqs = ();

my $file = io("$ARGV[0]");
   $file->autoclose(0);
  while (my $line = $file->getline || $file->getline){
	chomp $line; 
	my @data = split /\t/,$line; 
	my $g = shift @data;

	my $seq = "";
	my $el = ""; 
	foreach $el (@data){
		$el =~ s/\s+//g;
		$el =~ s/\+.//g;
		$el =~ s/\-//g;
		$el =~ s/\{+//g;
		$el =~ s/\}+//g; 
		$el =~ s/\d+//g; 
		$el =~ s/\*//g;
		$seq .= $el;
		}
		$seqs{$g} = $seq;
	$seq = "";
	}

my $s = ""; 
foreach $s (keys %seqs){
	my @splittedbyThree = splitByThree($seqs{$s});
	my $translated = translated(@splittedbyThree);
	say ">$s";
	say $translated;
}
	
## sub
sub translated {
	my (@groupedBythree) = @_; 
	
	my @array1 = (); 
	my $sequence = ""; 

	my $i = ""; 
	foreach $i (@groupedBythree){
		if ($IUPCA{$i}){
			my $code =  $IUPCA{$i};
			$sequence .= $code;	
		} else {
		#	warn"####\t$i\n";
		}	

	}
	return($sequence);
}
sub splitByThree {
	my ($in) = @_; 
	
	my @bythree = (); 
	
	my @tmp = split //, $in; 
	my $d = 0; 
	
	while ($d <= @tmp){
		my $n = substr $in,$d,3;
		push(@bythree, $n);	
		$d += 3;
	}
	return @bythree;
}
