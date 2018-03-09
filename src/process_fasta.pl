#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $in_file = $ARGV[0];
my $start_pos = $ARGV[1];
my $end_pos = $ARGV[2];

my $in = Bio::SeqIO->new ( -file => $in_file, -format => 'fasta');
my $out = Bio::SeqIO->new( -file => ">$in_file.out", -format => 'fasta');
while (my $seq = $in->next_seq() ) {
    $seq->display_id( $seq->display_id() . "_$start_pos-$end_pos" );
	if (($seq->length) > $end_pos){
    		$out->write_seq( $seq->trunc($start_pos, $end_pos) );
	} else {
		$out->write_seq( $seq->trunc($start_pos,$seq->length));
	}
}
