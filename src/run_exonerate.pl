#!/usr/bin/env perl

use strict; 
use warnings;
use feature 'say'; 
use Carp; 
use FGMP; 

my ($candidateRegion,$protein) = @ARGV;
FGMP::run("exonerate --model protein2genome --percent 25 -q $protein --showtargetgff Y -t $candidateRegion --showvulgar F --showalignment T --ryo \'%qi,%ql,%qab,%qae,%ti,%tl,%tab,%tae,%et,%ei,%es,%em,%r,%pi,%ps,%C\n\' > $candidateRegion.p2g");

