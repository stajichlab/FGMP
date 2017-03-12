#!/opt/linux/centos/7.x/x86_64/pkgs/perl/5.20.2/bin/perl -w

use strict; 
use feature 'say'; 
use Carp; 
use Fgmp; 

my ($candidateRegion,$protein) = @ARGV;
Fgmp::run("exonerate --model protein2genome --percent 25 -q $protein --showtargetgff Y -t $candidateRegion --showvulgar F --showalignment T --ryo \'%qi,%ql,%qab,%qae,%ti,%tl,%tab,%tae,%et,%ei,%es,%em,%r,%pi,%ps,%C\n\' > $candidateRegion.p2g");

