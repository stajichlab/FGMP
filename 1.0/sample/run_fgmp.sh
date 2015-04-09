export FGMP=/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0
export PERL5LIB="/opt/perl/5.16.3/packages/lib/perl5:/opt/perl/5.16.3/packages/lib/site_perl:/opt/arch_independent/share/perl/5.8.8:/opt/illumina-GAPipeline-1.3.2/perl/lib:/rhome/ocisse/tcoffee/Version_10.00.r1613/perl:/rhome/ocisse/tcoffee/Version_10.00.r1613/perl:/bigdata/ocisse/Project3_cema/Version1/src/fgmp/1.0/lib:$FGMP/lib"
cd ~/bigdata/Project3_cema/Version1/src/fgmp/1.0/sample
perl ../src/fgmp.pl -g sample.dna -p sample.prot -T 9
