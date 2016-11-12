use 5.006_001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Fgmp',
    VERSION_FROM      => 'lib/Fgmp.pm’, # finds $VERSION
    # These versions are the lowest ones I've test the module with.
    PREREQ_PM         => {'Test::More'          => '0.45',
                          'Filter::Simple'      => '0.78',
                          'Text::Tabs'          => '98.112801',
                          #'Bio::SeqIO'          =>  '1.007001',
                          #'IO::All'             => '0.86',
                          'IPC::Open3'           =>  '1.13'
                            },
);