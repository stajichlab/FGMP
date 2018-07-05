#!/bin/bash

cat $1  | perl -ne '
if (/Target range:/){
        <>;<>;
        while ($_ !~ /^hreads/){
        <>;
        chomp($target_aa .= <>);
        chomp($target_na .= <>);
        <>;
        chomp($_=<>);
        }
        $target_na =~ s/[a-z\s\:\d\-\.]//g;
        $target_na =~ s/{.+?}//g;
        $target_na =~ s/[\s\-\+]//g;
        $target_na =~ s/{.+?}//g;
        print "$_\t$target_aa\t$target_na\n";
        $target_aa = $target_na = "";
}
' > $1.aa

