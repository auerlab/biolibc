#!/bin/sh -e

cd ..
./cave-man-install.sh
cd GFF-test
cc -I.. -o gff-test gff-test.c -L.. -lbiolibc -L../../local/lib -lxtend
cat << EOM

Terminal output should be:

gene:ENSDARG00000070713 prss60.1
gene:ENSDARG00000055644 prss60.2
gene:ENSDARG00000070710 prss60.3
gene:ENSDARG00000077308 gpr84
gene:ENSDARG00000031647 stat2
gene:ENSDARG00000090980 apof
gene:ENSDARG00000097596 tac3b
gene:ENSDARG00000055561 c1galt1b
gene:ENSDARG00000020301 os9

EOM
for file in good; do
    ./gff-test < $file.gff3 > out.gff3
    echo ''
    if diff $file.gff3 out.gff3; then
	printf "GFF3 test: Processed $file.gff3 OK.\n"
    else
	printf "GFF3 test: Failure on $file.gff3.\n"
    fi
done
