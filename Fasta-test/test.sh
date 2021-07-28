#!/bin/sh -e

if [ $0 != ./test.sh ]; then
    printf "Must be run as ./test.sh.\n"
    exit 1
fi

cd ..
./cave-man-install.sh
cd Fasta-test

printf "FASTA test:\n\n"
cc -o fasta-test fasta-test.c -I../../local/include \
    -L../../local/lib -Wl,-rpath,../../local/lib -lbiolibc -lxtend
./fasta-test < test.fasta > out.fasta
if diff test.fasta out.fasta; then
    printf "No differences found, test passed.\n"
else
    printf "Differences found, test failed.\n"
fi
rm -f fasta-test out.fasta
