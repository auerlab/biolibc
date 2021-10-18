#!/bin/sh -e

if [ $0 != ./test.sh ]; then
    printf "Must be run as ./test.sh.\n"
    exit 1
fi

cd ..
./cave-man-install.sh
cd ORF-test

printf "ORF test:\n\n"
cc -o orf-test orf-test.c -I../../local/include \
    -L../../local/lib -Wl,-rpath,../../local/lib -lbiolibc -lxtend
./orf-test < orf.txt > out.txt
if diff correct.txt out.txt; then
    printf "No differences found, test passed.\n"
else
    printf "Differences found, test failed.\n"
fi
rm -f orf-test out.txt
