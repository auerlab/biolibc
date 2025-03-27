#!/bin/sh -e

if [ $0 != ./test.sh ]; then
    printf "Must be run as ./test.sh.\n"
    exit 1
fi

cd ..
./cave-man-install.sh
cd Fastq-test

printf "FASTQ test:\n\n"
cc -o fastq-test fastq-test.c -I../../local/include \
    -L../../local/lib -Wl,-rpath,../../local/lib -lbiolibc -lxtend

printf "Min qual 0...\n"
./fastq-test 0 < correct.fastq > out.fastq
if diff correct.fastq out.fastq; then
    printf "No differences found, test passed.\n"
else
    printf "Differences found, test failed.\n"
fi

min_qual=20
printf "\n===\nMin qual $min_qual...\n"
./fastq-test $min_qual < low-qual.fastq > out.fastq
cat out.fastq
if diff low-qual-correct.fastq out.fastq; then
    printf "No differences found, test passed.\n"
else
    printf "Differences found, test failed.\n"
fi

rm -f fastq-test out.fastq
