#!/bin/sh -e

cd ..
#./cave-man-install.sh
cd Bed-test
cc -I.. -I../../local/include -o bed-test bed-test.c \
    -L.. -lbiolibc -L../../local/lib -lxtend
for file in good; do
    ./bed-test < $file.bed > out.bed
    if diff $file.bed out.bed; then
	printf "BED test: Processed $file.bed OK.\n"
    else
	printf "BED test: Failure on $file.bed.\n"
    fi
done
