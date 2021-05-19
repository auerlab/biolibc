#!/bin/sh -e

for file in good; do
    ./bed-test < $file.bed > out.bed
    if diff $file.bed out.bed; then
	printf "BED test: Processed $file.bed OK.\n"
    else
	printf "BED test: Failure on $file.bed.\n"
    fi
done
