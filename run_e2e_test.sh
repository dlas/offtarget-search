#!/usr/bin/bash

./test.py --fasta_reference_path=../data/tiny_chr21.fa.gz --vcf_augmentation_path=../data/tiny_vcf.vcf.bgz --contig=chr21 --pam=AT --spacer=TTGTGACTGAAGGGC --pdf_output="foo.pdf" > temp_test_file.txt

cmp temp_test_file.txt test_data/expected_test_output.txt
if [ $? == 0 ];  then
	echo TEST PASS
else
	echo TEST FAIL
fi

