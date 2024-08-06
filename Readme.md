
Implementation of off-target CRISPR search:

Example:
./main.py `--fasta_reference_path=../data/chr21.fa.gz --vcf_augmentation_path=../data/gnomad_af01_snps_chr21.vcf.bgz --contig=chr21 --pam=CTAN --spacer=ACACGTGTA --pdf_output=op1.pdf`


Limitations:
(1) This program only operates on one chromosome at a time.
(2) This program uses an inefficient in-memory DNA representation (python strings) causing high memory usage and slow performance. 
(3) This program does not explicitly normalize lower-case base names in a fasta or VCF file. 
(4) This program using a simple base-by-base search with no clever alignment or index structure, also causing poor performance. 
(5) No input validation is performed.
