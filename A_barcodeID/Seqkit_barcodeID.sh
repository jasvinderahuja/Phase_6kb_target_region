#!/bin/bash
#SBATCH --job-name seqkit                                                                             
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --gres=lscratch:100
#SBATCH -t 1:00:00
#SBATCH -o R_%j.out                                                                              
#SBATCH -e R_%j.err


echo "command: " $0 $1 $2 $3 $4

in_fastq=$1
fwd_barcodes=$2
rev_barcodes=$3

# module load seqkit
## Flag --degenerate is on because patterns contain degenerate bases.
## seqkit locate --degenerate --ignore-case --pattern-file enzymes.fa viral.1.1.genomic.fna.gz

cat $in_fastq \
	| seqkit locate --ignore-case --pattern-file $fwd_barcodes > Fwd_barcode_scan.tsv
cat $in_fastq \
	| seqkit locate --ignore-case --pattern-file $rev_barcodes > Rev_barcode_scan.tsv

