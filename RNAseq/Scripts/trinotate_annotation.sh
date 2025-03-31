#!/bin/bash
fasta_file=$1
expression_data=$2
output_dir=$3

mkdir -p $output_dir
Trinotate Trinotate.sqlite init --gene_trans_map $fasta_file --transcript_fasta $fasta_file --transdecoder_pep $fasta_file.transdecoder.pep
Trinotate Trinotate.sqlite report > $output_dir/trinotate_annotation_report.txt
