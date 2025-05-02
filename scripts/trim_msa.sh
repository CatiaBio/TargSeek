#!/bin/bash

input_folder="results/gram_negative/fasta_msa"
output_folder="results/gram_negative/fasta_msa_trimmed"

mkdir -p "$output_folder"

for file in "$input_folder"/*.fasta
do
    filename=$(basename "$file")
    trimal -in "$file" -out "$output_folder/$filename" -automated1
done
