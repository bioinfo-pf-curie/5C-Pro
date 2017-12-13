#!/bin/bash

set -o pipefail                     
set -o errexit

## Generate the annotation file
../bin/prepare_reference.sh -i ./Galupa_et_al/Xic3_MAIN_primer_pool_mm9.tsv -c ./config.txt -o annotations

## Run the pipeline on the test dataset
../bin/5C_pro.sh -i Galupa_et_al/input.fastq.gz -c config.txt -o test-op-res
