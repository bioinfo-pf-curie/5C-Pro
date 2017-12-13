Author: Nicolas Servant  
Date: 11-07-2016  
Version: 0.0.1  

# 5C Pro Quickstart Guide

<!-- This page is a quick start guide, please read the full `online manual <link>`_ for more information. -->

See NEWS for information about changes in this and previous versions

See LICENSE for license information


## What is it ?

The pipeline was designed to process 5C data. It is currently designed for single ends sequencing reads.
The following workflow is applied ;
* Adapter trimming of input fastq files
* Alignment on the pseudo-reference
* Count of 3C product based on the alignment results
* Build raw contact maps at the restriction fragment resolution
* Quality controls


## Input Files

1. Configuration file

All the parameters and paths are defined in the configuration file.
Edit the 'config.txt' file and define the different parameters.

## Build annotations

Prepare all annotations files from the primer design.

Primers input files must contain the following information :
PRIMER_NAME CHROMOSOME P_STARTPOS P_ENDPOS TYPE P_SPECIFICPRIMER_SEQUENCE

./bin/prepare_reference.sh -i test_op/Galupa_et_al/Xic3_MAIN_primer_pool_mm9.tsv -c config.txt -o test-op/annotations

The prepare_reference.sh script extract the primer sequences, and build the bowtie2 indexes for the mapping

## How to use it ?

2. Run the pipeline to generate the 5C maps

./bin/5C_pro.sh -i ./test-op/Galupa_et_al/input.fastq.gz -c ./test-op/config.txt -o ./test-op/res


## Output files

The outputs of each step are stored in a dedicated folder, i.e. trimming, mapping, qc, maps.


