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

## How to use it ?

1. Prepare all annotations files from the primer design

./prepare_reference.sh -i test_data/primer_pool_mm9.tsv -c config.txt -o annotations

2. Run the pipeline to generate the 5C maps

./5C_pro.sh -i test_data/input.fastq.gz -c config.txt -o res


## Output files

The outputs of each step are stored in a dedicated folder, i.e. trimming, mapping, qc, maps.

## Run it in command line for all the samples

for i in `seq 1 28`
do 
SAMPLE=$(awk -F"," -v i=${i} '$1==i{print $2}' /bioinfo/users/nservant/projects_analysis/NGS/5C_galupa/data_ln/samples.txt); 
INPUT=$(awk -F"," -v i=${i} '$1==i{print $3}' /bioinfo/users/nservant/projects_analysis/NGS/5C_galupa/data_ln/samples.txt); 
PROJPATH=/data/kdi_prod/.kdi/project_workspace_0/1066/acl/01.00/
echo $SAMPLE
bash ${PROJPATH}/scripts/5C_pro.sh -i ${PROJPATH}/data_ln/${INPUT} -c ${PROJPATH}/scripts/config_centos.txt -o  ${PROJPATH}/results/${SAMPLE}
done
