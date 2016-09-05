#!/bin/bash

## Nicolas Servant
## Copyright (c) 2016 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

set -o pipefail
set -o errexit 

read_config()
{
    eval "$(sed -e '/^$/d' -e '/^#/d' -e 's/ =/=/' -e 's/= /=/' $1 | awk -F"=" '{printf("%s=\"%s\"; export %s;\n", $1, $2, $1)}')"
}

## Usage
function usage {
    echo -e "Usage : ./5C_pro.sh"
    echo -e "-i"" <input>"
    echo -e "-c"" <config>"
    echo -e "-o"" <output directory/prefix>"
    echo -e "-h"" <help>"
    exit
}

mode=global
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF=$2; shift;;
	(-i) INPUT=$2; shift;;
	(-o) ODIR=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

## Read configuration files
read_config $CONF

prefix=$(basename $INPUT | sed -e 's/.fastq\(.gz\)*//')
cinput=${ODIR}/trimming/${prefix}_trim.fastq


###################################
## 2- Clean input data
###################################
if [[ $DO_TRIM == 1 ]]; then

    if [[ -z ${T3_ADAPTER} || -z ${T7_ADAPTER} ]]; then
	echo "Please specify T3 and T7 adapter sequence. Exit"
	exit
    fi
    
    echo "Remove T3/T7 adapter ..."
    
    mkdir -p ${ODIR}/trimming
    
    T3_ADAPTER_RC=$(echo ${T3_ADAPTER} | rev | tr ATGC TACG)
    T7_ADAPTER_RC=$(echo ${T7_ADAPTER} | rev | tr ATGC TACG)
    
    cmd="${CUTADAPT_PATH}/cutadapt -g ${T3_ADAPTER} -a ${T7_ADAPTER_RC} -n 2 -m ${MIN_LENGTH} -o ${cinput} ${INPUT} &> ${ODIR}/trimming/cutadapt_${prefix}.log"
    eval $cmd
fi

###################################
## 3- Run bowtie mapping
###################################
if [[ $DO_MAPPING == 1 ]]; then
    if [[ -z ${BOWTIE2_IDX_PATH} ]]; then
	echo "Please specify bowtie2 index or fasta file for indexing. Exit"
	exit
    fi

    ## Check mapping options
    if [[ -z ${BOWTIE2_OPTIONS} ]]; then
	echo "Mapping options not defined. Exit"
	exit 1
    fi
    
    if [[ -z ${BOWTIE2_IDX_PATH} ]]; then
	echo "Bowtie2 Index is not defined. Exit"
	exit 1
    fi
    
    if [[ ! -e ${cinput} ]]; then
	echo "$cinput not find for mapping. Exit"
	exit 1
    fi

    prefix=$(basename ${INPUT} | sed -e 's/.fastq\(.gz\)*//')

    ## Output
    mkdir -p ${ODIR}/mapping
    echo "Align reads on primers indexes ..."

    ## Run bowtie
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX_PATH} -U ${cinput} 2> ${ODIR}/mapping/bowtie_${prefix}_bowtie2.log | ${SAMTOOLS_PATH}/samtools view -bS - > ${ODIR}/mapping/${prefix}_bwt2.bam"
    eval $cmd
fi


########################################################
## Build maps
########################################################

if [[ $DO_BUILD == 1 ]]; then
    
    echo "Build maps from BAM file ..."
    mkdir -p ${ODIR}/maps

    cmd="${PYTHON_PATH}/python ./scripts/bam2maps.py -i ${ODIR}/mapping/${prefix}_bwt2.bam -o ${ODIR}/maps/${prefix}_bwt2_rf.matrix -q ${MIN_MAPQ} -v > ${ODIR}/maps/${prefix}_bwt2_bam2maps.log"
    eval $cmd
fi


########################################################
## First QC
########################################################

if [[ $DO_QC == 1 ]]; then
    
    echo "Run QC ..."
    mkdir -p ${ODIR}/qc
    mkdir -p ${ODIR}/data

    cmd="${R_PATH}/R --no-save --no-restore -e \"input='${ODIR}/maps/${prefix}_bwt2_rf.matrix'; fbed='${FORWARD_BED}'; rbed='${REVERSE_BED}'; chr='chrX'; blist='${PRIMERS_BLACKLIST}'; odir='${ODIR}'; rmarkdown::render('./scripts/qualityControl.Rmd', output_file='${ODIR}/qc/${prefix}_qc.html')\"" 
    echo $cmd
    eval $cmd

fi
