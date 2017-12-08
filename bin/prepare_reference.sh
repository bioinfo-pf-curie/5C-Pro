## Nicolas Servant
## Copyright (c) 2016 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.  
## See the LICENCE file for details

## The iput file is a tab delemited file (.tsv) with the following information ;
## PRIMER_NAME REGION TYPE ASSEMBLY CHROMOSOME FRAGMENT_ID PRIMER_ID P_STARTPOS P_ENDPOS P_SPECIFIC P_SPECIFIC_SIZE P_FILLER P_FILLER_SIZE P_TMP_GC F_STARTPOS F_ENDPOS F_SIZE P_MER P_BLAST BARCODE_NUM BARCODE_SEQ PRIMER_SEQUENCE
## 5C_938_XIC-3_FOR_3 Xic FOR mm9_dna chrX 3 1 98837477 98837506 CGTGTTTTATATTATTGCAGCATTTAAAAG 30 0 48.33624 26.67 98834145 98837506 3361 487 1 966 CGGCTATAATACGACTCACTATAGCCCGGCTACGTGTTTTATATTATTGCAGCATTTAAAAG


set -o pipefail
set -o errexit

read_config()
{
    eval "$(sed -e '/^$/d' -e '/^#/d' -e 's/ =/=/' -e 's/= /=/' config.txt | awk -F"=" '{printf("%s=\"%s\"; export %s;\n", $1, $2, $1)}')"
}

## Usage           
function usage {
    echo -e "Usage : ./prepare_reference.sh"
    echo -e "-i"" <input design - .tsv"
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
        (-i) INPUT_DESIGN=$2; shift;;
        (-o) ODIR=$2; shift;;
        (-h) usage;;
        (--) shift; break;;
        (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
        (*)  break;;
    esac
    shift
done

if [ $# -lt 1 ]
then
    usage
    exit
fi

## Read configuration files
read_config $CONF
prefix=$(basename $INPUT_DESIGN | sed -e 's/.tsv$//')

if [[ ! -e $INPUT_DESIGN ]]; then
    echo "Input file not found"
    exit
fi

mkdir -p ${ODIR}
echo "extract FOR/REV ..."

## Extract FOR and REV primers sequences into fasta file
awk -v org=${ORG} 'NR>1{print ">"$1"_"org"_"$5"_"$15"_"$16} $3=="FOR"{print $21$10} $3=="REV"{print $10$21}' ${INPUT_DESIGN} | grep -A1 "FOR" | grep -v "\-\-" >  ${ODIR}/${prefix}_FOR.fasta 
awk -v org=${ORG} 'NR>1{print ">"$1"_"org"_"$5"_"$15"_"$16} $3=="FOR"{print $21$10} $3=="REV"{print $10$21}' ${INPUT_DESIGN} | grep -A1 "REV" | grep -v "\-\-" >  ${ODIR}/${prefix}_REV.fasta 

echo "make BED ..."

## Make BED of FOR/REV fragments
awk -v org=${ORG} '$1~">"{gsub(">","",$0);split($0,sl,"_"org"_"); split(sl[2],coord,"_"); print coord[1]"\t"coord[2]"\t"coord[3]"\t"sl[1]}' ${ODIR}/${prefix}_FOR.fasta > ${ODIR}/${prefix}_FOR.bed
awk -v org=${ORG} '$1~">"{gsub(">","",$0);split($0,sl,"_"org"_"); split(sl[2],coord,"_"); print coord[1]"\t"coord[2]"\t"coord[3]"\t"sl[1]}' ${ODIR}/${prefix}_REV.fasta > ${ODIR}/${prefix}_REV.bed

## Create fasta of all pairwise information
function mergeFasta {
    while read  name1; do
	read  seq1
	while read  name2; do
	    read  seq2
	    if [[ $name1 =~ ">" && $name2 =~ ">" ]]; then
		header=$(echo $name1"|"$name2 | sed -e 's/>//g')
		echo ">"$header
		echo $seq1$seq2
	    fi
	done < "$2"
    done < "$1"
}

echo "make paired fasta reference ..."
mergeFasta ${ODIR}/${prefix}_FOR.fasta ${ODIR}/${prefix}_REV.fasta > ${ODIR}/${prefix}_paired.fasta
    
## Build Index
echo "Build Bowtie index ..."
mkdir -p ${ODIR}/bwt2
ppool=$(basename ${prefix}_paired.fasta | sed -e 's/.fa\(sta\)*//') 
cmd="${BOWTIE2_PATH}/bowtie2-build ${ODIR}/${prefix}_paired.fasta ${ODIR}/bwt2/$ppool 2> bwt2_build_${ppool}.log "
eval $cmd
    
