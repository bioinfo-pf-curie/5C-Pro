## Nicolas Servant
## Pre-process 5C primers information

## Extract FOR and REV primers sequences

awk 'NR>1{print ">"$1"_mm9_"$5"_"$15"_"$16} $3=="FOR"{print $21$10} $3=="REV"{print $10$21}' archive/Xic3_MAIN_primer_pool_mm9.tsv | grep -A1 "FOR">  Xic3_MAIN_primer_pool_mm9_FOR.fasta 
awk 'NR>1{print ">"$1"_mm9_"$5"_"$15"_"$16} $3=="FOR"{print $21$10} $3=="REV"{print $10$21}' archive/Xic3_MAIN_primer_pool_mm9.tsv | grep -A1 "REV">  Xic3_MAIN_primer_pool_mm9_REV.fasta 

## Make BED of FOR/REV primers
awk '$1~">"{gsub(">","",$0);split($0,sl,"_"); print sl[5]"\t"sl[6]"\t"sl[7]"\t"sl[1]"-"sl[2]"-"sl[3]}' Xic3_MAIN_primer_pool_mm9_FOR.fasta > Xic3_MAIN_primer_pool_mm9_FOR.bed
awk '$1~">"{gsub(">","",$0);split($0,sl,"_"); print sl[5]"\t"sl[6]"\t"sl[7]"\t"sl[1]"-"sl[2]"-"sl[3]}' Xic3_MAIN_primer_pool_mm9_REV.fasta > Xic3_MAIN_primer_pool_mm9_REV.bed


## Create fasta of all pairwise information

function mergeFasta {
    while read  name1; do
	read  seq1
	while read  name2; do
	    read  seq2
	    if [[ $name1 != "" && $name2 != "" ]]; then
		echo $name1"-"$name2 | sed -e 's/>//g'
		echo $seq1$seq2
	    fi
	done < "$2"
    done < "$1"
}

mergeFasta Xic3_MAIN_primer_pool_mm9_FOR.fasta Xic3_MAIN_primer_pool_mm9_REV.fasta > Xic3_MAIN_primer_pool_paired.fasta
perl -pi -e 's/^Xic/>Xic/g' Xic3_MAIN_primer_pool_paired.fasta