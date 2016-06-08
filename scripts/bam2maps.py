#!/usr/bin/env python

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details


"""
Script to pair 2 SAM/BAM files into one PE BAM
- On 03/05/16 Ferhat made changes starting from ~/bin/HiC-Pro_2.7.2b/scripts/mergeSAM.py 
to make singletons possible to be reported
"""

import getopt
import sys
import os
import re
import pysam
from itertools import izip

def usage():
    """Usage function"""
    print "Usage : python mergeSAM.py"
    print "-i/--input <read mapped file>"
    print "-o/--output <Output file. Default is stdin>"
    print "[-q/--qual] <minimum reads mapping quality>"
    print "[-v/--verbose] <Verbose>"
    print "[-h/--help] <Help>"
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:o:q:vh",
            ["input=",
             "output=", "qual=", 
             "verbose", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts


def is_unique_bowtie2(read):
    ret = False
    if not read.is_unmapped and read.has_tag('AS'):
        if read.has_tag('XS'):
            primary =  read.get_tag('AS')
            secondary = read.get_tag('XS')
            if (primary > secondary):
                ret = True
        else:
            ret = True
    
    return ret

## Remove everything after "/" in read's name
def get_read_name(read):
    name = read.qname
    return name.split("/",1)[0]
    

if __name__ == "__main__":
    ## Read command line arguments
    opts = get_args()
    inputFile = None
    outputFile = None
    mapq = None
    verbose = False


    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            inputFile = arg
        elif opt in ("-o", "--output"):
            outputFile = arg
        elif opt in ("-q", "--qual"):
            mapq = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    ## Verbose mode
    if verbose:
        print "## filterBAM.py"
        print "## input=", inputFile
        print "## output=", outputFile
        print "## min mapq=", mapq
        print "## verbose=", verbose

    ## Initialize variables
    if outputFile is None: 
        outputFile = re.sub('\.bam$', '.validPairs', os.path.basename(inputFile))

    reads_counter = 0
    multi_reads_counter = 0
    uniq_reads_counter = 0
    unmapped_reads_counter = 0 
    lowq_reads_counter = 0
    rd = None
    count = {}
    
    ## Reads are 0-based too (for both SAM and BAM format)
    ## Loop on all reads
    if verbose:
        print "## Reading input file ..."
  
    with  pysam.Samfile(inputFile, "rb") as hr : 
        handle_out = open(outputFile, 'w')

        for rd in hr.fetch(until_eof=True):
            reads_counter +=1

            if (reads_counter % 1000000 == 0 and verbose):
                print "##", reads_counter
                
            ## both unmapped
            if rd.is_unmapped == True:
                unmapped_reads_counter += 1
                continue

            ## both mapped
            else:
                ## quality
                if mapq != None and (rd.mapping_quality < int(mapq)):
                    lowq_reads_counter += 1
                    continue
                 
                ## Unique mapping
                if is_unique_bowtie2(rd) == True :
                    uniq_reads_counter += 1
                else:
                    multi_reads_counter += 1
                    continue
            
            ##readname/chr/start/strand/chr/start/strand/length/fname/fname/mapq/mapq
            chrom = hr.getrname(rd.tid)
            schrom=chrom.split("-")
            info1=schrom[0].split("_")
            info2=schrom[1].split("_")

            pname1=info1[0]+"-"+info1[1]+"-"+info1[2]
            pname2=info2[0]+"-"+info2[1]+"-"+info2[2]
            pname=pname1+"|"+pname2

            ## Build dict
            if  count.has_key(pname):
                count[pname]=count[pname] + 1
            else:
                count[pname] = 1

            #handle_out.write(rd.qname+"\t"+info1[3]+"\t"+info1[4]+"\t+\t"+info2[3]+"\t"+info2[4]+"\t+\t"+str(0)+"\t"+pname1+"\t"+pname2+"\t"+str(rd.mapping_quality)+"\t"+str(rd.mapping_quality)+"\n")

    hr.close()
    

    ## Write maps file
    for key in count:
        ks=key.split("|")
        handle_out.write(ks[0] + "\t" + ks[1] + "\t" + str(count[key])+"\n")
    handle_out.close()


    ## stat
    statfile = re.sub('\.bam$', '.pairstat', os.path.basename(inputFile))
    handle_stat = open(statfile, 'w')
            
    handle_stat.write("Total_reads_processed\t" + str(reads_counter) + "\t" + str(round(float(reads_counter)/float(reads_counter)*100,3)) + "\n")
    handle_stat.write("Unmapped_reads\t" + str(unmapped_reads_counter) + "\t" + str(round(float(unmapped_reads_counter)/float(reads_counter)*100,3)) + "\n")
    handle_stat.write("Low_qual_reads\t" + str(lowq_reads_counter) + "\t" + str(round(float(lowq_reads_counter)/float(reads_counter)*100,3)) + "\n")
    handle_stat.write("Unique_reads_alignments\t" + str(uniq_reads_counter) + "\t" + str(round(float(uniq_reads_counter)/float(reads_counter)*100,3)) + "\n")
    handle_stat.write("Multiple_reads_alignments\t" + str(multi_reads_counter) + "\t" + str(round(float(multi_reads_counter)/float(reads_counter)*100,3)) + "\n")
    handle_stat.close()


