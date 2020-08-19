#!/usr/bin/env python
# coding:utf-8
# Author:  Richard Corbett
# Purpose: Parse Strelka output to determine the fraction of the bases that are non reference for each call 
# Created: 11/04/2013

from optparse import OptionParser
import sys
import os
import re
import subprocess
import fnmatch
import numpy

usage = " %prog -v passed.somatic.vcf "

def main():
    
    parser = OptionParser(usage)
    parser.add_option("-v", dest="vcfFile",
                      action="store", type="string",
                      help="path to strelka vcf file.")
    
    #get at the arguments
    (options, args) = parser.parse_args()
    vcfFile = options.vcfFile
    if (vcfFile == None):
        parser.print_help()
        sys.exit()
    
    #Check that our vcf file exists
    if (not os.path.isfile(vcfFile)):
        sys.stderr.write(vcfFile + " is not a valid file\n")
        sys.exit()
    
    #Loops over the lines and pull out the called allele counts
    #Uses the following fields:
    # DP (sample based) "Read depth for tier1 (used+filtered)"
    # FDP "Number of basecalls filtered from original read depth for tier1"
    # AU,CU,TU,GU "Description="Number of 'A|C|T|G' alleles used in tiers 1,2""
    # So using the tier1 ratios: AF=A[0](DP-FDP)
    AFs = []
    base_pos = {'A':0, 'C':1, 'G':2, 'T':3 }
    for line in open(vcfFile):
        #Skip any header lines
        if(line[0] == "#"):
            continue
        
        tokens = line.split("\t")
        ref = tokens[3]
        alt = tokens[4]
        
        #skip any lines with multiple non-ref alleles
        if("," in alt):
            continue
        
        tumour_counts = tokens[10].rstrip()
        DP = float(tumour_counts.split(":")[0])
        FDP = float(tumour_counts.split(":")[1])
        base_counts = tumour_counts.split(":")[4:]
        called_base_t1_count = float(base_counts[base_pos[alt]].split(",")[0]) 

        #print " ".join((ref, alt, str(DP), str(FDP), str(called_base_t1_count), str(called_base_t1_count/(DP-FDP))))
        AFs.append(called_base_t1_count/(DP-FDP))

    #create a histogram of the AFs and print it out
    hist, edges = numpy.histogram(
        AFs,
        bins=20,
        range=(0, 1),
        density=False)
    print "Count\tRange"
    for i, count in enumerate(hist):
        if(edges[i+1]<1):
            print "%s\t%1.2f<=AF<%1.2f" % (count, edges[i], edges[i+1])
        else:
            print "%s\t%1.2f<=AF<=%1.2f" % (count, edges[i], edges[i+1])

   
        
if __name__ == "__main__":
    main()
