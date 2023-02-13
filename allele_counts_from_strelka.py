#!/usr/bin/env python
# coding:utf-8
# Author:  Richard Corbett
# Purpose: Parse Strelka output to determine the fraction of the bases that are non reference for each call
# Created: 11/04/2013
import os
import re
import sys
from optparse import OptionParser


def main():
    parser = OptionParser(" %prog -v passed.somatic.vcf ")
    parser.add_option(
        "-v", dest="vcfFile", action="store", type="string", help="path to strelka vcf file."
    )

    # Get at the arguments
    options, _ = parser.parse_args()
    vcfFile = options.vcfFile
    if not vcfFile:
        parser.print_help()
        sys.exit()

    # Check that our bam file exists
    if not os.path.isfile(vcfFile):
        sys.stderr.write(vcfFile + " is not a valid file\n")
        sys.exit()

    # Loops over the lines and pull out the called allele counts
    basePos = {"A": 0, "C": 1, "G": 2, "T": 3}
    for line in open(vcfFile):
        # Skip any header lines
        if line[0] == "#":
            continue

        tokens = line.split("\t")
        ref = tokens[3]
        alt = tokens[4]
        infos = tokens[7]

        # Skip any lines with multiple non-ref alleles
        if "," in alt:
            continue

        # Get the QSS score from the info field
        matches = re.search("QSS=([0-9]+);", infos)
        QSS = matches.group(1)

        # Calculate the allele fraction as described
        # https://github.com/Illumina/strelka/issues/109
        tumourCounts = tokens[10].rstrip()
        t1Depth = tumourCounts.split(":")[0]
        tdata = tumourCounts.split(":")[4:]
        calledBaseT1Count = tdata[basePos[alt]].split(",")[0]

        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
        # DP:FDP:SDP:SUBDP:AU:CU:GU:TU
        print("\t".join([tokens[0], tokens[1], ref, alt, QSS, t1Depth, calledBaseT1Count]))


if __name__ == "__main__":
    main()
