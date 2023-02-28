#!/usr/bin/env python
# coding:utf-8
# Author:  Richard Corbett
# Purpose: Split the Strelka SNV file into smaller files split by af range.
# Created: 19/08/2020
import os
import sys
from optparse import OptionParser


def main():
    parser = OptionParser(" %prog -v passed.somatic.vcf ")
    parser.add_option(
        "-v", dest="vcfFile", action="store", type="string", help="path to strelka vcf file."
    )

    # get at the arguments
    options, _ = parser.parse_args()
    vcfFile = options.vcfFile
    if vcfFile is None:
        parser.print_help()
        sys.exit()

    # Check that our vcf file exists
    if not os.path.isfile(vcfFile):
        sys.stderr.write(vcfFile + " is not a valid file\n")
        sys.exit()

    # Loops over the lines and pull out the called allele counts
    basename = os.path.basename(vcfFile)
    af_0p0_0p05_file = open("%s_%s.vcf" % (basename, "af_0p0_0p05"), "w")
    af_0p05_0p1_file = open("%s_%s.vcf" % (basename, "af_0p05_0p1"), "w")
    af_0p1_0p15_file = open("%s_%s.vcf" % (basename, "af_0p1_0p15"), "w")
    af_0p15_0p2_file = open("%s_%s.vcf" % (basename, "af_0p15_0p2"), "w")
    af_0p2_0p25_file = open("%s_%s.vcf" % (basename, "af_0p2_0p25"), "w")
    af_0p25_0p3_file = open("%s_%s.vcf" % (basename, "af_0p25_0p3"), "w")
    af_0p3_0p35_file = open("%s_%s.vcf" % (basename, "af_0p3_0p35"), "w")
    af_0p35_0p4_file = open("%s_%s.vcf" % (basename, "af_0p35_0p4"), "w")
    af_0p4_0p45_file = open("%s_%s.vcf" % (basename, "af_0p4_0p45"), "w")
    af_0p45_0p5_file = open("%s_%s.vcf" % (basename, "af_0p45_0p5"), "w")
    af_0p5_0p55_file = open("%s_%s.vcf" % (basename, "af_0p5_0p55"), "w")
    af_0p55_0p6_file = open("%s_%s.vcf" % (basename, "af_0p55_0p6"), "w")
    af_0p6_0p65_file = open("%s_%s.vcf" % (basename, "af_0p6_0p65"), "w")
    af_0p65_0p7_file = open("%s_%s.vcf" % (basename, "af_0p65_0p7"), "w")
    af_0p7_0p75_file = open("%s_%s.vcf" % (basename, "af_0p7_0p75"), "w")
    af_0p75_0p8_file = open("%s_%s.vcf" % (basename, "af_0p75_0p8"), "w")
    af_0p8_0p85_file = open("%s_%s.vcf" % (basename, "af_0p8_0p85"), "w")
    af_0p85_0p9_file = open("%s_%s.vcf" % (basename, "af_0p85_0p9"), "w")
    af_0p9_0p95_file = open("%s_%s.vcf" % (basename, "af_0p9_0p95"), "w")
    af_0p95_1p0_file = open("%s_%s.vcf" % (basename, "af_0p95_1p0"), "w")

    for line in open(vcfFile):
        # write header lines to all the files
        if line[0] == "#":
            print >> af_0p0_0p05_file, line.rstrip()
            print >> af_0p05_0p1_file, line.rstrip()
            print >> af_0p1_0p15_file, line.rstrip()
            print >> af_0p15_0p2_file, line.rstrip()
            print >> af_0p2_0p25_file, line.rstrip()
            print >> af_0p25_0p3_file, line.rstrip()
            print >> af_0p3_0p35_file, line.rstrip()
            print >> af_0p35_0p4_file, line.rstrip()
            print >> af_0p4_0p45_file, line.rstrip()
            print >> af_0p45_0p5_file, line.rstrip()
            print >> af_0p5_0p55_file, line.rstrip()
            print >> af_0p55_0p6_file, line.rstrip()
            print >> af_0p6_0p65_file, line.rstrip()
            print >> af_0p65_0p7_file, line.rstrip()
            print >> af_0p7_0p75_file, line.rstrip()
            print >> af_0p75_0p8_file, line.rstrip()
            print >> af_0p8_0p85_file, line.rstrip()
            print >> af_0p85_0p9_file, line.rstrip()
            print >> af_0p9_0p95_file, line.rstrip()
            print >> af_0p95_1p0_file, line.rstrip()
        else:
            af = calculate_af(line)
            if af >= 0 and af < 0.05:
                print >> af_0p0_0p05_file, line.rstrip()
            elif af >= 0.05 and af < 0.1:
                print >> af_0p05_0p1_file, line.rstrip()
            elif af >= 0.1 and af < 0.15:
                print >> af_0p1_0p15_file, line.rstrip()
            elif af >= 0.15 and af < 0.2:
                print >> af_0p15_0p2_file, line.rstrip()
            elif af >= 0.2 and af < 0.25:
                print >> af_0p2_0p25_file, line.rstrip()
            elif af >= 0.25 and af < 0.3:
                print >> af_0p25_0p3_file, line.rstrip()
            elif af >= 0.3 and af < 0.35:
                print >> af_0p3_0p35_file, line.rstrip()
            elif af >= 0.35 and af < 0.4:
                print >> af_0p35_0p4_file, line.rstrip()
            elif af >= 0.4 and af < 0.45:
                print >> af_0p4_0p45_file, line.rstrip()
            elif af >= 0.45 and af < 0.5:
                print >> af_0p45_0p5_file, line.rstrip()
            elif af >= 0.5 and af < 0.55:
                print >> af_0p5_0p55_file, line.rstrip()
            elif af >= 0.55 and af < 0.6:
                print >> af_0p55_0p6_file, line.rstrip()
            elif af >= 0.6 and af < 0.65:
                print >> af_0p6_0p65_file, line.rstrip()
            elif af >= 0.65 and af < 0.7:
                print >> af_0p65_0p7_file, line.rstrip()
            elif af >= 0.7 and af < 0.75:
                print >> af_0p7_0p75_file, line.rstrip()
            elif af >= 0.75 and af < 0.8:
                print >> af_0p75_0p8_file, line.rstrip()
            elif af >= 0.8 and af < 0.85:
                print >> af_0p8_0p85_file, line.rstrip()
            elif af >= 0.85 and af < 0.9:
                print >> af_0p85_0p9_file, line.rstrip()
            elif af >= 0.9 and af < 0.95:
                print >> af_0p9_0p95_file, line.rstrip()
            elif af >= 0.95 and af <= 1.0:  # this last one includes 1.0
                print >> af_0p95_1p0_file, line.rstrip()


def calculate_af(line):
    """
    Based on the manual from Strelka, af can be calculated as:
    refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
    altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
    tier1RefCounts = First comma-delimited value from $refCounts
    tier1AltCounts = First comma-delimited value from $altCounts
    Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    """
    base_pos = {"A": 0, "C": 1, "G": 2, "T": 3}

    tokens = line.split("\t")
    ref_base = tokens[3]
    alt_base = tokens[4]

    # counts for the tumour bases are in the 10th field
    tumour_counts = tokens[10].rstrip()
    base_counts = tumour_counts.split(":")[4:]
    alt_t1_count = float(base_counts[base_pos[alt_base]].split(",")[0])
    ref_t1_count = float(base_counts[base_pos[ref_base]].split(",")[0])
    af = alt_t1_count / (alt_t1_count + ref_t1_count)

    return af


if __name__ == "__main__":
    main()
