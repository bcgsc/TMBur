#!/usr/bin/env python
# Author:  Richard Corbett
# coding:utf-8
# Purpose: Run mutation spectrum code to create figures
# Created: 07/12/2020
# Requires: SigProfilerMatrixGenerator and SigProfilerPlotting
"""
This short program will run the mutation spectrum plotting routines available from https://github.com/AlexandrovLab/SigProfilerPlotting
"""
import argparse
import os
import sys

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen


def main():
    # Initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--vcf_file_folder", help="path to folder containing small variant VCF files"
    )
    parser.add_argument(
        "-n", "--name", help="string to associate with the files/figures", default="mut_spec"
    )
    parser.add_argument(
        "-r", "--reference", help="Must be one of GRCh38,GRCh37,mm10,mm9,etc.", default="GRCh37"
    )

    # get at the arguments
    args = parser.parse_args()
    if args.vcf_file_folder is None or not os.path.exists(args.vcf_file_folder):
        parser.print_help()
        sys.exit()

    # set up to use GRCh37
    # genInstall.install('GRCh37', rsync=False, bash=True)  # should be set up when container is built
    matGen.SigProfilerMatrixGeneratorFunc(
        args.name, args.reference, args.vcf_file_folder, plot=True
    )


# Required stuff to run a script with funtions in it
if __name__ == "__main__":
    main()
