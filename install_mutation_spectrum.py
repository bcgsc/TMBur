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

from SigProfilerMatrixGenerator import install as genInstall

def main():
    # Initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r", "--reference", help="Must be one of GRCh38,GRCh37,mm10,mm9,etc.", default="GRCh37"
    )

    # get at the arguments
    args = parser.parse_args()

    # set up to use GRCh37
    genInstall.install(args.reference, rsync=False, bash=True)  # should be set up when container is built


# Required stuff to run a script with funtions in it
if __name__ == "__main__":
    main()
