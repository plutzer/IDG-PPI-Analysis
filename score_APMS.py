#!/usr/bin/env python
import sys
import os
import argparse
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups

SAINT_DIR = os.path.dirname(os.path.realpath(__file__)) + "/build/"
SAINT_v2_INT_DIR = SAINT_DIR + "saint-int-ctrl"
SAINT_v2_SPC_DIR = SAINT_DIR + "saint-spc-ctrl"
SAINT_EXPRESS_INT_DIR = SAINT_DIR + "SAINTexpress-int"
SAINT_EXPRESS_SPC_DIR = SAINT_DIR + "SAINTexpress-spc"

########################################################################################################################
# Command line argument parsing
########################################################################################################################

description = "This is the entry point to the program. It will execute the requested tasks."

# initialize the parser
parser = argparse.ArgumentParser(description=description)


# generic arguments
parser.add_argument("proteinGroups",
                    help="path to MaxQuant ProteinGroups.txt")

parser.add_argument("experimentalDesign",
                    help="path to experimental design file")

# output arguments
parser.add_argument("outputPath",
                    help="path for the output directory. If it already exists it will be overwritten.")

parser.add_argument("-q", "--quantification",
                    help="quantification values to use",
                    choices=["spc", "LFQ", "intensity"],
                    default="LFQ")

parser.add_argument("-i", "--imputation",
                    help="method to impute missing quantification data",
                    choices=["none", "perseus"],
                    default="none")

# SAINT arguments
parser.add_argument("-s", "--SAINT",
                    help="SAINT version to use",
                    choices=["v2", "express"])

parser.add_argument("--nburnin",
                    help="SAINT v2: number of burn-in iterations in MCMC. Default: 2000",
                    type=int,
                    default=2000)

parser.add_argument("--niter",
                    help="SAINT v2: number of main iterations in MCMC. Default: 10000",
                    type=int,
                    default=10000)

parser.add_argument("--lowMode",
                    help="SAINT v2: exclude extremely high counts in the model.\n" +
                         " - If baits are densely connected or dataset is small (few baits), use 1.\n" +
                         " - otherwise, use 0. Default: 0",
                    choices=["0", "1"],
                    default="0")

parser.add_argument("--minFold",
                    help="SAINT v2: forcing separation between true and false distributions. \n" +
                         " - If user wishes to allow typical contaminants with significant \n" +
                         " differential enrichment over control purifications, use 0. \n" +
                         "- otherwise, use 1. Default: 0",
                    choices=["0", "1"],
                    default="0")

parser.add_argument("--normalize",
                    help="SAINT v2: divide the counts by the total spectral counts in each IP. Default: 1",
                    choices=["0", "1"],
                    default="1")


args = parser.parse_args()


########################################################################################################################
# Check arguments
########################################################################################################################

# output directory
if os.path.exists(args.outputPath):
    print("WARNING: Output directory: '" + args.outputPath + "' already exists. Results will be overwritten.")
else:
    os.makedirs(args.outputPath)

# proteinGroups
if not os.path.exists(args.proteinGroups) or not os.path.isfile(args.proteinGroups):
    print("ERROR: proteinGroups file not found: " + args.proteinGroups)
    exit(1)

# experimentalDesign
if not os.path.exists(args.experimentalDesign) or not os.path.isfile(args.experimentalDesign):
    print("ERROR: experimental design file not found: " + args.experimentalDesign)
    exit(1)

if args.SAINT == "v2":
    if args.quantification == "intensity" or args.quantification == "LFQ":
        if not os.path.exists(SAINT_v2_INT_DIR):
            print("ERROR: " + SAINT_v2_INT_DIR + " not found")
            exit(1)
    elif args.quantification == "spc":
        if not os.path.exists(SAINT_v2_SPC_DIR):
            print("ERROR: " + SAINT_v2_SPC_DIR + " not found")
            exit(1)

if args.SAINT == "express":
    if args.quantification == "intensity" or args.quantification == "LFQ":
        if not os.path.exists(SAINT_EXPRESS_INT_DIR):
            print("ERROR: " + SAINT_EXPRESS_INT_DIR + " not found")
            exit(1)
    elif args.quantification == "spc":
        if not os.path.exists(SAINT_EXPRESS_SPC_DIR):
            print("ERROR: " + SAINT_EXPRESS_SPC_DIR + " not found")
            exit(1)


# confirm they aren't trying to impute spc
if args.imputation and args.imputation == "1" and args.quantification == "spc":
    print("ERROR: Imputation of spectral counts is not supported.")
    exit(1)

# parse experimental design
experimental_design = ExperimentalDesign(args.experimentalDesign)

# process MaxQuant proteinGroups
protein_groups = ProteinGroups(experimental_design, args.proteinGroups, args.quantification)

# impute missing values if requested
if args.imputation == "1":
    protein_groups.impute()

if args.SAINT == "v2":
    protein_groups.to_SAINTv2(args.outputPath)
