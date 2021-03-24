#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups
#import tab_to_JSON

SAINT_DIR = os.path.dirname(os.path.realpath(__file__)) + "/build/"
SAINT_v2_INT_DIR = SAINT_DIR + "saint-int-ctrl"
SAINT_v2_SPC_DIR = SAINT_DIR + "saint-spc-ctrl"
SAINT_EXPRESS_INT_DIR = SAINT_DIR + "SAINTexpress-int"
SAINT_EXPRESS_SPC_DIR = SAINT_DIR + "SAINTexpress-spc"
SAINT_Q_DIR = SAINT_DIR + "SAINTq"

########################################################################################################################
# Command line argument parsing
########################################################################################################################

description = "This is the entry point to the program. It will execute the requested tasks."

# initialize the parser
parser = argparse.ArgumentParser(description=description)


# generic arguments
parser.add_argument("--proteinGroups",
                    help="path to MaxQuant ProteinGroups.txt",
                    default="proteinGroups.txt")

parser.add_argument("--experimentalDesign",
                    help="path to experimental design file",
                    default="ExperimentalDesign.csv")

# output arguments
parser.add_argument("--outputPath",
                    help="path for the output directory. If it already exists it will be overwritten.",
                    default=".")

parser.add_argument("-qs", "--quantification-saint",
                    help="quantification values to use for SAINT",
                    choices=["spc", "LFQ", "intensity"],
                    default="spc")

parser.add_argument("-qc", "--quantification-comppass",
                    help="quantification values to use for CompPASS",
                    choices=["spc", "LFQ", "intensity"],
                    default="spc")

#parser.add_argument("-i", "--imputation",
#                    help="method to impute missing quantification data",
#                    choices=["none", "perseus"],
#                    default="none")

# SAINT arguments
parser.add_argument("-s", "--SAINT",
                    help="SAINT version to use",
                    choices=["v2", "express", "q"],
                    default="express")

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

# SAINTq arguments
parser.add_argument("--normalize-control",
                    help="SAINTq: normalize control intensities by multiplying a constant to all control intensities\n" +
                         "so that the average observed test intensities is equal to the average control intensities.\n" +
                         " Default: 0",
                    choices=["0", "1"],
                    default="0")

parser.add_argument("--compress-n-ctrl",
                    help="SAINTq: the number of control baits used in calculations, with priority for baits with\n" +
                         "greatest intensities. Setting this number to a large number makes the program use all\n" +
                         "available control data (recommend in cases with at most several controls. Default: 1000",
                    type=int,
                    default=1000)

parser.add_argument("--compress-n-rep",
                    help="SAINTq: the number of test bait replicates used for scoring, with priotiy given to the\n" +
                         "baits with higher probability scores. If this number is greater than or equal to the number\n" +
                         "of available replicates, then the scores will use the data from all replicates. Otherwise, \n" +
                         "the highest scoring replicate scores wil lbe averaged to yield the final probability score.\n" +
                         "Default: 1000",
                    type=int,
                    default=1000)

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

if args.SAINT == "q":
    if not os.path.exists(SAINT_Q_DIR):
        print("ERROR: " + SAINT_Q_DIR + " not found")
        exit(1)

if args.SAINT == "v2":
    if args.quantification_saint == "intensity" or args.quantification_saint == "LFQ":
        if not os.path.exists(SAINT_v2_INT_DIR):
            print("ERROR: " + SAINT_v2_INT_DIR + " not found")
            exit(1)
    elif args.quantification_saint == "spc":
        if not os.path.exists(SAINT_v2_SPC_DIR):
            print("ERROR: " + SAINT_v2_SPC_DIR + " not found")
            exit(1)

if args.SAINT == "express":
    if args.quantification_saint == "intensity" or args.quantification_saint == "LFQ":
        if not os.path.exists(SAINT_EXPRESS_INT_DIR):
            print("ERROR: " + SAINT_EXPRESS_INT_DIR + " not found")
            exit(1)
    elif args.quantification_saint == "spc":
        if not os.path.exists(SAINT_EXPRESS_SPC_DIR):
            print("ERROR: " + SAINT_EXPRESS_SPC_DIR + " not found")
            exit(1)

if args.SAINT == "q":
    if args.quantification_saint == "spc":
        print("ERROR: SAINTq does not support spectral counts.")
        exit(1)

# confirm they aren't trying to impute spc
#if args.imputation and args.imputation == "1" and args.quantification == "spc":
#    print("ERROR: Imputation of spectral counts is not supported.")
#    exit(1)

# parse experimental design
experimental_design = ExperimentalDesign(args.experimentalDesign)

# process MaxQuant proteinGroups
protein_groups = ProteinGroups(experimental_design, args.proteinGroups,
                               args.quantification_saint, args.quantification_comppass)

# impute missing values if requested
#if args.imputation == "1":
#    protein_groups.impute()

# output SAINT formatted files
if args.SAINT == "q":
    protein_groups.to_SAINTq(args.outputPath, args.normalize_control, args.compress_n_ctrl, args.compress_n_rep)
else:
    protein_groups.to_SAINT(args.outputPath)

# output ComPASS formatted files
#prey2bait2comppass = protein_groups.calc_CompPASS()
protein_groups.to_CompPASS(args.outputPath)

# execute SAINT
if args.SAINT == "q":
    SAINT_path = SAINT_Q_DIR
    print("\nExecuting SAINT:", SAINT_path)
    p = subprocess.run([SAINT_path,
                        os.path.join(args.outputPath, "params.txt")])


if args.SAINT == "v2":
    if args.quantification_saint == "intensity" or args.quantification_saint == "LFQ":
        SAINT_path = SAINT_v2_INT_DIR
    else:
        SAINT_path = SAINT_v2_SPC_DIR

    print("\nExecuting SAINT:", SAINT_path)
    p = subprocess.run([SAINT_path,
                        os.path.join(args.outputPath, "output/interaction.txt"),
                        os.path.join(args.outputPath, "output/prey.txt"),
                        os.path.join(args.outputPath, "output/bait.txt"),
                        str(args.nburnin),
                        str(args.niter),
                        str(args.lowMode),
                        str(args.minFold),
                        str(args.normalize)])

if args.SAINT == "express":
    if args.quantification_saint == "intensity" or args.quantification_saint == "LFQ":
        SAINT_path = SAINT_EXPRESS_INT_DIR
    else:
        SAINT_path = SAINT_EXPRESS_SPC_DIR

    print("\nExecuting SAINT:", SAINT_path)
    p = subprocess.run([SAINT_path,
                        "-L",
                        "2",
                        os.path.join(args.outputPath, "interaction.txt"),
                        os.path.join(args.outputPath, "prey.txt"),
                        os.path.join(args.outputPath, "bait.txt")],
                       cwd="output\\")


#if args.SAINT == "express":
#    protein_groups.align_scores(os.path.join(args.outputPath,  "output\list.txt"),
#                                prey2bait2comppass,
#                                os.path.join(args.outputPath, "candidates.tsv"))
#else:
#    protein_groups.align_scores(os.path.join(args.outputPath, "RESULT", "unique_interactions"),
#                                prey2bait2comppass,
#                                os.path.join(args.outputPath, "candidates.tsv"))

#Run R Script for CompPASS
q = subprocess.run(["Rscript",
                    "compPASS.R",
                    os.path.join(args.outputPath,"output/to_CompPASS.csv"),
                    os.path.join(args.outputPath, "output/compPASS.csv")])

# Start R Script to merge CompPASS and SAINT
merge = subprocess.run(["Rscript",
                    "Merge_CompPASS_SAINT.R",
                    os.path.join(args.outputPath, "output/compPASS.csv"),
                        os.path.join(args.outputPath, "output/list.txt"),
                    os.path.join(args.outputPath, "output/Merge_CompPASS_SAINT.csv")])

# Run R Script to annotate the merged files
annotate = subprocess.run(["Rscript",
                    "annotate_filter.R",
                    os.path.join(args.outputPath, "output/Merge_CompPASS_SAINT.csv"),
                     os.path.join(args.outputPath, "output/Annotated_Merge_filtered_ALL.csv"),
                  os.path.join(args.outputPath, "output/Annotated_Merge_filtered_DKK.csv")])