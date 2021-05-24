# IDG Interactome PPI Analysis Pipeline
The Interactome PPI Analysis Pipeline uses two different protein interaction identifying software – 
SAINT and CompPASS – to form a consensus prediction on protein-protein interactions. The goal of this 
pipeline is to better identify protein-protein interactions and filter through contaminants in 
Affinity Purification Mass Spectrometry (AP/MS) data. It currently adds GO and BioGrid annotations.

#Installation
- MaxQuant
- Python 3.9
- R 4.0.3
- Copy precompiled binary file for Saint Express from “external” folder to “bin” folder
- Install cRomppass branch from the smarasolo github repository using:
```console
library("devtools")
devtools::install_github("smarasolo/cRomppass")
library("cRomppass") 
```
- Install org.Hs.eg.db through the BiocManger package
```console
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
```

#Input Files
experimentalDesign.csv
- 5 columns
  - Experiment name – the protein name and run number of the experiment
    - ** Has to match experiment name that was used in MaxQuant because it is matched with the MS/MS count column from MaxQuant for further analysis **
    - You can make sure the experiment name is the same by checking the experiment_name MS/MS count column in the MaxQuant results
- Type – T, for test/experiments, or C, for controls
- Bait – name of bait
- Replicate – replicate number
- Bait ID – Uniprot ID for Bait
  - Exactly matches the identifier in the MaxQuant results so if a bait matches two Protein IDS, it is recorded in the Bait Id

proteinGroups.txt

#Usage
Basic usage:
1) Open up a terminal (powershell, cmd, or cmder)
2) Navigate to the directory with your proteinGroups.txt and experimentDesign.csv file
```console
>cd "C:\Users\majorlab.3225-WD-00001\Box\CellBio-MajorLab\Users\Dhaval\MS data searches\IDG\MiniTurbo-TurboID\2019022"
```
3) Execute the program. By default it will look for a proteinGroups.txt file and ExperimentalDesign.csv file. If you changed the names, you can specifiy their path when calling the program:
```console
>score_APMS.exe --experimentalDesign ED_MA.csv
```

Full usage options:
```console
>score_APMS.exe --help
usage: score_APMS.exe [-h] [--proteinGroups PROTEINGROUPS] [--experimentalDesign EXPERIMENTALDESIGN] [--outputPath OUTPUTPATH] [-qs {spc,LFQ,intensity}]
                      [-qc {spc,LFQ,intensity}] [-s {v2,express,q}] [--nburnin NBURNIN] [--niter NITER] [--lowMode {0,1}] [--minFold {0,1}] [--normalize {0,1}]
                      [--normalize-control {0,1}] [--compress-n-ctrl COMPRESS_N_CTRL] [--compress-n-rep COMPRESS_N_REP]

This is the entry point to the program. It will execute the requested tasks.

optional arguments:
  -h, --help            show this help message and exit
  --proteinGroups PROTEINGROUPS
                        path to MaxQuant ProteinGroups.txt
  --experimentalDesign EXPERIMENTALDESIGN
                        path to experimental design file
  --outputPath OUTPUTPATH
                        path for the output directory. If it already exists it will be overwritten.
  -qs {spc,LFQ,intensity}, --quantification-saint {spc,LFQ,intensity}
                        quantification values to use for SAINT
  -qc {spc,LFQ,intensity}, --quantification-comppass {spc,LFQ,intensity}
                        quantification values to use for CompPASS
  -s {v2,express,q}, --SAINT {v2,express,q}
                        SAINT version to use
  --nburnin NBURNIN     SAINT v2: number of burn-in iterations in MCMC. Default: 2000
  --niter NITER         SAINT v2: number of main iterations in MCMC. Default: 10000
  --lowMode {0,1}       SAINT v2: exclude extremely high counts in the model. - If baits are densely connected or dataset is small (few baits), use 1. - otherwise, use 0.
                        Default: 0
  --minFold {0,1}       SAINT v2: forcing separation between true and false distributions. - If user wishes to allow typical contaminants with significant differential
                        enrichment over control purifications, use 0. - otherwise, use 1. Default: 0
  --normalize {0,1}     SAINT v2: divide the counts by the total spectral counts in each IP. Default: 1
  --normalize-control {0,1}
                        SAINTq: normalize control intensities by multiplying a constant to all control intensities so that the average observed test intensities is equal
                        to the average control intensities. Default: 0
  --compress-n-ctrl COMPRESS_N_CTRL
                        SAINTq: the number of control baits used in calculations, with priority for baits with greatest intensities. Setting this number to a large number
                        makes the program use all available control data (recommend in cases with at most several controls. Default: 1000
  --compress-n-rep COMPRESS_N_REP
                        SAINTq: the number of test bait replicates used for scoring, with priotiy given to the baits with higher probability scores. If this number is
                        greater than or equal to the number of available replicates, then the scores will use the data from all replicates. Otherwise, the highest scoring
                        replicate scores wil lbe averaged to yield the final probability score. Default: 1000
```

#Output files
Saint Input:
- bait.txt – list of baits
- interaction.txt – list of interactions
- prey.txt – list of prey

Saint Output:
- list.txt 
   - Output columns in list:
      - Spec
      - SpecSum
      - AvgSpec
      - NumReplicates
      - ctrlCounts – the number of spectra in the pry across all controls
      - AvgP – the average probability score for a bait and prey pair
      - MaxP – the largest probability score for a bait and prey pair across all replicate purifications
      - TopoAvgP
      - SaintScore - the probability the bait/prey interaction is a true interaction
      - logOddsScore
      - FoldChange
      - BFDR – A Bayesian false discovery rate SAINT calculates from a combined probability score from independent scoring of each replicate
      - Boosted_by

CompPASS Input:
- to_CompPASS.csv – formats experimentalDesign to fit input needed for compPASS

CompPASS Output:
- compPASS.csv – runs compPASS on data from MaxQuant
   - Columns in compPASS
      - AvePSM
      - SumAPSM
      - Mean
      - SD
      - Little.N – Number of baits(including control) that had at least 1 spectral count
      - Little.P – Number of replicates for this particular bait that had at least 1 spectral count
      - Z -  Z-score for each bait and prey pair
      - WD – WD-Score for each bait and prey pair as defined by Sowa et al.
      - Entropy – The Shannon Entropy for each spectral count

Final Output:
- Merge_CompPASS_Saint.csv – Merges the compPASS and Saint outputs
- Annotated_Merge_NO_FILTER.csv – Annotates the merged output with GO and Biogrid annotations, but no filter for Saint or CompPASS
- Annotated_Merge_Saint_filter.csv – Takes the merged and annotated file from before and filters Saint for BFDR <= 0.05 and AvgP >= 0.7
- Annotated_Merge_All_filtered.csv – Takes the merged, annotated, and filtered for Saint file and filters compPASS by taking the top 5% scored of each bait or the top 10 of each bait, if there are more than 10 in the top 5%
- [bait name].csv – Every bait will have one of these files that includes the Entrez ID of any prey-prey interaction seen in that bait

