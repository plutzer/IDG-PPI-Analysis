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
