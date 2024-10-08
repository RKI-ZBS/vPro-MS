# vPro-MS
Untargeted viral proteomics workflow for the identification of human-pathogenic viruses

# Overview
This script of the viral proteomics workflow (vPro-MS) enables identification of human-pathogenic viruses from patient samples by untargeted proteomics. vPro-MS is based on an in-silico derived peptide library covering the human virome in [UniProtKB](https://www.uniprot.org/) (331 viruses, 20,386 genomes, 121,977 peptides).  The script is intended to identify human-pathogenic viruses from DiaNN (https://github.com/vdemichev/DiaNN) outputs of either DIA or diaPASEF data. A scoring algorithm (vProID) assesses the confidence of virus identification and the results are finally summarized in a report table.

![vPro-MS workflow](./workflow-BioRender.png)
Created in BioRender. Grossegesse, M. (2024) BioRender.com/s15u220

# Installation
## System Requirements
### Hardware Requirements
The software was developed on tested on a virtual machine with: 2x AMD CPU 2,96GHz and 8 GB RAM running on a Windows 10 Enterprise operating system.
### Software Requirements
The software was developed and tested with
- R (version 4.3.0)
- R package tidyverse (version 2.0.0)
- R package xfun (version 0.4.0)

## Data Sources
### Library
The vPro Peptide Library can be downloaded from zenodo (https://zenodo.org/records/13832021). The zip folder contains 3 peptide FASTA files (Contaminants.fasta, Human.fasta, vPro.Virus.fasta), which were used to predict the spectral library (vPro-lib.predicted.speclib). Please note, that the additional commands “--cut” and “--duplicate-proteins” are needed to reprocess the prediction in DiaNN. This spectral library should be used to identify peptide sequences from samples of human origin using DiaNN. Furthermore, the folder contains the metadata file of the viral peptide sequences (vPro.Peptide.Library.txt) and a summary file of the virus taxonomy covered by the library (Taxonomy.Summary.txt). The metadata file is used by this vPro script to identify viruses from the DiaNN main report.

inkl. data format

### Output file Dia-NN
Generated using DIA-NN (https://github.com/vdemichev/DiaNN)


## Installation Guide
1. Install R and its packages
    Get and install R-version from https://www.r-project.org/
    ```
    CLI>R
    R>install.packages("tidyverse")
    R>install.packages("xfun")
    R>exit
    ```    
2. Get vPro-MS
    ```
    git clone https://github.com/RKI-ZBS/vPro-MS.git
    ```
3. Download necessary peptide libraries
    ```
    wget https://zenodo.org/records/13832021/files/vPro-MS%20Library%201.0.zip?download=1
    unzip vPro-MS%20Library%201.0.zip?download=1
    ```
  
## Demo und Tests
### Demo Input
### Demo Output
### Demo Runtime

# Usage
## Small Examples
```
library(vPro)
setwd("C:/Home/Users/ProjectA/assignSpecies")
assign_viral_species(file_peptides = "report_VirusID_Specificity.tsv",
                     file_virusDB = "Viral.Peptide.Library.txt",
                     file_export = "results_virusID.tsv",
                     nr_human_peptides = 591159,
                     fdr = 0.01,
                     min_pep_species = 2,
                     min_pep_subspecies = 2,
                     min_virIDScore = 2,
                     topn_precursor = 3,
                     filter_virIDScore = TRUE)
```
## Input and Arguments
Argument | Description | Example
--- | --- | ---
file_peptides | | report_VirusID_Specificity.tsv
file_virusDB | | Viral.Peptide.Library.txt
file_export | | results_virusID.tsv
nr_human_peptides | | 591159
fdr | | 0.01
min_pep_species | | 2
min_pep_subspecies | | 2
min_virIDScore | | 2
topn_precursor | | 3
filter_virIDScore | | TRUE

## Output
The vPro script summarizes virus identification results in a single report table (Results_vPro.txt). The results represent an independent analysis of each sample in the DiaNN main report. The virus identifications are already filtered according to the thresholds configured in the vPro script. The report consists of the following columns:

Column | Description | Example
--- | --- | ---
Run | Name of the sample | 
Species | Virus species identified (NA = no identification) | 
vProID.Score | vProID score | 
No.Peptide.Sequences | Number of peptide sequences unique to the virus species | 
Virus.Quantity | Quantity of the virus species calculated using the Top3 approach for absolute protein quantification | 
Peptide.Sequences | Viral peptide sequences | 
CScores | CScores of each virus peptide as reported by DiaNN | 
Subspecies | Virus subspecies identified | 
No.Peptide.Sequences.Subspecies | Number of peptide sequences unique to the virus subspecies | 
Proteomes | Uniprot proteome accession number(s) for the top ranked virus proteome(s) | 

# Reference
Grossegesse, M.; Horn, F.; Kurth, A.; Lasch, P.; Nitsche, A.; Doellinger, J. vPro-MS enables identification of human-pathogenic viruses from patient samples by untargeted proteomics. medRxiv 2024, https://doi.org/10.1101/2024.08.21.24312107

# License

# Support
Please post any questions, feedback, comments or suggestions on the GitHub Discussion board.
