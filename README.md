# vPro-MS
Untargeted viral proteomics workflow for the identification of human-pathogenic viruses

This script of the viral proteomics workflow (vPro-MS) enables identification of human-pathogenic viruses from patient samples by untargeted proteomics. vPro-MS is based on an in-silico derived peptide library covering the human virome in [UniProtKB](https://www.uniprot.org/) (331 viruses, 20,386 genomes, 121,977 peptides).  The script is intended to identify human-pathogenic viruses from DiaNN (https://github.com/vdemichev/DiaNN) outputs of either DIA or diaPASEF data. A scoring algorithm (vProID) assesses the confidence of virus identification and the results are finally summarized in a report table.

![vPro-MS workflow](./workflow-BioRender.png)

# Installation

# Library
The vPro Peptide Library can be downloaded from zenodo (https://zenodo.org/records/13832021). The zip folder contains 3 peptide FASTA files (Contaminants.fasta, Human.fasta, vPro.Virus.fasta), which were used to predict the spectral library (vPro-lib.predicted.speclib). Please note, that the additional commands “--cut” and “--duplicate-proteins” are needed to reprocess the prediction in DiaNN. This spectral library should be used to identify peptide sequences from samples of human origin using DiaNN. Furthermore, the folder contains the metadata file of the viral peptide sequences (vPro.Peptide.Library.txt) and a summary file of the virus taxonomy covered by the library (Taxonomy.Summary.txt). The metadata file is used by this vPro script to identify viruses from the DiaNN main report.

# Usage

# Output
The vPro script summarizes virus identification results in a single report table (Results_vPro.txt). The results represent an independent analysis of each sample in the DiaNN main report. The virus identifications are already filtered according to the thresholds configured in the vPro script. The report consists of the following columns:

    • Run:                                 Name of the sample 
    • Species:                             Virus species identified (NA = no identification) 	
    • vProID.Score:                        vProID score  
    • No.Peptide.Sequences:                Number of peptide sequences unique to the virus species	
    • Virus.Quantity:                      Quantity of the virus species calculated using the Top3 approach for absolute protein quantification	
    • Peptide.Sequences: 	               Viral peptide sequences 
    • CScores:                             CScores of each virus peptide as reported by DiaNN	
    • Subspecies:                          Virus subspecies identified	
    • No.Peptide.Sequences.Subspecies:     Number of peptide sequences unique to the virus subspecies	
    • Proteomes:                           Uniprot proteome accession number(s) for the top ranked virus proteome(s)

# Reference
Grossegesse, M.; Horn, F.; Kurth, A.; Lasch, P.; Nitsche, A.; Doellinger, J. vPro-MS enables identification of human-pathogenic viruses from patient samples by untargeted proteomics. medRxiv 2024, https://doi.org/10.1101/2024.08.21.24312107

# Support
Please post any questions, feedback, comments or suggestions on the GitHub Discussion board.
