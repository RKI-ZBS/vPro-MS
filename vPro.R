#packages#
library(tidyverse)
library(xfun, include.only = "normalize_path")

####Functions####

#' Load a DIA-NN report file
#'
#' @param f_name A path to the DIA-NN main report file.
#' @param min_cscore Filter entries below this CScore. (default 0.95)
#'
#' @return Detected Peptides (data frame).
#'
#' @examples
#' diann <- read_peptides(f_name = "diann-report.tsv", min_cscore = 0.9)
#' diann <- read_peptides("diann-report.tsv")
read_peptides <- function(f_name, min_cscore = 0.95) {
  pep <- read_tsv(file = f_name,
                   show_col_types = FALSE,
                   col_select = c("Run", "Stripped.Sequence",
                                  "Precursor.Normalised", "CScore")) %>%
    filter(CScore >= min_cscore) %>%
    rename("run" = "Run", "seqs" = "Stripped.Sequence",
           "precursorNormalised" = "Precursor.Normalised", "cscore" = "CScore")
  return(pep)
}

#' Load the vPro Peptide Library
#'
#' @param f_name A path to the vPro Peptide Library.
#'
#' @return vPro Peptide Library (data frame)
#'
#' @examples
#' virLib <- read_virus_db(f_name = "vPro.Peptide.Library.txt")
#' virLib <- read_virus_db("vPro.Peptide.Library.txt")
read_virus_db <- function(f_name) {
  pepLib <- read_tsv(file = f_name,
                     quote = "",
                     show_col_types = FALSE,
                     col_select = c("Species", "Subspecies", "Proteomes",
                                    "PeptideSequence")) %>%
    separate_rows(Proteomes, sep = "; ") %>%
    mutate(Proteomes = str_remove(Proteomes, ": Genome")) %>%
    distinct() %>%
    rename("species" = "Species", "subspecies" = "Subspecies",
           "proteomes" = "Proteomes", "pepSeqs" = "PeptideSequence")
  return(pepLib)
}

#' Annotate identified peptides with database information and
#' get information about the occurrences of the peptides in the
#' measurement and the peptide library. Calculate the virID score.
#'
#' @param pep measured peptides (data frame)
#' @param pep_lib vPro Viral Peptide Library (data frame)
#' @param nr_human_peps number of human peptides in Library
#' @param peptide_fdr false discovery rate that should be applied
#'
#' @return annotated identified peptides (data frame)
#'
#' @examples
#' annotatedPeptides <- annotate_peptides(diann, virLib)
annotate_peptides <- function(pep,
                              pepLib,
                              nr_human_peps = 591159,
                              peptide_fdr = 0.01) {
  #taxonomic and UniProt annotation
  annPep <- merge.data.frame(pep, pepLib, by.x = "seqs", by.y = "pepSeqs")

  #get occurrences
  ##number of peptide sequences per run
  countPeptides <- pep %>%
    group_by(run) %>%
    mutate(nrPeptides = n_distinct(seqs)) %>%
    select(run, nrPeptides) %>%
    distinct() %>%
    ungroup()
  annPep <- merge.data.frame(annPep, countPeptides, by = "run")

  ##number of identified peptides per DB entry
  virusDBCountProteomes <- pepLib %>%
    group_by(proteomes) %>%
    mutate(nrPeptidesProteome = n_distinct(pepSeqs)) %>%
    select(proteomes, nrPeptidesProteome) %>%
    distinct() %>%
    ungroup()
  annPep <- merge.data.frame(annPep, virusDBCountProteomes, by = "proteomes")

  #get occurences part 2
  ##number of distinct sequences per run and proteome
  annPep <- annPep %>%
    group_by(run, proteomes) %>%
    mutate(nrSpecificPeptides = n_distinct(seqs)) %>%
    ungroup()

  #calculate vProIDScore
  ##number of all peptides assigned to a species
  virusDBCountSpecies <- pepLib %>%
    group_by(species) %>%
    mutate(nrPeptidesSpecies = n_distinct(pepSeqs)) %>%
    select(species, nrPeptidesSpecies) %>%
    distinct() %>%
    ungroup()
  virusDBCountTotal <-
    sum(virusDBCountSpecies$nrPeptidesSpecies)
  if (virusDBCountTotal <= 0) {
    stop("No Peptides found.")
  }

  ##expected number of false positive peptides
  annPep <- annPep %>%
    mutate(predFalsePeptides = ((nrPeptides * peptide_fdr * 
                                   nrPeptidesProteome) /
                                  (nr_human_peps + virusDBCountTotal)))
  ##vProIDScore
  annPep <- annPep  %>%
    mutate(vProIDScore = log10(nrSpecificPeptides /
                                predFalsePeptides))
  #remove unneccessary data
  annPep <- annPep %>%
    select(-nrPeptidesProteome, -nrPeptides, -predFalsePeptides)
  return(annPep)
}

#' Assign species taxonomy to all measurements. Apply filters and aggregate
#' functions.
#'
#' @param annPep annotated measured peptides and their vProID Score (data frame)
#'
#' @return species information for each measurement (data frame)
#'
#' @examples
#' run2species <- assign_species(annotatedPeptides)
assign_species <- function(annPep,
                           min_peptides_per_species = 2,
                           min_vProIDScore = 2,
                           topn_precursor = 3,
                           filter_vProIDScore = TRUE) {
  #get number of peptides found for a species
  #used for output and ranking of different taxonomic assignments
  procPeptides <- annPep %>%
    group_by(run, species) %>%
    mutate(nrPepSeqs = n_distinct(seqs)) %>%
    ungroup()

  #remove duplicated peptides originating from different precursors
  #keep precursor with the highest quantity
  procPeptides <- procPeptides %>%
    select(-subspecies) %>%
    group_by(run, seqs) %>%
    arrange(desc(precursorNormalised)) %>%
    distinct(run, seqs, proteomes, .keep_all = TRUE) %>%
    ungroup()

  #get mean of  <topn_precursor> most abundant precursors and aggregate
  #entries for each protein per run
  procPeptides <- procPeptides %>%
    group_by(run, proteomes) %>%
    arrange(desc(precursorNormalised)) %>%
    slice_max(precursorNormalised, n = topn_precursor, with_ties = FALSE) %>%
    mutate(meanPrecursorTop = mean(precursorNormalised)) %>%
    mutate(vProIDScoreTop = max(vProIDScore)) %>%
    select(-precursorNormalised, -vProIDScore) %>%
    distinct() %>%
    ungroup()

  #minimum number of peptides per species
  #remove proteomes with low vProIDScore
  procPeptides <- procPeptides %>%
    filter(nrSpecificPeptides >= min_peptides_per_species) %>%
    select(-nrSpecificPeptides)
  if (filter_vProIDScore) {
    procPeptides <- procPeptides %>%
      filter(vProIDScoreTop >= min_vProIDScore)
  }

  #Aggregate proteome information (i.e. per run and per species)
  procPeptides <- procPeptides %>%
    group_by(run, species) %>%
    distinct(run, species, proteomes, .keep_all = TRUE) %>%
    arrange(desc(vProIDScoreTop)) %>%
    slice_max(vProIDScoreTop, n = 1, with_ties = TRUE) %>%
    mutate(proteomes = paste(proteomes, collapse = "; ")) %>%
    ungroup()
  
  #select top entry with highest vProIDScore
  procPeptides <- procPeptides %>%
    group_by(run, species) %>%
    arrange(desc(vProIDScoreTop)) %>%
    slice_max(vProIDScoreTop, n = 1, with_ties = FALSE) %>%
    distinct() %>%
    ungroup()

  #for each species per run: aggregate sequences and corresponding cscores
  allSeqs <- annPep %>%
    group_by(run, species) %>%
    distinct(run, species, seqs, .keep_all = TRUE) %>%
    arrange(seqs) %>%
    mutate(allSeqs = paste(seqs, collapse = ";")) %>%
    mutate(allCscores = paste(cscore, collapse = ";")) %>%
    ungroup() %>%
    select(run, species, allSeqs, allCscores) %>%
    distinct()
  procPeptides <- merge.data.frame(procPeptides,
                                   allSeqs,
                                   by = c("run", "species"))
  #remove columns that are not needed
  procPeptides <- procPeptides %>%
    select(-seqs, -cscore)

  #foreach species per run: get mean of top precursor quantity
  quantityPrecursor <- annPep %>%
    group_by(run, species) %>%
    distinct(run, species, precursorNormalised, .keep_all = TRUE) %>%
    arrange(desc(precursorNormalised)) %>%
    slice_max(precursorNormalised, n = 3, with_ties = FALSE) %>%
    mutate(meanPrecursorTop = mean(precursorNormalised)) %>%
    ungroup() %>%
    select(run, species, meanPrecursorTop) %>%
    distinct()
  procPeptides <- procPeptides %>%
    select(-meanPrecursorTop)
  procPeptides <- merge.data.frame(procPeptides,
                                   quantityPrecursor,
                                   by = c("run", "species"))

  return(procPeptides)
}

#' Assign subspecies taxonomy.
#'
#' @param annPep annotated measured peptides and their virID Score (data frame)
#' @param annPep_tax species information for each measurement (data frame)
#'
#' @return species and subspecies information for each measurement (data frame)
#'
#' @examples
#' run2species <- assign_subspecies(annotatedPeptides, run2species)
assign_subspecies <- function(ann_pep,
                              ann_pep_tax,
                              min_peptides_per_subspecies = 2) {
  subspec <- ann_pep %>%
    select(run, species, seqs, subspecies) %>%
    distinct()
  #remove all entries without subspecies information
  subspec <- subspec[!is.na(subspec$subspecies), ]
  #get number of proteins supporting subspecies assignment and filter for it
  subspec <- subspec %>%
    group_by(run, species, subspecies) %>%
    mutate(nrPeptideSeqsSubspecies = n_distinct(seqs)) %>%
    select(-seqs) %>%
    distinct() %>%
    filter(nrPeptideSeqsSubspecies >= min_peptides_per_subspecies)
  #merge with original dataset
  res <- merge(ann_pep_tax,
               subspec,
               by = c("run", "species"),
               all.x = TRUE)
  return(res)
}

#' Write results to a file.
#'
#' @param results Results including species information (data frame)
#' @param f_name File name for output (file path)
#'
#' @return None
#'
#' @examples
#' export_results(run2species, f_name = "results_vPro.tsv")
export_results <- function(results, f_name = file_export) {
  #only keep interesting columns and order them
  results <- results %>%
    select(run, species, vProIDScoreTop, nrPepSeqs, meanPrecursorTop,
           allSeqs, allCscores, subspecies, nrPeptideSeqsSubspecies,
           proteomes) %>%
    relocate(run, species, subspecies, nrPepSeqs, nrPeptideSeqsSubspecies,
             vProIDScoreTop, meanPrecursorTop, allSeqs, allCscores,  proteomes)
  #order by run und vProIDScore
  results <- results %>%
    group_by(run) %>%
    arrange(run, desc(vProIDScoreTop)) %>%
    ungroup()
  #rename columns
  colnames(results) <- c("Run", "Species", "Subspecies", "No.Peptide.Sequences",
                     "No.Peptide.Sequences.Subspecies", "vProID.Score",
                     "Virus.Quantity", "Peptide.Sequences",
                     "CScores", "Proteomes")

  write.table(results, file = f_name, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Results written to", f_name)
}

#' Assign taxonomy to each run using a virusID Score
#'
#' @param file_peptides DiaNN main output file containing the identified precursors
#' * tab-delimited text file
#' * header included
#' * needs to contain columns: Run, Stripped.Sequence, Precursor.Normalised,
#'      CScore
#' @param file_virusDB vPro Peptide Library
#' * tab-delimited text file
#' * header included
#' * needs to contain columns: Species, Subspecies, Proteomes, PeptideSequence
#' @param file_export Output file of vPro-MS Virus Identification
#' * tab-delimited text file
#' * header included
#' * contains columns: Run, Species, Subspecies, No.Peptide.Sequences, 
#'  No.Peptide.Sequences.Subspecies,	vProID.Score,
#' 	Virus.Quantity, Peptide.Sequences, CScores,	Proteomes
#' @param min_cScore minimal CScore threshold for precursor identification
#' @param nr_human_peptides Number of human peptides in the library
#' @param fdr FalseDiscoveryRate for peptide identification
#' @param min_pep_species minimal number of identified peptides necessary for a
#'    species assignment
#' @param min_pep_subspecies minimal number of identified peptides necessary for
#'    a subspecies assignment
#' @param min_vProIDScore minimal vProID score necessary for taxonomic assignment
#' @param topn_precursor number of most abundant precursor that are
#'    used for quantification
#' @param filter_vProIDScore should the results be filtered according to vProIDScore
#'    (Boolean)
#'
#' @return
#' @export
#'
#' @examples
#' assign_viral_species("diann_output.tsv", "vPro.Peptide.Library.tsv", "results_vPro.tsv")
assign_viral_species <- function(file_peptides,
                                 file_virusDB,
                                 file_export,
                                 min_cScore = 0.95,
                                 nr_human_peptides = 591159,
                                 fdr = 0.01,
                                 min_pep_species = 2,
                                 min_pep_subspecies = 2,
                                 min_vProIDScore = 2,
                                 topn_precursor = 3,
                                 filter_vProIDScore = TRUE) {

  #Assert script conditions
  stopifnot(file_test("-f", file_peptides))
  stopifnot(file_test("-f", file_virusDB))
  if (file_test("-f", file_export)) {
    stop(paste("Output file", file_export, "already exists."))
  }
  dirExport <- dirname(normalize_path(file_export))
  if (!file.access(dirExport, 2) == 0) {
    stop(paste("Output directory", dirExport, "is not writeable."))
  }

  peptides <- read_peptides(file_peptides, min_cScore)
  virusDB <- read_virus_db(file_virusDB)
  annPeptides <- annotate_peptides(peptides,
                                   virusDB,
                                   nr_human_peps = nr_human_peptides,
                                   peptide_fdr = fdr)
  runSpecies <- assign_species(annPeptides,
                               min_pep_species,
                               min_vProIDScore,
                               topn_precursor,
                               filter_vProIDScore)
  runSpecies <- assign_subspecies(annPeptides,
                                  runSpecies,
                                  min_pep_subspecies)
  export_results(runSpecies, file_export)
}

###MAIN###
setwd("C:/Users/Doellingerj/Desktop/VirusID/Virus Calling v4")
assign_viral_species(file_peptides = "report_VirusID_Specificity.tsv",
                     file_virusDB = "vPro.Peptide.Library.txt",
                     file_export = "results_vPro.tsv",
                     nr_human_peptides = 591159,
                     fdr = 0.01,
                     min_pep_species = 2,
                     min_pep_subspecies = 2,
                     min_vProIDScore = 2,
                     topn_precursor = 3,
                     filter_vProIDScore = TRUE)
