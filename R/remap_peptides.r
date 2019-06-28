#' Remap and merge MaxQuant Evidence file into "tidy" data frame
#'
#' @description Given a protein database, this function will remap and annotate
#' phosphorylated peptide sequence obtained from a MaxQuant evidence.txt file
#' ready for further statistical analysis
#' @param evidence_file MaxQuant evidence.txt file
#' @param annotation_file Table with annotations of the columns to extract
#' @param fasta_file Protein fasta file used for original mapping with MaxQuant
#' @param window_size The size of the AA window to be used for the remapping
#' around the phosphosite. E.g. XXXSXXX, would correspond to a window size of 7.
#' Default = 15
#' @param min_prob The minimum probability observed for any phosphosite in a
#' peptide. Recommended to use 0.75, but assess with distribution as per 
#' vignette for example. Default = 0.75
#' @param max_sites This will filter out peptides with number of phosphosites 
#' greater than this value. Options are numeric values (e.g., 1, 2, 3) or "all".
#' E.g. 1 will only retain phosphopeptides with a single phosphosite, 2 with 2 
#' and "all" keeps all peptides. Default="all"
#' @param filter_site_method Two methods for filtering peptides by site
#' probability are available. "site" will remove all "phosphosites" identified
#' with a probability lower that "min_prob". "peptide" will remove all peptides
#' that contain a set combination of sites which do not contain a minimum 
#' probability, min_prob, and following the rules of "filter_peptide". Options 
#' are "site" or "peptide". Default = "peptide"
#' @param filter_peptide Multiple methods for filtering peptide probabilities. 
#' "method_6" will keep all phosphopeptides if there is evidence for existence 
#' in any "experiment" as per the "annotation_file" for any major site above 
#' 0.5, and at least one site is above the "min_prob". Currently only "method_6"
#' is implemented. Default = "method_6"
#' @param return_intensity If TRUE will evaluate the MaxQuant evidence results
#' and return a table of intensities for further analysis. If FALSE, the
#' evidence file will not be "tidied" and the intensity object returned by 
#' tidyEvidence will be NA. Default = TRUE
#' @param return_mapping_table If TRUE a mapping table will be generated from 
#' the fasta file input and returned. If FALSE the mapping_table returned by 
#' tidyEvidence will be NA. Default = TRUE
#' @param return_site_probability If TRUE will evaluate the MaxQuant evidence 
#' file and return a vector of all the site probabilities identified in all 
#' phosphorylated peptides. If FALSE, the site_probability vector returned by
#' tidyEvidence will be NA. Default = TRUE.
#' @param return_site_numbers If TRUE will evaluate the MaxQuant evidence 
#' file and return a vector of the number of sites identified for each peptide. 
#' If FALSE, the site_probability vector returned by tidyEvidence will be NA. 
#' Default = TRUE.
#' @param verbose Print progress to screen. Default = FALSE
#'
#' @examples
#'## Return the annotated data with extracted peptides:
#'## load example data in phosphoProcessR
#' data(human_fasta_example)
#' data(malaria_evidence_example)
#' data(malaria_annotation_example)
#'
#'## Take the first 100 for an example
#' evidence_file_example <- head(malaria_evidence_example, 100)
#'
#'## convert evidence file into useful table for further analyses
#' evidence_tidy <- tidyEvidence(evidence_file = evidence_file_example,
#'                               annotation_file = malaria_annotation_example,
#'                               fasta_file = human_fasta_example,
#'                               window_size = 15,
#'                               min_prob = 0.75,
#'                               filter_site_method = "peptide")
#'
#' ## View results of table ready for statistical analysis
#' head(evidence_tidy$intensity, 5)
#'
#' ## Access site probabilities
#' head(evidence_tidy$site_probability)
#'
#' @return List of 4 objects. 1) xx$intensity is a reformatted table of
#' intensities, now suitable for further statistical analyses, 2)
#' xx$mapping_table is the evidence file with additional annotations for the 
#' remapped peptides, 3) xx$site_probability is a vector of all the site 
#' probabilities identified (includes all major and minor sites) and 4) 
#' xx$site_numbers is a vector the number of sites that were identified 
#' (includes all major and minor sites)
#'
#' @export tidyEvidence
#'
#' @importFrom seqinr read.fasta
#' @import BiocParallel
#' @importFrom stringi stri_length
#' @importFrom stats median
#' @importFrom plyr join_all

tidyEvidence <- function(evidence_file = NULL,
                         annotation_file = NULL,
                         fasta_file = NULL,
                         window_size = 15,
                         min_prob = 0.75,
                         max_sites = "all",
                         filter_site_method = "peptide",
                         filter_peptide = "method_6",
                         return_intensity = TRUE,
                         return_mapping_table = TRUE,
                         return_site_probability = TRUE,
                         return_site_numbers = TRUE,
                         verbose = FALSE
                         ){
#----------------------------------------------
# format checks
if ((filter_site_method %in% c("site", "peptide")) == FALSE)
  stop("filter_site_method must be either: site [or] peptide")
if ((filter_peptide %in% c("method_6")) == FALSE)
  stop("probability_evidence must be a method either: method_6")
if ((min_prob >= 0 & min_prob <= 1) == FALSE)
  stop("min_prob must be between zero[0] and one[1]")
if((max_sites > 0 | max_sites =="all") == FALSE)
  stop("max_sites must be greater than zero[0] or \"all\"")
if (is.null(evidence_file)) 
  stop("evidence_file not provided; you must provide an evidence_file table")
# CHECK FORMAT OF TABLE HERE
if (is.null(fasta_file)) 
  stop("fasta_file not provided; you must provide an fasta_file table")
if (!is.list(fasta_file))
  stop("fasta_file is not a list format; something has gone wrong. Make sure
      fasta_file is formatted correctly and has been imported using
      read.fasta from seqinr as per the vignette.")
# check all seqs in fasta_file are the right class
class_check <- sapply(seq_len(length(fasta_file)), 
                      function(i)
                      class(fasta_file[[i]]) == "SeqFastaAA")
if(length(fasta_file) != length(class_check[class_check == TRUE]))
  stop("fasta_file contains AA sequences with incorrect class; something has
      gone wrong. Make sure fasta_file is formatted correctly and has been
      imported using read.fasta from seqinr as per the vignette.")
annot_check <- sapply(seq_len(length(fasta_file)), 
                      function(i)
                      "Annot" %in% names(attributes(fasta_file[[i]])))
if(length(fasta_file) != length(annot_check[annot_check == TRUE]))
  stop("fasta_file contains AA sequences without Annotations! Something has
    gone wrong. Make sure fasta_file is formatted correctly and has been
    imported using read.fasta from seqinr as per the vignette.")
# check format of annotation_file
if("samples" %in% colnames(annotation_file) == FALSE){
  stop("A annotation_file must be supplied with a column labelled 
       \"samples\"")
}
if("labels" %in% colnames(annotation_file) == FALSE){
  stop("A annotation_file must be supplied with a column labelled \"labels\"")
}
if("group" %in% colnames(annotation_file) == FALSE){
  stop("A annotation_file must be supplied with a column labelled \"group\"")
}
# cross reference experiment if provided in annotation_file, with evidence
if("experiment" %in% colnames(annotation_file) == FALSE && 
   "Experiment" %in% colnames(evidence_file) == FALSE){
  warning("Annotation_file does not contain a column labelled \"experiment\" &
          evidence_file does not contain a column labelled \"Experiment\". All
          measurements are assumed to have originated from the same 
          experiment.")
  # add experiment column to annotation_file for downstream compatability
  annotation_file$experiment <- as.factor(1)
  # force evidence file Experiment column to a single experiment
  evidence_file$Experiment <- as.factor(1)
  # still specify uniq_experiments for use later
  uniq_experiments <- unique(evidence_file$Experiment)
}
if("experiment" %in% colnames(annotation_file) == FALSE && 
   "Experiment" %in% colnames(evidence_file) == TRUE){
      uniq_experiments <- unique(evidence_file$Experiment)
      # determine number of experiments in file
      if(length(uniq_experiments) == 1){
        annotation_file$experiment <- uniq_experiments
        warning("Annotation_file did not contain a column labelled 
                \"experiment\". However, evidence_file did contain a column 
                labelled \"Experiment\" and has one experiment named ", 
                uniq_experiments, " analysis proceeding as a single 
                experiment")
      }
      if(length(uniq_experiments) > 1){
        stop("Annotation_file did not contain a column labelled \"experiment\".
              and more than one experiment was identified in the evidence_file. ",
              paste(uniq_experiments, collapse=","), " <== were found in the 
              annotation_file")
      }
}
if("experiment" %in% colnames(annotation_file) == TRUE){
  uniq_experiments <- unique(annotation_file$experiment)
  uniq_evidence_file_experiments <- unique(evidence_file$Experiment)
  # check that the Experiments exist:
  set_diff <- setdiff(uniq_experiments, uniq_evidence_file_experiments)
  if(length(set_diff) > 0)
    stop("Additional experiments have been found that differ between the 
         annotation_file and evidence_file. Check the files.\n",
         paste(uniq_experiments, collapse=","), " <== were found in the annotation_file &\n",
         paste(uniq_evidence_file_experiments, collapse=",") , " <== were found in the evidence_file")
} 

# check that the column names used in the annotation file for the samples
# correspond to column names in the evidence file input. Unique is used here
# as multiple experiments can be assigned to the same column names
cols_of_interest <- intersect(colnames(evidence_file), 
                              unique(annotation_file$samples))
set_diff <- setdiff(annotation_file$sample, colnames(evidence_file))
if(length(set_diff) > 0)
  stop("There is not a matching number of samples in the evidence_file from
        the samples in the annotation_file. Check file names match or
        correct annotation_file provided. ", set_diff, " not found!")
# check evidence file column names are the right format
cols_of_interest <- c("Proteins", "Sequence", "Modified.sequence",
                      "Phospho..STY..Probabilities", "Experiment")
set_diff <- setdiff(cols_of_interest, colnames(evidence_file))
if(length(set_diff) > 0)
  stop("Column names of evidence file differ to those expected. ",
       set_diff, " not found!. Check evidence file labels.")
# check that window_size is an odd number.
if (((window_size %% 2) == 0) == TRUE)
  stop("window_size must be an odd number! I.e. centered sequence of ++X++")
if (verbose) {
  message("window_size is ", window_size)
}
# window_size now becomes number of AAs to add to each side of site
window_size <- (window_size-1)/2

if(return_intensity == TRUE){
  return_mapping_table <- TRUE
  return_site_probability <- TRUE
  if (verbose) {
    message("return_intensity is TRUE. Over-riding all setting for 
            return_mapping_table and return_site_probability and setting both
            to TRUE")
  }
}

#----------------------------------------------
# END OF CHECKS
#----------------------------------------------
# set up number of analysis steps and counter
analysis_steps <- 6
analysis_start <- 0
# 1. annotate evidence file
evidence_annotate <- evidence_file
# 2. remove contaminants
# reverse peptides will contain "__", contaminants: "CON__", REV: "REV__"
cont_mapped <- grep("__", evidence_annotate$Proteins)
# 3. extract peptides with phospho-sites
contains_ph <- grep("(ph)", evidence_annotate$Modified.sequence)
# 4. remove unwanted rows from evidence file
keep <- setdiff(contains_ph, cont_mapped)
evidence_annotate <- evidence_annotate[keep,]
# 5. remove cases where there are no probabilities assigned
# UPDATE HERE - keep NAs
evidence_annotate$Phospho..STY..Probabilities[evidence_annotate$Phospho..STY..Probabilities == ""] <- NA
#evidence_annotate <- evidence_annotate[!is.na(evidence_annotate$Phospho..STY..Probabilities),]

# 6. determine probability distribution and site nunbers from all site probabilities
# includes minor and major probabilities
phospho_prob_out <- extract_probability(evidence_annotate$Phospho..STY..Probabilities)[[1]]
# needs updating - as the site numbers should be drawn from "evidence_annotate$Modified.sequence"
site_numbers_out <- extract_site_numbers(evidence_annotate$Phospho..STY..Probabilities)
# 7. extract data from the fasta file (accessible as table and list now)
fasta_data <- extract_fasta_data(fasta_file)

if (verbose) {
  analysis_start <- analysis_start+1
  message("[", analysis_start, "/", analysis_steps, "] Data extracted.")
}

if(return_mapping_table == TRUE){
  # this will differ with dimethyl option/accross experiment
  # "site" level filtering will remove sites < probability
  # "peptide" level filtering considers all sites
  if (filter_site_method == "site"){
    # filter out sites above a minimium, min_prob, probability
    # return to BPLAPPY
    filter_peptides <- bplapply(seq_along(evidence_annotate$Phospho..STY..Probabilities),
                                  function(i)
                                  filter_site_prob(input_sequence = evidence_annotate$Phospho..STY..Probabilities[i], 
                                                   threshold = min_prob,
                                                   method = "greater_equal"))
    # continue with the filtered phospho peptides
    evidence_annotate$prob_filter <- unlist(filter_peptides)
    # count number of (ph) sites in mod_seq - determines mono/di etc.
    phospho_no <- gsub("(ph)", "1/", evidence_annotate$prob_filter, fixed = TRUE)
    phospho_no <- gsub("([[:alpha:]])", "", phospho_no)
    phospho_no <- gsub("_|\\(|\\)", "", phospho_no)
    phospho_no <- sapply(strsplit(phospho_no, '/'), function(x) sum(as.numeric(x)))
    evidence_annotate$phospho_no <- phospho_no
    # because it is filtered at site level, the number of sites passing min_prob
    # is already contained in the phospho_no columns now.
    evidence_annotate$prob_count_filter <- prob_count_n
    # not interested in less than 0.5 for site level filtering
    evidence_annotate$prob_count_05 <- NA
  }

  # if site is not selected
  if (filter_site_method == "peptide"){
    # count site numbers
    phospho_no <- gsub("(ph)", "1/", evidence_annotate$Modified.sequence, fixed = TRUE)
    phospho_no <- gsub("([[:alpha:]])", "", phospho_no)
    phospho_no <- gsub("_|\\(|\\)", "", phospho_no)
    phospho_no <- sapply(strsplit(phospho_no, '/'), function(x) sum(as.numeric(x)))
    evidence_annotate$phospho_no <- phospho_no
    # determine number of p>="min_prob" sites (but filter later)
    prob_no <- gsub("([[:alpha:]])", "", evidence_annotate$Phospho..STY..Probabilities)
    # change BACK to bpapply? benchmark
    prob_count <- lapply(seq_along(prob_no), function(x)
                      regmatches(prob_no[x], 
                                 gregexpr("(?<=\\().*?(?=\\))", 
                                          prob_no[x], 
                                          perl=T))[[1]])
    prob_count_n <- sapply(seq_along(prob_count), function(x)
                          length(prob_count[[x]][prob_count[[x]] >= min_prob])
                          )
    # number of sites passing the min_prob
    evidence_annotate$prob_count_filter <- prob_count_n
    # Number of total sites (major and minor), that pass a minimum of (p > 0.5).
    # Used later for removal of multi-site assignments to a major with less a 
    # mimunum of 0.5.
    # change BACK to bpapply? benchmark
    mapping_table_ph_05 <- lapply(seq_along(evidence_annotate$Phospho..STY..Probabilities),
                                  function(i)
                                    filter_site_prob(input_sequence = evidence_annotate$Phospho..STY..Probabilities[i], 
                                                     threshold = 0.5,
                                                     method = "greater"))
    phospho_no <- gsub("(ph)", "1/", unlist(mapping_table_ph_05), fixed = TRUE)
    phospho_no <- gsub("([[:alpha:]])", "", phospho_no)
    phospho_no <- sapply(strsplit(phospho_no, '/'), function(x) sum(as.numeric(x)))
    evidence_annotate$prob_count_05 <- phospho_no
    # clean up the peptide of other modifications.
    # change BACK to bpapply? benchmark
    modifications <- lapply(seq_along(evidence_annotate$Modified.sequence), function(x)
                            regmatches(evidence_annotate$Modified.sequence[x], 
                                 gregexpr("(?<=\\().*?(?=\\))", 
                                          evidence_annotate$Modified.sequence[x], 
                                          perl=T))[[1]])
    
    modifications <- unique(unlist(modifications))
    modifications <- setdiff(modifications, "ph")
    modifications <- paste(modifications, collapse="|")
    eval_cmd <- paste("gsub(", "'", modifications, "'", ',"",', "evidence_annotate$Modified.sequence)", sep="")
    mod_seq <- eval(parse(text=eval_cmd))
    mod_seq <- gsub("\\()|_", "", mod_seq)
    evidence_annotate$prob_filter <- mod_seq
    }
}

if (verbose) {
  analysis_start <- analysis_start+1
  message("[", analysis_start, "/", analysis_steps, "] Peptides and probabilities extracted")
}

# obtain site positions and centered sequences for each site from fasta
# sites with no phospho will be empty
# this needs to be vectorised!
# library(stringi)
# CHANGE BACK TO BPLAPPY
centred_sites <- bplapply(seq_len(nrow(evidence_annotate)), function(x)
                          map_sites(evidence_annotate,
                                    fasta_file,
                                    fasta_data,
                                    window_size,
                                    which_row = x))

centred_sites <- unlist(centred_sites)
# add annotation to evidence file
evidence_annotate$annotation <- centred_sites

# return annotated evidence file - before filtering and merging
mapped_out_return <- evidence_annotate

# option to keep only sites with a defined number of phosphosites
if(max_sites != "all"){
  keep <- evidence_annotate$phospho_no <= max_sites
  evidence_annotate <- evidence_annotate[keep,]
}

if (verbose) {
  analysis_start <- analysis_start+1
  message("[", analysis_start, "/", analysis_steps, "] Centered sequences extracted. Site Annotation complete.")
}

if(return_intensity == TRUE){
  # experiment wide-filtering
  # filter out peptides without a minimum probability
  # onyl filters out peptides that dont match prob_count_filter and phospho_no 
  # for site level. I.e. all sites must pass
  # there is a bpapply in tidy_sample function.
  clean_evidence <- tidy_samples(evidence_file_in = evidence_annotate,
                                 annotation_file = annotation_file,
                                 filter_site_method = filter_site_method,
                                 min_prob = min_prob)
  
  if (verbose) {
    analysis_start <- analysis_start+1
    message("[", analysis_start, "/", analysis_steps, "] Experiment level filtering complete.")
  }
  
  # merge intensities by experiment
  experiment_lists <- lapply(seq_along(uniq_experiments), function(x)
                              clean_evidence[clean_evidence$experiment == 
                                             uniq_experiments[x],])
  names(experiment_lists) <- uniq_experiments

  # change to bplappy
  clean_evidence <- bplapply(seq_along(experiment_lists), function(x)
                             tidy_experiments(experiment_lists[[x]],
                                          annotation_file,
                                          names(experiment_lists)[x]))
  if (verbose) {
    analysis_start <- analysis_start+1
    message("[", analysis_start, "/", analysis_steps, "] Merging of intensity values complete.")
  }
  
  # merge experiments into one output, these already have unique names
  # check length to determine if merge is required.
  if(length(clean_evidence) == 1){
    # if just one experiment, only take the first object in the list
    med_out <- clean_evidence[[1]]
  }
  if(length(clean_evidence) > 1){
    # ideally query a database here, keys will be all mapped sites.
    med_out <- lapply(1:length(clean_evidence), function(x)
                      data.frame(clean_evidence[[x]],
                                 "ID"=as.character(
                                   rownames(data.frame(clean_evidence[[x]])))))
    #library(plyr)
    med_out <- join_all(med_out, by="ID", type="full")
    rownames(med_out) <- med_out$ID
    # remove ID column
    med_out <- med_out[ , !(colnames(med_out) %in% "ID")]
  }
}

if (verbose) {
  analysis_start <- analysis_start+1
  message("[", analysis_start, "/", analysis_steps, "] Data clean-up complete.")
}
  
# Define how returns are selected
# return_site_probability
if(return_intensity == FALSE){
  med_out <- NA
}
if(return_mapping_table == FALSE){
  mapped_out_return <- NA
}
if(return_site_probability == FALSE){
  phospho_prob_out <- NA
}
if(return_site_numbers == FALSE){
  site_numbers_out <- NA
}


return(list("intensity" = med_out,
       "mapping_table" = mapped_out_return,
       "site_probability" = phospho_prob_out,
       "site_numbers" = site_numbers_out))
}

# helper function for mapping the phospho-positions &
# extracting the recentred sequences
map_sites <- function(evidence_annotate,
                      fasta_file,
                      fasta_data,
                      window_size,
                      which_row){

  fasta_protein_name <- fasta_data[[6]]
  fasta_gene_name <- fasta_data[[7]]
  fasta_evidence <- fasta_data[[8]]
  fasta_database <- fasta_data[[4]]
  fasta_length <- fasta_data[[5]]
  fasta_sequences <- fasta_data[[3]]
    
  # 1. determine if is multimapped
  lead_p_mapped <- unlist(strsplit(as.character(evidence_annotate$Proteins[which_row]), ";"))
  # check that they are not the same protein!
  lead_p_mapped <- unique(lead_p_mapped)
  
  # multi - mapped? I.e. greater than one lead protein?
  if(length(lead_p_mapped) != 0L){
    # 1. select best proteins
    best_proteins <- best_protein(lead_p_mapped)
    # which fasta database rows match the proteins
    row_matches <- sapply(seq_along(best_proteins), function(x)
      which(fasta_protein_name == best_proteins[x]))      
    # 2. extract genes here:
    gene_symbols <- sapply(seq_along(row_matches), function(x)
      fasta_gene_name[row_matches[x]])
    # 3. Get vector for the best gene criteria
    best_gene_rows <- c(sapply(seq_along(unique(gene_symbols)), function(i)
             best_gene(unique(gene_symbols)[i],
                       row_matches[which(unique(gene_symbols)[i] == gene_symbols)],
                       fasta_evidence,
                       fasta_database, 
                       fasta_length)))
    best_gene_rows <- unlist(best_gene_rows)
    # update annotations
    protein_annotation <- paste(fasta_protein_name[best_gene_rows], collapse=":")
    gene_symbols  <- paste(fasta_gene_name[best_gene_rows], collapse=":") 
    protein_database  <- paste(fasta_database[best_gene_rows], collapse=":") 
    protein_evidence  <- paste(fasta_evidence[best_gene_rows], collapse=":") 
    protein_length  <- paste(fasta_length[best_gene_rows], collapse=":") 
    # obtain stripped peptide sequence
    mod_peptide <- as.character(evidence_annotate$prob_filter[which_row])
    # extract sites
    # update row_matches to the best matches
    row_matches <- best_gene_rows
    # map to all proteins - end up with length of 1 or > for multi-mapped
    phospho_site_list <- lapply(seq_along(row_matches), function(x)
      phosho_positions(mod_peptide = mod_peptide, 
                       seq = fasta_sequences[row_matches[x]],
                       window_size = window_size))
    # i.e. just one protein, with potentially multiple sites
    if(length(phospho_site_list) == 1){
      phospho_sites <- paste(phospho_site_list[[1]]$sites, collapse=";")
      centered_seqs <- paste(phospho_site_list[[1]]$seqs, collapse=";")
    }
    # i.e. there are multiple proteins - i.e. MULTI-MAPPED
    if(length(phospho_site_list) > 1){
      phospho_sites <- c(sapply(seq_len(length(phospho_site_list)), function(x)
        paste(phospho_site_list[[x]]$sites, collapse=";")
      ))
      phospho_sites <- paste(phospho_sites, collapse=":")
      centered_seqs <- c(sapply(seq_len(length(phospho_site_list)), function(x)
        paste(phospho_site_list[[x]]$seqs, collapse=";")
      ))
      centered_seqs <- paste(centered_seqs, collapse=":")
    }
  return(paste(gene_symbols, 
               protein_annotation, 
               phospho_sites, 
               centered_seqs, 
               protein_database, 
               protein_evidence, 
               protein_length, 
               sep="|")
         )
  }
  
  # 2. if it is not annotated? It will have a length of zero - return(NO_MATCH)
  if(length(lead_p_mapped) == 0L){
    return("NO_MATCH")
  }

}

best_gene <- function(each_gene, 
                      each_gene_rows,
                      fasta_evidence, 
                      fasta_database, 
                      fasta_length){
  # If the genes are the same - determine best ACCESSION for EACH unique gene
  # 1) select gene with best (smallest number) Protein evidence (PE)
  each_gene_evidence_in <- sapply(seq_along(each_gene_rows), function(x)
    fasta_evidence[each_gene_rows[x]])
  # update vector of row matches
  # this will only select the minimum...
  each_gene_rows <- each_gene_rows[each_gene_evidence_in >= min(each_gene_evidence_in)]
  # If database is the same - select "sp" over "tr"
  each_gene_database_in <- sapply(seq_along(each_gene_rows), function(x)
    fasta_database[each_gene_rows[x]])
  tr_or_sp <- "sp" %in% each_gene_database_in
  # update the gene rows
  if(tr_or_sp == TRUE){
    each_gene_rows <- each_gene_rows[which(each_gene_database_in == "sp")]
  }
  # If multiple proteins remain - select longest protein
  each_gene_length <- sapply(seq_along(each_gene_rows), function(x)
    fasta_length[each_gene_rows[x]])
  # remove genes above the minimum / update vector of row matches
  each_gene_rows <- each_gene_rows[!each_gene_length > min(each_gene_length)]
return(each_gene_rows)   
}


best_protein <- function(proteins_in){
  # vectorise this
  # determine any matches?
  # anything without a -X now gets a -0 - used for ranking variants
  variants <- proteins_in
  variants[which(!grepl("-", variants))] <- 
    paste(variants[which(!grepl("-", variants))], "-0", sep="")
  variants <- as.numeric(gsub(".*-", "", variants))
  # strip all the variants to define unique proteins
  proteins <- gsub("-.*", "", proteins_in)
  uniq_proteins <- unique(proteins)
  # determine the best for the full set of proteins and return:
  best_proteins <- c(sapply(seq_along(uniq_proteins), function(i)
    proteins_in[intersect(which(variants == 
                                  min(variants[which(uniq_proteins[i] == 
                                                       proteins)])), 
                          which(uniq_proteins[i] == proteins))]
  ))
return(best_proteins)
}

# helper function for extracting phospho-positions
phosho_positions <- function(mod_peptide, seq, window_size){
  count1_2 <- strsplit(seq, gsub("[a-z]|[(]|[)]|_", "", mod_peptide))
  #1. determine the site of the peptide in the protein.
  frag_lengths <- stri_length(unlist(count1_2))
  #2. site positions.
  phospho_sites <- mod_peptide
  #this why the underscore is used to preserve the tailing phosphosites.
  phospho_sites <- gsub("_", "", phospho_sites)
  phospho_count <- unlist(strsplit(phospho_sites, "ph"))
  stripped_count <- gsub("[a-z]|[(]|[)]|_", "", phospho_count)
  # this is the length of the first fragment (i.e. start distance):
  phospho_position <- frag_lengths[1]
  # this now returns the position of the phosphosites in the protein sequence:
  phosho_sites <- sapply(seq_len(length(stripped_count)-1), function(i)
    phospho_position+sum(stri_length(stripped_count[1:i]))
  )
  # then get the centered sequences...
  # end the protein by the window size to include overhanging
  seq_extend <- gsub(" ", "", paste(paste(rep("X", window_size), collapse=""),
                      seq,
                      paste(rep("X", window_size), collapse=""), collapse=""))
  centred_seqs <- sapply(seq_len(length(phosho_sites)), function(i)
                    substr(seq_extend, phosho_sites[i], phosho_sites[i]+(window_size*2))
  )
return(list("sites"=phosho_sites, "seqs"=centred_seqs))
}

# helper function for filtering at site level
# method = "greater" [OR] "greater_equal"
filter_site_prob <- function(input_sequence, 
                             threshold,
                             method = "greater_equal"){
  # conditional within the brackets...
  # put into a function for cleaning up all individual sites!
  input_seq <- strsplit(as.character(input_sequence),"\\)|\\(", perl=TRUE)[[1]]

  if(method == "greater_equal"){
      input_seq <- sapply(1:length((input_seq)), function(i)
        ifelse(input_seq[i] < threshold, "", input_seq[i]))
      #replace the site scores with a (ph) when meeting criteria
      output_seq <- paste(sapply(1:length((input_seq)), function(i)
        ifelse(input_seq[i] >= threshold & input_seq[i] <= 1,
               "(ph)", input_seq[i])), collapse="")
  }
  if(method == "greater"){
    input_seq <- sapply(1:length((input_seq)), function(i)
      ifelse(input_seq[i] <= threshold, "", input_seq[i]))
    #replace the site scores with a (ph) when meeting criteria
    output_seq <- paste(sapply(1:length((input_seq)), function(i)
      ifelse(input_seq[i] > threshold & input_seq[i] <= 1,
             "(ph)", input_seq[i])), collapse="")
  }

return(output_seq)
}


# extract probability of each site, given a vector of sequences
extract_probability <- function(phospho_input){
  #phospho_input <- evidence_annotate$Phospho..STY..Probabilities
  phospho_prob <- gsub("([[:alpha:]])", "", phospho_input)
  # phospho_out: output the probabilities of all phosphosites
  phospho_prob_out <- gsub(")(", ", ", phospho_prob, fixed=TRUE)
  phospho_prob_out <- gsub("\\)|\\(", "", phospho_prob_out)
  phospho_prob_out <- paste(phospho_prob_out, collapse = ", ")
  phospho_prob_out <- strsplit(phospho_prob_out, ", ")[[1]]
  phospho_prob_out <- suppressWarnings(as.numeric(phospho_prob_out))
  phospho_prob_out <- phospho_prob_out[!is.na(phospho_prob_out)]
  # from here, fit a mixture model to automatically determine which to keep?
  # continue clean -up
  phospho_prob <- gsub(")(", ";", phospho_prob, fixed=TRUE)
  phospho_prob <- gsub("\\)|\\(", "", phospho_prob)
  phospho_prob[phospho_prob==1] <- "1;0"#if only one site
  # for instances where there are no probabilities, 
  # make a zero (these are multi-mapped)!
  phospho_prob[phospho_prob==""] <- "0;0"
  return(list(phospho_prob_out, phospho_prob))  
}

extract_site_numbers <- function(phospho_probabilities){
  phospho_no <- gsub("([[:alnum:]])|\\.", "", phospho_probabilities, 
                     fixed = FALSE)
  phospho_no <- gsub("()", "1/", phospho_no, 
                     fixed = TRUE)
  phospho_no <- sapply(strsplit(phospho_no, '/'), 
                       function(x) sum(as.numeric(x)))
return(phospho_no)
}
#evidence_file_in = evidence_annotate

tidy_samples <- function(evidence_file_in,
                         annotation_file,
                         filter_site_method,
                         min_prob){
  # Select columns of interest from annotation file
  samples_selected <- evidence_file_in[colnames(evidence_file_in) %in% annotation_file$samples]
  # replace all zeros for NA's - i.e. where an experiment had no measurement
  samples_selected[samples_selected == 0] <- NA
  # trimmed table
  clean_data <- data.frame("site_id" = evidence_file_in$annotation,
                           samples_selected,
                           #"max_prob" = evidence_file_in$phospho_prob_max,
                           "phospho_no" = evidence_file_in$phospho_no,
                           "prob_count_filter" = evidence_file_in$prob_count_filter,
                           "prob_count_05" = evidence_file_in$prob_count_05,
                           #"prob" = evidence_file_in$Phospho..STY..Probabilities,
                           "prob_filter" = evidence_file_in$Phospho..STY..Probabilities,
                           "experiment" = evidence_file_in$Experiment)
  uniq_experiments <- unique(clean_data$experiment)
  # remove duplicates (usually some duplication here)
  # should not have duplication intensity measurements at all
  clean_data <- unique(clean_data)
  # ONLY keep data that is annotated
  clean_data <- clean_data[(!is.na(clean_data$site_id)) ,]
  
  # for site just remove all that dont pass the filter
  if (filter_site_method == "site"){
    # contain a phospho site about the filter
    # AND there are greater number of sites above filter
    keep <- clean_data$prob_count_filter != 0 & 
      (clean_data$prob_count_filter >= clean_data$phospho_no)
    matches_keep2 <- clean_data[keep,]
  }
  if (filter_site_method == "peptide"){
    # remove peptides that dont pass the prob_count greater than 0.5
    # put in option here to cut at a minimum probability...
    #clean_data <- clean_data[clean_data$prob_count_05 >= clean_data$phospho_no, ]
    # uniq_sites remaining    
    uniq_site <- as.character(unique(clean_data$site_id))
    # determine which uniq_sites to keep
    # search accross all experiments for a minimum level of evidence
    # change back to bplappy
    matches <- bplapply(seq_len(length(uniq_site)), function(i)
                      keep_peptide_test(data = clean_data,
                                        site = uniq_site[i]))
    matches <- unlist(matches)
    uniq_site_keep <- uniq_site[matches]
    matches_keep2 <- clean_data[clean_data$site_id %in% uniq_site_keep,]
  }
return(matches_keep2) 
}
#data = clean_data
#site = uniq_site[1]
# check that the evidence does exist in the dataset accross all experiments
# this way there is just one subset and check
keep_peptide_test <- function(data, site){
  data_subset <- subset(data, site_id==site)
  #clean_data <- clean_data[clean_data$prob_count_05 >= clean_data$phospho_no, ]
  # evidence of some peptide have all sites being greater than a value
  # and some peptide having all sites greater than 0.5 (from same peptide)
  # this will implement strategy 3!
  #keep <- max(data_subset$prob_count_filter) >= max(data_subset$phospho_no)
  # strategy 6:
  # must contain a single site > 0.75
  # AND
  # all sites must be > 0.5
  keep <- (data_subset$prob_count_filter >= 1) & (data_subset$prob_count_05 >= data_subset$phospho_no)
  keep <- "TRUE" %in% keep
  
return(keep)
}

tidy_experiments <- function(evidence_file_in,
                         annotation_file,
                         experiment){
  # could vectorise this
  # update column annotations
  annotation_file_filter <- annotation_file[annotation_file$experiment == experiment,]
  samples_selected <- evidence_file_in[colnames(evidence_file_in) %in% annotation_file$samples]
  colnames(samples_selected) <- annotation_file_filter$labels[match(annotation_file_filter$samples,
                                                             colnames(samples_selected))]
  samples_selected$site_id <- evidence_file_in$site_id
  # determine unique peptide annotations
  uniq_site <- unique(samples_selected$site_id)
  # extract site matches in a list    
  matches <- lapply(seq_along(uniq_site), function(i)
    subset(samples_selected, site_id==uniq_site[i]))
  # update colnames
  matches <- lapply(seq_along(matches), function(i)
    matches[[i]][colnames(matches[[i]]) %in% annotation_file_filter$labels]) 
  # determine the median of the matches and return
  med <- lapply(seq_along(matches), function(i)
      apply(log2(matches[[i]]), 2, median, na.rm = TRUE))
  # collapse table
  med2 <- do.call("rbind", med)
  rownames(med2) <- uniq_site
  # clean up NAs
  na_remove <- as.vector(apply(is.na(med2), 1, sum))
  # remove all NAs (where a row is only NAs)
  na_remove <- c(na_remove < ncol(med2))
  
  med_out <- data.frame(med2[na_remove,])
  colnames(med_out) <- colnames(med2)
  return(med_out)
}

# reformat fasta data to be vectorised
extract_fasta_data <- function(fasta_in){
  #fasta_in <- fasta_file
  # reformat header
  fasta_headers <- sapply(seq_len(length(fasta_in)), function(i)
    attributes(fasta_in[[i]])$Annot)
  # obtain sequence
  fasta_sequences <- sapply(seq_len(length(fasta_in)), function(i)
    fasta_in[[i]][1])
  # which database
  # sp/tr etc...
  fasta_database <- sapply(seq_len(length(fasta_headers)), function(i)
    gsub("\\|.*|>", "", fasta_headers[i]))
  # fill in non sp/tr here?
  # length of sequence
  fasta_length <- c(sapply(seq_len(length(fasta_sequences)), function(i)
    nchar(fasta_sequences[i])))
  # gene name
  fasta_gene_name <- sapply(seq_len(length(fasta_headers)), function(i)
    gsub("[ |>].*", "", gsub(".*GN=","",fasta_headers[i]))
  )
  # some dont have gene names...
  # replace NAs with NA
  fasta_gene_name[fasta_gene_name==""] <- "NA"
  
  fasta_evidence <- sapply(seq_len(length(fasta_headers)), function(i)
    gsub("[ |>].*", "", gsub(".*PE=","",fasta_headers[i]))
  )
  # some dont have evidence...
  # replace NAs with NA
  fasta_evidence[fasta_evidence==""] <- "NA"
  
  fasta_protein <- sapply(seq_along(names(fasta_in)), function(i)
    sub(".*\\|(.*)\\|.*","\\1",names(fasta_in[i]))
  )
  
  fasta_annotation <- data.frame("header" = fasta_headers,
                                 "seq" = fasta_sequences,
                                 "database" = fasta_database,
                                 "length" = fasta_length,
                                 "protein" = fasta_protein,
                                 "gene" = fasta_gene_name,
                                 "evidence" = fasta_evidence)
  
return(list(fasta_annotation,
            fasta_headers,
            fasta_sequences,
            fasta_database,
            fasta_length,
            fasta_protein,
            fasta_gene_name,
            fasta_evidence))
}
