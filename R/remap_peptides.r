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
#' @param is_TMT TMT MaxQuant output does not include an experiment column, as
#' each experiment is a TMT label. Use "FALSE" if a di/tri methyl experiment,
#' which will extract experiment information from the "Experiment" column in the
#' evidence.txt MaxQuant table. Default = TRUE
#' @param min_prob The minimum probability observed for any phosphosite in a
#' peptide. Recommended to use 0.7, but first assess with distribution as per
#' output table. See vignette for example. Default = 0
#' @param filter_site_method Two methods for filtering peptides by site
#' probability are available. "site" will remove all phosphosites identified
#' with a probability lower that "min_prob". "peptide" will remove all peptides
#' that are predicted to contain a set combination of sites which do not conatin
#' a minimum probability of a site. Options are "site" or "peptide". Default =
#' "peptide"
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
#'                               is_TMT = TRUE,
#'                               min_prob = 0.7,
#'                               filter_site_method = "site")
#'
#' ## View results of table ready for statistical analysis
#' head(evidence_tidy$intensity, 5)
#'
#' ## Access site probabilities
#' head(evidence_tidy$site_probability)
#'
#' @return List of 3 objects. 1) xx$intensity is a reformatted table of
#' intensities, now suitable for further statistical analyses, 2)
#' xx$intensity is table of remapped peptides and 3) xx$site_probability is a
#' vector of all the site probabilities identified.
#'
#' @export tidyEvidence
#'
#' @importFrom seqinr read.fasta
#' @import BiocParallel
#' @importFrom stringi stri_length
#' @importFrom stats median

tidyEvidence <- function(evidence_file = NULL,
                         annotation_file = NULL,
                         fasta_file = NULL,
                         window_size = 15,
                         is_TMT = TRUE,
                         min_prob = 0,
                         filter_site_method = "peptide",
                         verbose = FALSE
                         ){
  #----------------------------------------------
  #format checks:

  if ((filter_site_method %in% c("site", "peptide")) == FALSE)
    stop("filter_site_method must be either site or peptide")
  if ((is_TMT %in% c("TRUE", "FALSE")) == FALSE)
    stop("is_TMT must be either TRUE or FALSE")
  if ((min_prob >= 0 & min_prob <= 1) == FALSE)
    stop("min_prob must be between zero[0] and one[1]")
  if (is.null(evidence_file)) stop("evidence_file not provided; you must provide
                                   an evidence_file table")
  if (is.null(fasta_file)) stop("fasta_file not provided; you must provide
                                an fasta_file table")
  if (!is.list(fasta_file))
    stop(
      "fasta_file is not a list format; something has gone wrong. Make sure
        fasta_file is formatted correctly and has been imported using
        read.fasta from seqinr as per the vignette.")

  # check all the seqs are the right class
  class_check <- sapply(seq_len(length(fasta_file)), function(i)
                        class(fasta_file[[i]]) == "SeqFastaAA")
  if(length(fasta_file) != length(class_check[class_check == TRUE]))
    stop(
      "fasta_file contains AA sequences with incorrect classes; something has
      gone wrong. Make sure fasta_file is formatted correctly and has been
      imported using read.fasta from seqinr as per the vignette.")

  annot_check <- sapply(seq_len(length(fasta_file)), function(i)
    "Annot" %in% names(attributes(fasta_file[[i]])))

  if(length(fasta_file) != length(annot_check[annot_check == TRUE]))
    stop(
      "fasta_file contains AA sequences without Annotations!; something has
      gone wrong. Make sure fasta_file is formatted correctly and has been
      imported using read.fasta from seqinr as per the vignette.")

  # check all formats - including column names.
  if("samples" %in% colnames(annotation_file) == FALSE){
    stop("A annotation_file must be supplied with a column labelled \"samples\"")
  }
  if("labels" %in% colnames(annotation_file) == FALSE){
    stop("A annotation_file must be supplied with a column labelled \"labels\"")
  }
  if("group" %in% colnames(annotation_file) == FALSE){
    stop("A annotation_file must be supplied with a column labelled \"group\"")
  }

  #check that the annotation_file$sample exist!
  cols_of_interest <- intersect(colnames(evidence_file), annotation_file$sample)
  set_diff <- setdiff(annotation_file$sample, colnames(evidence_file))

  if(length(cols_of_interest) != nrow(annotation_file))
    stop("There is not a matching number of samples in the evidence_file from
            the samples in the annotation_file. Check file names match or
            correct annotation_file provided. ", set_diff, " not found!")

  # check evidence file column names!
  cols_of_interest <- c("Proteins", "Sequence", "Modified.sequence",
                        "Phospho..STY..Probabilities")
  set_diff <- setdiff(cols_of_interest, colnames(evidence_file))

  if(length(set_diff) > 0)
    stop("Column names of evidence file differ to those expected. ",
         set_diff, " not found!")

  # check that window_size is an odd number.
  if (((window_size %% 2) == 0) == TRUE)
    stop("window_size must be an odd number! I.e. centered sequence
         of ++X++")
  if (verbose) {
    message("window_size is ", window_size)
  }

  # window.size now becomes number of AAs to add to each side of site
  window_size <- (window_size-1)/2

  #----------------------------------------------
  # 1. Build mapping table from evidence file
  # re-mapping and centering is performed on the "mod_seq" column
  mapping_table <- data.frame("protein" = evidence_file$Proteins,
                              "seq" = evidence_file$Sequence,
                              "mod_seq" = evidence_file$Modified.sequence,
                              "prob" = evidence_file$Phospho..STY..Probabilities
                              )

  mapping_table <- unique(mapping_table)

  # 2. extract peptides with phospho-site
  contains_ph <- grep("(ph)", mapping_table$mod_seq)
  # strip all the characters and lower case letters into a new column
  mapping_table$stripped <- gsub("[a-z]|[(]|[)]|_", "", mapping_table$mod_seq)

  # 3. identify contaminants
  # reverse peptides, will contain "__", contaminants: "CON__", REV: "REV__"
  cont_mapped <- grep("__", mapping_table$protein)

  # 4. remove unwanted rows from phospho mapped file
  keep <- setdiff(contains_ph, cont_mapped)
  mapping_table_ph <- mapping_table[keep,]

  # this will differ with dimethyl option/accross experiment
  # site level will remove sites < probability prior to merging
  if (filter_site_method == "site"){
    # filter out sites with a certain probability here
    # CANNOT REPLACE THE MODIFIED SEQUENCE - USED FOR JOIN LATER!
    mapping_table_ph$prob_filter <- sapply(1:length(mapping_table_ph$prob),
                                           function(i)
                                filter_site_prob(mapping_table_ph$prob[i], min_prob))
  }

  # also want to simply get all the probababilities!


  if (filter_site_method != "site"){
    mapping_table_ph$prob_filter <- mapping_table_ph$mod_seq
  }


  # 5. Extract some information from the fasta file
  fasta_headers <- sapply(seq_len(length(fasta_file)), function(i)
                          attributes(fasta_file[[i]])$Annot)
  fasta_sequences <- sapply(seq_len(length(fasta_file)), function(i)
                            fasta_file[[i]][1])
  # very slow:
  #library(stringi) #for stri_length function below:
  #sequence.lengths <- c(sapply(seq_len(length(fasta.sequences)), function(i)
  #  stri_length(fasta.sequences)))
  #test <- attributes(fasta.file[[1]])$Annot
  fasta_gene_name <- sapply(seq_len(length(fasta_headers)), function(i)
                      gsub("[ |>].*", "", gsub(".*GN=","",fasta_headers[i]))
                      )
  #fasta.protein.evidence <- sapply(seq_len(length(fasta.headers)), function(i)
  #  gsub("[ |>].*", "", gsub(".*PE=","",fasta.headers[i]))
  #)
  #fasta.proteins <- data.frame(Reduce(rbind,
  #strsplit(names(fasta.file), "|", fixed = TRUE)))

  # 6. Get centered position and centered sequence for each site from the database
  # this needs to be vectorised!
  centred_sites <- bplapply(seq_len(nrow(mapping_table_ph)), function(x)
                            map_sites(mapping_table_ph,
                                      fasta_file,
                                      fasta_gene_name,
                                      fasta_sequences,
                                      window_size,
                                      which_row = x))

  centred_sites <- unlist(centred_sites)

  # check match here.
  # i.e. or new data.frame may not match up if contains funky ids

  # 7. merge with evidence file annotation

  mapped_out <- data.frame("annotation"= centred_sites, mapping_table_ph)


  # X. - here is where to add the best annotation option.

  # 8. remerge with the evidence file:
  mapped_out <- data.frame("annotation" = mapped_out$annotation,
                           "Modified.sequence" = mapped_out$mod_seq)

  # Phospho..STY..Probabilities
  # merge in new evidence file:
  evidence_remap <- merge(mapped_out, evidence_file, by="Modified.sequence")



  # peptide level will not alter the phosphosite (but only remove collections
  # with less than minimum phosphosite probability)
    # filter out peptides without a minimum probability
    # 9. determine site probabilities
    phospho_prob <- gsub("([[:alpha:]])", "", evidence_remap$Phospho..STY..Probabilities)

    # add option here for phospho_out to output the probabilities of
    # all the phosphosites
    phospho_prob_out <- gsub(")(", ", ", phospho_prob, fixed=TRUE)
    phospho_prob_out <- gsub("\\)|\\(", "", phospho_prob_out)
    phospho_prob_out <- paste(phospho_prob_out, collapse = ", ")
    phospho_prob_out <- strsplit(phospho_prob_out, ", ")[[1]]
    phospho_prob_out <- as.numeric(phospho_prob_out)
    phospho_prob_out <- phospho_prob_out[!is.na(phospho_prob_out)]
    # from here, fit a mixture model to automatically determine which to keep.
    #

    # continue clean -up
    phospho_prob <- gsub(")(", ";", phospho_prob, fixed=TRUE)
    phospho_prob <- gsub("\\)|\\(", "", phospho_prob)
    phospho_prob[phospho_prob==1] <- "1;0"#if only one site
    phospho_prob[phospho_prob==""] <- "0;0"#for instances where there are no probabilities, make a zero (these are multi-mapped)!

    # 10. determine number of sites
    #phospho_no <- evidence_remap$Modified.sequence
    #phospho_no <- gsub("(ph)", "1", phospho_no, fixed = TRUE)
    #phospho_no <- gsub("([[:alpha:]])", "", phospho_no)
    #phospho_no <- gsub("\\)|\\(|\\_", "", phospho_no)
    #n_p <- as.numeric(phospho_no)
    # this is still in 1, 11, 111 format - clean up

    #take the max of peptide phospho scores and determine length

    # HERE IS WHERE MAX OR CLEAN BY MIN IS USED
    phospho_prob_max <- as.numeric(lapply(strsplit(phospho_prob,";"), max))
    #phospho_prob_length <- as.numeric(lapply(strsplit(phospho_prob,";"), length))

  #}



  # 11. Trim file and converge

  #REQUIRES COLUMNS SELECTION HERE:
  #Below is ok for TMT - ten-plex
  #note that there has been an update/difference with MaxQuant versions
  #if Reporter.intensity.not.corrected.0 - use Reporter.intensity.0
  #if Reporter.intensity.corrected.0 - use Reporter.intensity.corrected.0

  # Select columns of interest from annotation file
  samples_selected <- evidence_remap[,colnames(evidence_remap) %in% annotation_file$samples]
  colnames(samples_selected) <- annotation_file$labels

  # new trimmed table to work with
  clean_data <- data.frame("site_id" = evidence_remap$annotation,
                            samples_selected,
                            "max_prob" = phospho_prob_max)

  # For compatability with TMT.
  # CAN DO A CHECK HERE FOR THE COLUMN - IF EXIST?
  if(is_TMT == TRUE){
    clean_data$Experiment="1"
  }
  # remove duplicates (usually some duplication here)
  clean_data <- unique(clean_data)
  # replace all zeros for NA's - i.e. where an experiment had no measurement
  clean_data[clean_data == 0] <- NA
  # ONLY keep data that is annotated and complete
  # clean_data <- clean_data[(!is.na(clean_data$max_prob)) & !is.na(clean_data$site_id) ,]
  clean_data <- clean_data[(!is.na(clean_data$site_id)) ,]

  # 12. merge duplicate site values
  # key points:
    #1: identify common phoshposites within each experiment
    #2: collapse into the median value given all measurements recorded
    #3: if a phosphosite does not have an identification of >0.75 remove.

  # identify the unique sites/site combinations
  uniq_site <- as.character(unique(clean_data$site_id))

  # built a list of all the unique sites
  matches <- lapply(seq_len(length(uniq_site)), function(i)
                    subset(clean_data, site_id==uniq_site[i]))
  names(matches) <- uniq_site
  match_names <- names(matches)

  # keep peptides with a probability greater than min defined
  if (filter_site_method == "peptide"){
    keep <- sapply(seq_len(length(uniq_site)), function(i)
                  max(matches[[i]]$max_prob, na.rm = TRUE) >= min_prob)
    matches <- matches[keep]
    match_names <- names(matches)
  }


  # 13. Tidy up the list columns - reformat to sample values (intensity)
  matches <- lapply(seq_len(length(matches)), function(i)
                    matches[[i]][, colnames(matches[[i]]) %in% annotation_file$labels])
  names(matches) <- match_names

  # 14. merge the values, when there are mulitple values:
  # ADD IF HERE BASED ON METHOD (currently only median utilised)
  med <- lapply(seq_len(length(matches)), function(i)
                apply(log2(matches[[i]]), 2, median)) #,na.rm = TRUE))
  # 15. collapse to meaninful table:
  med2 <- do.call("rbind", med)
  rownames(med2) <- match_names
  # clean up NA's
  na_remove <- as.vector(apply(is.na(med2), 1, sum))
  # remove all columns are NA
  na_remove <- c(na_remove < ncol(med2))
  med_out <- med2[na_remove,]

return(list("intensity" = med_out,
       "mapping_table" = mapped_out,
       "site_probability"=phospho_prob_out))
}

# helper function for mapping the phospho-positions &
# extracting the recentred sequences
map_sites <- function(mapping_table_ph,
                      fasta_file,
                      fasta_gene_name,
                      fasta_sequences,
                      window_size,
                      which_row){
  #print(which.row)
  #which.row <- 317

  # 1. determine if is multimapped
  lead_p_mapped <- unlist(strsplit(as.character(mapping_table_ph$protein[which_row]), ";"))

  # check is there are mapped proteins (if there are - get data)
  if(length(lead_p_mapped) != 0L){
    # check that they are not the same protein!
    lead_p_mapped <- unique(lead_p_mapped)
    # collapse the protein annotation:
    protein_annotation <- paste(lead_p_mapped, collapse=":")
    # obtain the stripped peptide sequence:
    peptide <- as.character(mapping_table_ph$stripped[which_row])
    mod_peptide <- as.character(mapping_table_ph$prob_filter[which_row])

    # which fasta database rows match the protein:
    row_matches <- sapply(seq_len(length(lead_p_mapped)), function(x)
      grep(paste(lead_p_mapped[x], "[|]", sep=""), names(fasta_file)))

    # extract the genes here:
    gene_symbols <- sapply(seq_len(length(row_matches)), function(x)
                          fasta_gene_name[row_matches[x]])
    gene_symbols <- paste(gene_symbols, collapse=":")
    # extract the sites:
    #x=1
    phospho_site_list <- lapply(seq_len(length(row_matches)), function(x)
      phosho_positions(peptide, mod_peptide, seq = fasta_sequences[row_matches[x]],
                       window_size = window_size))
    #this is a list now


    # i.e. just one protein, with potentially multiple sites
    if(length(phospho_site_list) == 1){
      phospho_sites <- paste(phospho_site_list[[1]]$sites, collapse=";")
      centered_seqs <- paste(phospho_site_list[[1]]$seqs, collapse=";")

    }
    # i.e. there are multiple proteins
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
    #"A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN"
    return(paste(gene_symbols, protein_annotation, phospho_sites, centered_seqs, sep="|"))
  }

  # 2. if it is not annotated? It will have a length of zero - retunr(NO_MATCH)
  if(length(lead_p_mapped) == 0L){
    #lead_p_mapped <- "NO_MATCH"
    return("NO_MATCH")
  }

}


# helper function for extracting phospho-positions
phosho_positions <- function(peptide, mod_peptide, seq, window_size){
  #print(peptide)
  #str(seq)
  #print(seq)
  #print(mod_peptide)
  count1_2 <- strsplit(seq, peptide)
  #print(count1_2)
  #1. determine the site of the peptide in the protein.
  frag_lengths <- stri_length(unlist(count1_2))
  #print(frag_lengths)
  #2. site positions.
  #DOES IT CONTAIN PHOSPHO SITES and determine the sites:
  phospho_sites <- mod_peptide
  #this why the underscore is used to preserve the tailing phosphosites.
  phospho_sites <- gsub("_", "", phospho_sites)
  phospho_count <- unlist(strsplit(phospho_sites, "ph"))
  stripped_count <- gsub("[a-z]|[(]|[)]|_", "", phospho_count)
  # this is the length of the first fragment (i.e. start distance):
  phospho_position <- frag_lengths[1]
  #print(phospho_position)
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

  # deals with multiple sites:
  #phospho_sites <- print(phospho_sites, collapse=";")
return(list("sites"=phosho_sites, "seqs"=centred_seqs))
}

# helper function for filtering at site level
filter_site_prob <- function(input_sequence, threshold){
  # conditional within the brackets...
  # put into a function for cleaning up all individual sites!
  input_seq <- strsplit(as.character(input_sequence),"\\)|\\(", perl=TRUE)[[1]]

  input_seq <- sapply(1:length((input_seq)), function(i)
    ifelse(input_seq[i] < threshold, "", input_seq[i]))

  #replace the site scores with a (ph) when meeting criteria
  output_seq <- paste(sapply(1:length((input_seq)), function(i)
    ifelse(input_seq[i] >= threshold & input_seq[i] <= 1,
           "(ph)", input_seq[i])), collapse="")

return(output_seq)
}
