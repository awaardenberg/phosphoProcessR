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
#' peptide. Recommended to use 0.7, but first assess with distribution as per
#' output table. See vignette for example. Default = 0
#' @param filter_site_method Two methods for filtering peptides by site
#' probability are available. "site" will remove all phosphosites identified
#' with a probability lower that "min_prob". "peptide" will remove all peptides
#' that are predicted to contain a set combination of sites which do not conatin
#' a minimum probability of a site. Options are "site" or "peptide". Default =
#' "site"
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
#' xx$mapping_table is table of remapped peptides and 3) xx$site_probability is a
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
                         min_prob = 0,
                         filter_site_method = "site",
                         return_intensity = TRUE,
                         return_mapping_table = TRUE,
                         return_site_probability = TRUE,
                         verbose = FALSE
                         ){
  #----------------------------------------------
  # format checks
  if ((filter_site_method %in% c("site", "peptide")) == FALSE)
    stop("filter_site_method must be either: site [or] peptide")
  if ((min_prob >= 0 & min_prob <= 1) == FALSE)
    stop("min_prob must be between zero[0] and one[1]")
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
  
  # 1. Build mapping table from evidence file
  # update to be a database here for speed up
  mapping_table <- data.frame("protein" = evidence_file$Proteins,
                              "seq" = evidence_file$Sequence,
                              "mod_seq" = evidence_file$Modified.sequence,
                              "prob" = evidence_file$Phospho..STY..Probabilities
                              )
  # could be a little quicker here?
  mapping_table <- unique(mapping_table)

  # 2. extract peptides with phospho-sites
  contains_ph <- grep("(ph)", mapping_table$mod_seq)
  # strip all the characters and lower case letters into a new column
  mapping_table$stripped <- gsub("[a-z]|[(]|[)]|_", "", mapping_table$mod_seq)

  # 3. identify contaminants
  # reverse peptides will contain "__", contaminants: "CON__", REV: "REV__"
  cont_mapped <- grep("__", mapping_table$protein)

  # 4. remove unwanted rows from mapping_table
  keep <- setdiff(contains_ph, cont_mapped)
  mapping_table_ph <- mapping_table[keep,]
  
  # determine probability distribution here from all sites
  phospho_prob_out <- extract_probability(mapping_table_ph$prob)[[1]]

  # 5. extract data from the fasta file (accessible as table and list now)
  fasta_data <- extract_fasta_data(fasta_file)
  
  # filter sites
  if(return_mapping_table == TRUE){
    # this will differ with dimethyl option/accross experiment
    # site level will remove sites < probability prior to merging and not
    # consider evidence for this site from other experiments
    if (filter_site_method == "site"){
      # filter out sites above a minimium, min_prob, probability
      # speed up with parallel here - definately needs improvement
      #mapping_table_ph_filtered <- bplapply(seq_along(mapping_table_ph$prob),
      #                                       function(i)
      #                                         filter_site_prob(mapping_table_ph$prob[i], min_prob))
      
      mapping_table_ph$prob_filter <- sapply(seq_along(mapping_table_ph$prob),
                                             function(i)
                                  filter_site_prob(mapping_table_ph$prob[i], min_prob))
      keep <- grep("(ph)", mapping_table_ph$prob_filter)
      # keep only peptides that still contain a (ph), there above min_prob
      mapping_table_ph <- mapping_table_ph[keep,]
    }
  
    # if site is not selected
    if (filter_site_method != "site"){
      mapping_table_ph$prob_filter <- mapping_table_ph$mod_seq
    }

    # 6. Get centered position and centered sequence for each site from the database
    # this needs to be vectorised!
    # CHANGE BACK TO BPLAPPY
    centred_sites <- bplapply(seq_len(nrow(mapping_table_ph)), function(x)
                              map_sites(mapping_table_ph,
                                        fasta_file,
                                        fasta_data,
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
    
    # if best mapping:

  }
  
  if(return_intensity == TRUE){
    evidence_remap <- merge(mapped_out, evidence_file, by="Modified.sequence")
    # peptide level will not alter the phosphosite (but only remove collections
    # with less than minimum phosphosite probability)
    # filter out peptides without a minimum probability HERE

    # 9. determine site probabilities
    phospho_prob <- extract_probability(evidence_remap$Phospho..STY..Probabilities)[[2]]
    evidence_remap$phospho_prob_max <- as.numeric(lapply(strsplit(phospho_prob,";"), max))
    # 10. Trim file and converge
    # put each experiment into its own list.. i.e. subset
    experiment_lists <- lapply(seq_along(uniq_experiments), function(x)
                              evidence_remap[evidence_remap$Experiment == 
                                               uniq_experiments[x],])
    names(experiment_lists) <- uniq_experiments
    
    clean_evidence <- bplapply(seq_along(experiment_lists), function(x)
                               tidy_samples(experiment_lists[[x]],
                                            annotation_file,
                                            names(experiment_lists)[x],
                                            filter_site_method = filter_site_method,
                                            min_prob = min_prob))

    # merge experiments into one output, these already have unique names
    # check length to determine if merge is required.
    if(length(clean_evidence) == 1){
      # if just one experiment, only take the first object in the list
      med_out <- clean_evidence[[1]]
    }
    if(length(clean_evidence) > 1){
      med_out <- data.frame(clean_evidence[[1]])
      # check over character on merging lists...
      med_out$ID <- as.character(rownames(med_out))
      # ideally put into a database here
      merge_out <- bplapply(2:length(clean_evidence), function(x)
                          merge(med_out, 
                                data.frame(clean_evidence[[x]],
                                           "ID"=as.character(
                                             rownames(clean_evidence[[x]]))),
                                all=TRUE))
      # now merged, take first object in list
      merge_out <- merge_out[[1]]
      rownames(merge_out) <- merge_out$ID
      # remove ID column
      merge_out <- merge_out[ , !(colnames(merge_out) %in% "ID")]
      med_out <- merge_out
    }
  }

  # Define how returns are selected
  # return_site_probability
  if(return_intensity == FALSE){
    med_out <- NA
  }
  if(return_mapping_table == FALSE){
    mapped_out <- NA
  }
  if(return_site_probability == FALSE){
    phospho_prob_out <- NA
  }
  
return(list("intensity" = med_out,
       "mapping_table" = mapped_out,
       "site_probability" = phospho_prob_out))
}



# helper function for mapping the phospho-positions &
# extracting the recentred sequences
map_sites <- function(mapping_table_ph,
                      fasta_file,
                      fasta_data,
                      window_size,
                      which_row){
  
  
  #mapping_table_ph
  #fasta_file
  #fasta_data
  #which_row
  
  
  #print(which_row)
  #which_row <- 87
  #which_row <- 387
  #which_row <- 27473
  
  
  fasta_protein_name <- fasta_data[[6]]
  fasta_gene_name <- fasta_data[[7]]
  fasta_evidence <- fasta_data[[8]]
  fasta_database <- fasta_data[[4]]
  fasta_length <- fasta_data[[5]]
  fasta_sequences <- fasta_data[[3]]
    
  # 1. determine if is multimapped
  lead_p_mapped <- unlist(strsplit(as.character(mapping_table_ph$protein[which_row]), ";"))
  # check that they are not the same protein!
  lead_p_mapped <- unique(lead_p_mapped)
  
  # multi - mapped? I.e. greater than one lead protein?
  if(length(lead_p_mapped) != 0L){

    # collapse the protein annotation:
    # protein_annotation <- paste(lead_p_mapped, collapse=":")
    #1. select best proteins
    best_proteins <- best_protein(lead_p_mapped)
    # which fasta database rows match the proteins
    row_matches <- sapply(seq_along(best_proteins), function(x)
      which(fasta_protein_name == best_proteins[x]))      
    
    #2. # extract the genes here:
    gene_symbols <- sapply(seq_along(row_matches), function(x)
      fasta_gene_name[row_matches[x]])
###########################################    
    #DO THE FOLLOWING FOR EACH UNIQUE GENE:
    
    #3. Get the vector for the best gene criteria
    best_gene_rows <- c(sapply(seq_along(unique(gene_symbols)), function(i)
             best_gene(unique(gene_symbols)[i],
                       row_matches[which(unique(gene_symbols)[i] == gene_symbols)],
                       fasta_evidence,
                       fasta_database, 
                       fasta_length)))
    
    #need to unlist here as well.
    best_gene_rows <- unlist(best_gene_rows)
    
    # Annotation can be updated now!
    protein_annotation <- paste(fasta_protein_name[best_gene_rows], collapse=":")
    gene_symbols  <- paste(fasta_gene_name[best_gene_rows], collapse=":") 
    protein_database  <- paste(fasta_database[best_gene_rows], collapse=":") 
    protein_evidence  <- paste(fasta_evidence[best_gene_rows], collapse=":") 
    protein_length  <- paste(fasta_length[best_gene_rows], collapse=":") 

    # obtain the stripped peptide sequence (from the original row)
    peptide <- as.character(mapping_table_ph$stripped[which_row])
    mod_peptide <- as.character(mapping_table_ph$prob_filter[which_row])

    # extract the sites:
    # now update the row_matches to the best matches
    row_matches <- best_gene_rows
    # this will map to all the proteins - end up with length of 1 or > for multi-mapped
    phospho_site_list <- lapply(seq_along(row_matches), function(x)
      phosho_positions(peptide, mod_peptide, seq = fasta_sequences[row_matches[x]],
                       window_size = window_size))
    #this is a list now
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
    #EG  "A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN|tr|1|179"
    return(paste(gene_symbols, protein_annotation, phospho_sites, centered_seqs, protein_database, protein_evidence, protein_length, sep="|"))
  }
  
  # 2. if it is not annotated? It will have a length of zero - retunr(NO_MATCH)
  if(length(lead_p_mapped) == 0L){
    #lead_p_mapped <- "NO_MATCH"
    return("NO_MATCH")
  }

}


best_gene <- function(each_gene, each_gene_rows, fasta_evidence, fasta_database, fasta_length){
  #i=3
  #each_gene <- unique(gene_symbols)[i]
  #each_gene_rows <- row_matches[which(unique(gene_symbols)[i] == gene_symbols)]
  #all_genes <- NA
  #fasta_evidence <- fasta_evidence
  
  #each_gene <- unique(gene_symbols)[1]
  #each_gene_rows <- row_matches[which(each_gene == gene_symbols)]
  #fasta_evidence <- fasta_gene_name
  
  # 3. If the genes are the same - determine best ACCESSION for EACH unique gene
  #     1) select gene with best (smallest number) evdience - which gene has best evidence?
  # or removoal all with evidence greater than min...
  each_gene_evidence_in <- sapply(seq_along(each_gene_rows), function(x)
    fasta_evidence[each_gene_rows[x]])
  # remove genes with evidence above the minimum...
  # update vector of row matches
  each_gene_rows <- each_gene_rows[!each_gene_evidence_in > min(each_gene_evidence_in)]
  
  # 4. If database is the same
  #     1) select "tr" over "sp"
  each_gene_database_in <- sapply(seq_along(each_gene_rows), function(x)
    fasta_database[each_gene_rows[x]])
  tr_or_sp <- "tr" %in% each_gene_database_in
  # update the gene rows if tr are present, otherwise leave alone
  if(tr_or_sp == TRUE){
    each_gene_rows <- each_gene_rows[which(each_gene_database_in == "tr")]
  }
  # 5. If both are "sp"
  #     1) select ID without "-" (i.e. not a variant)
  # THIS IS REDUNDANT - caught earlier now
  # 6. If both are variants, i.e. contain "-"
  #     1) select the variant with smallest number, e.e GENE-1 over GENE-2
  # THIS IS REDUNDANT - caught earlier now
  
  # 7. If multiple proteins left
  #     1) select the longest protein
  each_gene_length <- sapply(seq_along(each_gene_rows), function(x)
    fasta_length[each_gene_rows[x]])
  # remove genes above the minimum...
  # update vector of row matches
  each_gene_rows <- each_gene_rows[!each_gene_length > min(each_gene_length)]
  
  # test each gene symbol
  
  # return row matches...
  return(each_gene_rows)   
}


best_protein <- function(proteins_in){
  # vectorise this
  # determine any matches?
  # anything without a -X now gets a -0
  variants <- proteins_in
  variants[which(!grepl("-", variants))] <- 
    paste(variants[which(!grepl("-", variants))], "-0", sep="")
  variants <- as.numeric(gsub(".*-", "", variants))
  
  # strip all the variants to define unique proteins
  proteins <- gsub("-.*", "", proteins_in)
  uniq_proteins <- unique(proteins)
  
  # now determine the best for the full set of proteins and return:
  best_proteins <- c(sapply(seq_along(uniq_proteins), function(i)
    proteins_in[intersect(which(variants == 
                                  min(variants[which(uniq_proteins[i] == 
                                                       proteins)])), 
                          which(uniq_proteins[i] == proteins))]
  ))
  return(best_proteins)
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


# extract probability of each site, given a vector of sequences
extract_probability <- function(phospho_input){
  #
  phospho_prob <- gsub("([[:alpha:]])", "", phospho_input)
  # phospho_out: output the probabilities of all phosphosites
  phospho_prob_out <- gsub(")(", ", ", phospho_prob, fixed=TRUE)
  phospho_prob_out <- gsub("\\)|\\(", "", phospho_prob_out)
  phospho_prob_out <- paste(phospho_prob_out, collapse = ", ")
  phospho_prob_out <- strsplit(phospho_prob_out, ", ")[[1]]
  phospho_prob_out <- as.numeric(phospho_prob_out)
  phospho_prob_out <- phospho_prob_out[!is.na(phospho_prob_out)]
  # from here, fit a mixture model to automatically determine which to keep?
  # continue clean -up
  phospho_prob <- gsub(")(", ";", phospho_prob, fixed=TRUE)
  phospho_prob <- gsub("\\)|\\(", "", phospho_prob)
  phospho_prob[phospho_prob==1] <- "1;0"#if only one site
  # for instances where there are no probabilities, 
  # make a zero (these are multi-mapped)!
  phospho_prob[phospho_prob==""] <- "0;0"
  # determine number of sites
  #phospho_no <- evidence_remap$Modified.sequence
  #phospho_no <- gsub("(ph)", "1", phospho_no, fixed = TRUE)
  #phospho_no <- gsub("([[:alpha:]])", "", phospho_no)
  #phospho_no <- gsub("\\)|\\(|\\_", "", phospho_no)
  #n_p <- as.numeric(phospho_no)
  # this is still in 1, 11, 111 format - clean up
  return(list(phospho_prob_out, phospho_prob))  
}

# function to select tables based on column names provided
tidy_samples <- function(evidence_file_in,
                         annotation_file,
                         experiment,
                         filter_site_method,
                         min_prob){
  
  #evidence_file_in <- experiment_lists[[1]]
  #experiment <- names(experiment_lists)[1]
  # Select columns of interest from annotation file
  annotation_file_filter <- annotation_file[annotation_file$experiment == experiment,]
  
  samples_selected <- evidence_file_in[,colnames(evidence_file_in) %in% annotation_file_filter$samples]
  # this is now experiment specific
  colnames(samples_selected) <- annotation_file_filter$labels
  
  # new trimmed table to work with
  clean_data <- data.frame("site_id" = evidence_file_in$annotation,
                           samples_selected,
                           "max_prob" = evidence_file_in$phospho_prob_max)
  
  
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
  
  ### REVIEW WHERE THIS GOES....
  # keep peptides with a probability greater than min defined
  if (filter_site_method == "peptide"){
    keep <- sapply(seq_len(length(uniq_site)), function(i)
      max(matches[[i]]$max_prob, na.rm = TRUE) >= min_prob)
    matches <- matches[keep]
    match_names <- names(matches)
  }
  
  
  
  # 13. Tidy up the list columns - reformat to sample values (intensity)
  matches <- lapply(seq_len(length(matches)), function(i)
    matches[[i]][, colnames(matches[[i]]) %in% annotation_file_filter$labels])
  
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
  
  fasta_annotation <- data.frame("header"=fasta_headers,
                                 "seq"=fasta_sequences,
                                 "database"=fasta_database,
                                 "length"=fasta_length,
                                 "protein"=fasta_protein,
                                 "gene"=fasta_gene_name,
                                 "evidence"=fasta_evidence)
  
return(list(fasta_annotation,
            fasta_headers,
            fasta_sequences,
            fasta_database,
            fasta_length,
            fasta_protein,
            fasta_gene_name,
            fasta_evidence))
}
