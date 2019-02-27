#' Perform Differential Phosphorylation Analysis (DPA)
#'
#' @description Given a remapped evidence file from a MaxQuant output, this
#' function performs statistical analysis of phosphorylation data. It performs
#' normalisation, missing value imputation, batch correction utilising SVA and
#' statistical analysis of differential phosphorylation using limma. Currently
#' it only supports automated pairwise analysis between two groups. Will not run
#' limma if more than two groups detected. See limma_analysis parameter.
#' @param phospho_input Data.frame of intensity values
#' @param annotation_file A data.frame describing the experimental design
#' @param normalisation_method Options are "none" or "quantile". Default =
#' "quantile"
#' @param impute Select from the following methods: "FALSE", "knn"
#' @param correct Performs batch correction using surrogate variable analysis
#' (SVA). Recommended to assess utility through PCA. Default = FALSE
#' @param limma_analysis Perform DPE using limma. Default = TRUE
#' @param adjust_method Method used for multiple comparison adjustment of
#' p-values. Options are: "holm", "hochberg", "hommel", "bonferroni", "BH",
#' "BY", "fdr" or "none". See ?p.adjust.methods for a details description and
#' references. Default = "BH"
#' @param input_percent Percentage of data that must be complete prior to
#' performing missing value imputation. Default = 0.75
#' @param pairing Default = "unpaired"
#' @param verbose Print progress to screen. Default = FALSE
#' @examples
#' ## Use the outputs from the first example to demonstrate
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
#' ## Run statistical analysis on tidy and remapped evidence file
#' DPA_out <- phosphoDE(phospho_input = evidence_tidy$intensity,
#'                      annotation_file = malaria_annotation_example,
#'                      correct = TRUE,
#'                      adjust_method = "none")
#'
#' ## View the corrected results
#' head(DPA_out$corrected_data)
#'
#' @return A List of 5 objects. 1) "results" are the limma results, 2)
#' "kinswingr_input" is a data.frame formatted for KinSwingR analysis, 3)
#' "corrected_data" is normalised + SVA corrected intensities 4)
#' "normalised_data" is normalised, before SVA and 5) "imputed_data" is imputed
#' data. Where an option was not selected, will result in NA.
#'
#' @export phosphoDE
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom impute impute.knn
#' @importFrom sva fsva sva
#' @importFrom stats model.matrix p.adjust
#' @importFrom limma topTable lmFit contrasts.fit eBayes makeContrasts
#' @importFrom utils combn

phosphoDE <- function(phospho_input = NULL,
                      annotation_file = NULL,
                      normalisation_method = "quantile",
                      impute = FALSE,
                      correct = FALSE,
                      limma_analysis = TRUE,
                      adjust_method = "BH",
                      input_percent = 0.75,
                      pairing = "unpaired",
                      verbose = FALSE){

#----------------------------------------------
#format checks:
# check the p-value adjustment method is allowable
if((adjust_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                           "fdr", "none")) == FALSE){
  stop("adjust_method must be one of the following methods: \"holm\",
       \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\" or
       \"none\"\n")
  }

if (is.null(phospho_input))
  stop("phospho_input not provided; you must provide an phospho_input table")
if (is.null(annotation_file))
  stop("annotation_file not provided; you must provide a annotation_file table")
if ((impute %in% c("TRUE", "FALSE")) == FALSE)
  stop("impute must be either TRUE or FALSE")
if ((correct %in% c("TRUE", "FALSE")) == FALSE)
  stop("correct must be either TRUE or FALSE")
if ((pairing %in% c("unpaired")) == FALSE)
  stop("pairing must only be unpaired - only method currently supported")

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

  
# annotation_file needs to be updated.
# for DE analysis only need to keep: "group" and "labels"
  
  
cols_of_interest <- intersect(colnames(phospho_input), annotation_file$labels)
set_diff <- setdiff(annotation_file$labels, colnames(phospho_input))

if(length(cols_of_interest) != nrow(annotation_file))
  stop("There is not a matching number of samples in the phospho_input from
        the labels in the annotation_file. Check file names match or
        correct annotation_file provided. ", set_diff, " not found!")

# check length of samples is only two (only two currently supported)
if((limma_analysis == TRUE && length(unique(annotation_file$group)) > 2))
  stop("phosphoDE currently only supports comparison between two groups. There
       are ", length(unique(annotation_file$group)), " groups detected.")

if (verbose) {
  message("Verbosity has been selected!")
}
#----------------------------------------------

# Filter NA's
# If not imputing - only consider complete data

# as exporting potentially empty files for each step, setup as NA
results <- NA
my_data <- NA
phospho_sva_out <- NA
phospho_impute_out <- NA
phospho_norm_out <- NA

if(impute == FALSE){
  # clean up NA's
  na_remove <- as.vector(apply(is.na(phospho_input), 1, sum))
  # keep only rows with no NA's
  keep <- c(na_remove == 0)
  phospho_input <- phospho_input[keep,]
}

# Normalise (currently only Quantile supported)
# option to plug in different normalisation methods here.
phospho_norm_out <- phospho_norm(phospho_input, normalisation_method)
phospho_input <- phospho_norm_out

# missing value imputation
if(impute == TRUE){
  # update function below
  phospho_input <- impute_data(phospho_input,
                               input_percent)
  phospho_impute_out <- phospho_input
}

# Surrogate variable analysis - for correction if required
if(correct == TRUE){
  phospho_input <- correct_data(phospho_input,
                                annotation_file)
  phospho_sva_out <- phospho_input
}

if(limma_analysis == TRUE){
  # perform limma analysis
  # build matrix as per the consensusDE package
  # extract information from annotation_file
  
  group <- factor(annotation_file$group)
  design_table <- data.frame("file"=as.character(annotation_file$labels),
                             "group"=group)
  
  # currently on unpaired is supported
  if(pairing=="unpaired"){
    design <- model.matrix(~group)
  }
  
  # paired, factor pairing
  if(pairing=="paired"){
  #  pairs <- factor(se$pairs)
  #  design.table <- data.frame(design.table,
  #                             "pairs"=se$pairs)
  #  design <- model.matrix(~pairs + group)
  }
  
  # format model.matrix
  rownames(design) <- as.character(annotation_file$labels)
  colnames(design) <- gsub("group", "", colnames(design))
  colnames(design)[1] <- "Intercept"
  
  fit <- lmFit(phospho_input, design)
  contrast_matrix <- build_contrast_matrix(design_table, design)
  all_fit <- contrasts.fit(fit = fit,
                          contrasts = t(data.frame(t(contrast_matrix[,1]),
                                      row.names = colnames(contrast_matrix)[1])))
  
  fit2 <- eBayes(all_fit)
  results <- topTable(fit2, number=nrow(fit2$coefficients))
  
  # p_value correction
  # reformat for kinswingr
  if(adjust_method == "none"){
    my_data <- data.frame("annotation"=rownames(results),
                          "peptide" = NA,
                          "fc" = results$logFC,
                          "pval" = results$P.Value)
  }
  # if a p_adjust method is selected:
  if(adjust_method != "none"){
    my_data <- data.frame("annotation"=rownames(results),
                          "peptide" = NA,
                          "fc" = results$logFC,
                          "pval" = p.adjust(results$P.Value, method = adjust_method))
  }
}

return(list("results" = results,
            "kinswingr_input" = my_data,
            "corrected_data" = phospho_sva_out,
            "normalised_data" = phospho_norm_out,
            "imputed_data" = phospho_impute_out))
}

# helper functions for normalisation method
phospho_norm <- function(data_in, normalisation_method){
  # no normalisation
  if(normalisation_method == "none"){
    return(data_in)
  }
  # quantile normalisation
  if(normalisation_method == "quantile"){
    data_out <- normalize.quantiles(as.matrix(data_in))
    colnames(data_out) <- colnames(data_in)
    rownames(data_out) <- rownames(data_in)
    data_out <- data.frame(data_out)
    return(data_out)
  }
}

impute_data <- function(data_in, input_percent){
  # uses a default seed
  data_out <- data.frame(impute.knn(as.matrix(data_in))$data)
return(data_out)
}

correct_data <- function(data_in, annotation_file){
  # build matrix for model from annotation.file provided
  mod_matrix <- data.frame(
    "sample" = as.factor(as.character(annotation_file$group)))
  # ensures order is correct
  rownames(mod_matrix) <- annotation_file$labels
  mod = model.matrix(~as.factor(sample), 
                     data = mod_matrix)
  mod0 <- model.matrix(~1,
                       data=mod_matrix)
  # fit the SVA iteratively reweighted model, using a model of variables and null of none.
  svaobj <- invisible(sva(as.matrix(data_in), 
                          mod, 
                          mod0, 
                          method="irw"))
  mod_sv <- cbind(mod,
                  svaobj$sv)
  mod0_sv <- cbind(mod0,
                   svaobj$sv)
  # reweight data
  v_2 <- invisible(fsva(as.matrix(data_in), 
                        mod = mod, 
                        sv = svaobj, 
                        newdat = data_in, 
                        method = "exact"))
  data_out <- v_2$new
return(data_out)
}


# This function is directly from consensusDE - call and reuse
# function to automate contrast matrix creation of all pairs
build_contrast_matrix <- function(table_design, design){
  # base will already be levelled - so take the first level
  base_intercept <- levels(table_design$group)[1]
  # what is remaining in the groups (other pairs)
  remaining <- setdiff(levels(table_design$group), base_intercept)
  # if more than two comparisons
  if(length(levels(table_design$group)) > 2){
    # all possible combinations (of remaining pairs):
    all_com <- combn(remaining, 2)
    remaining_com <- c(sapply(seq_len(ncol(all_com)), function(i)
      paste(all_com[1,i], all_com[2,i], sep="-")))
    # names for the comparisons
    real_com <- c(sapply(seq_len(length(remaining)), function(i)
      paste(remaining[i], base_intercept, sep="-")))
    all_c <- paste(c(remaining, remaining_com), collapse=",")
    contrast_cmd <- paste("makeContrasts(", all_c, ",levels=design)",
                          sep="")
    # evaluate command string
    contrast_matrix <- eval(parse(text=contrast_cmd))
    colnames(contrast_matrix) <- c(real_com, remaining_com)
  }
  # if exactly two comparisons
  if(length(levels(table_design$group)) == 2){
    all_c <- levels(table_design$group)[2]
    contrast_cmd <- paste("makeContrasts(", all_c, ",levels=design)",
                          sep="")
    contrast_matrix <- eval(parse(text=contrast_cmd))
    colnames(contrast_matrix) <- paste(levels(table_design$group)[2],
                                       levels(table_design$group)[1], sep="-")
  }
return(contrast_matrix)
}
