#' QC/diagnostic plotting
#
#' @description Wrappers for a series of plots to be used as diagnostics in
#' in conjunction with phosphoProcessR analyses. Currently 3 plots are possible
#' using this function: 1) Principle Component Analyis (PCA), 2) Hierarchical
#' clustering, 3) density distribution of phosphosite probability. See
#' descriptions of each in the parameter options below and for format
#' specification. See vignette for more information and examples.
#
#' @param data_in Either a numeric vector for plotting density distribution of
#' site probabilities of a data.frame that contains intensities.
#' @param annotation_file A data.frame describing the experimental design.
#' @param legend Include legend in plots? Legend is based on group data in
#' annotation_file Only used for PCA and hierarchical clustering. Options are:
#' TRUE, FALSE. Default = FALSE
#' @param label Include point labels in plots? Points are based on group data in
#' annotation_file. Options: TRUE, FALSE. Default = FALSE
#' @param title What to name the plot? Default = NULL
#' @param pca Perform unsupervised Principle Component Analysis (PCA) and plot
#' results. By default performs Singular Value Decomposition. Requires data_in
#' be a "data.frame" and utilises "label" meta-data from annotation_file for
#' colouring. Options:  TRUE, FALSE. Default = FALSE
#' @param hclust Performs unsupervised hierarchical clustering of samples.
#' Colours sample below plot according to group and numbered by inputs. Requires
#' data_in be a "data.frame" and utilises "label" meta-data from annotation_file
#' for colouring. Options:  TRUE, FALSE. Default = FALSE
#' @param density Plot density distributions of site probabilities. Input must
#' be a numerical vector. Options: TRUE, FALSE. Default = FALSE
#'
#' @examples
#' ## Load the example data set and attach
#'
#' @return Returns pretty plots
#'
#' @export plot_this
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pcaMethods prep pca R2cum scores
#' @importFrom dendextend colored_bars set %>%
#' @importFrom stats as.dendrogram

plot_this <- function(data_in = NULL,
                      annotation_file = NULL,
                      legend = TRUE,
                      label = TRUE,
                      title = NA,
                      pca = FALSE,
                      hclust = FALSE,
                      density = FALSE
                      ){

####///---- check inputs ----\\\###

  if(hclust == TRUE | pca==TRUE){
    if (is.null(annotation_file))
      stop("annotation_file not provided; you must provide a annotation_file
           table when plotting hclust")
    if("samples" %in% colnames(annotation_file) == FALSE){
      stop("A annotation_file must be supplied with a column labelled \"samples\"")
    }
    if("labels" %in% colnames(annotation_file) == FALSE){
      stop("A annotation_file must be supplied with a column labelled \"labels\"")
    }
    if("group" %in% colnames(annotation_file) == FALSE){
      stop("A annotation_file must be supplied with a column labelled \"group\"")
    }

    #check that the annotation_file$labels exist!
    cols_of_interest <- intersect(colnames(data_in), annotation_file$labels)
    set_diff <- setdiff(annotation_file$labels, colnames(data_in))
    if(length(cols_of_interest) != nrow(annotation_file))
      stop("There is not a matching number of samples in the data_in from
        the labels in the annotation_file. Check file names match or
        correct annotation_file provided. ", set_diff, " not found!")

    # establish colours - "Set3" is up to 12 distinct colours
    # must have a minimum of 3, this set-up as follows:
    colors <- brewer.pal(max(length(unique(annotation_file$group)), 3), "Set3")

    if(pca==TRUE){
      # this will plot 3 PCAs with and without labels
      plot_PCA_wrapper(data_in=data_in,
                       title=title,
                       comp1=1, comp2=2,
                       legend=legend,
                       label=label,
                       annotation_file=annotation_file,
                       colors=colors)
    }

    if(hclust == TRUE){
      plot_hclust_wrapper(data_in=data_in,
                          title=title,
                          colors=colors,
                          annotation_file=annotation_file,
                          legend=legend)
    }
  }

  if(density==TRUE){
    # check format
    if(!is.null(data_in) && !is.numeric(data_in) && !is.vector(data_in))
      stop("data_in is not numeric vector. A numeric vector must be supplied
           to plot the density of site probabilities. This can be derived from
           the tidyEvidence function")
    plot_density_wrapper(data_in=data_in,
                         title=title)
  }
}


plot_density_wrapper <- function(data_in=NULL,
                                 title=NA){
  # determine densities
  dens <- stats::density(data_in)
  # initialise plot
  graphics::plot(dens,
       main=paste("Density - ", title, sep=""),
       xlab="site probability")
}

plot_hclust_wrapper <- function(data_in=NULL,
                                title=NA,
                                colors=NA,
                                annotation_file=NA,
                                legend=TRUE
                                ){

  data_in <- t(data_in)
  # relabel hclust by sample number
  rownames(data_in) <- c(seq_len(nrow(data_in)))
  dd <- stats::dist(scale(data_in), method = "euclidean")
  hc <- stats::hclust(dd, method = "ward.D2")
  # labelling of nodes, by samples
  colors_hclust <- colors[annotation_file$group]
  # adding coloured bars
  the_bars <- colors_hclust
  colors_hclust <- sort(colors_hclust)[hc$order]

  # GENERATE DENDROGRAM
  dend <- data_in %>% scale %>% stats::dist(method = "euclidean") %>%
    stats::hclust(method = "ward.D2") %>% as.dendrogram
  dend %>% set("labels_col", value=c(colors_hclust))
  dend %>% graphics::plot(main = title,
                          sub="euclidean + ward(see legend for sample numbers)")
  colored_bars(colors = the_bars, dend = dend, sort_by_labels_order = TRUE)

  if(legend == TRUE){
    legend("topright",
           c(as.character(unique(annotation_file$group)),
             paste(seq_len(length(annotation_file$labels)),
                   annotation_file$labels, sep="-")),
           col=c(colors[seq_along(length(unique(annotation_file$group)))],
                 rep("black", length(annotation_file$labels))),
           pch = c(rep(19, length(unique(annotation_file$group))),
                   rep(0, length(annotation_file$labels))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }

}

plot_PCA_wrapper <- function(data_in=NULL,
                             title=NA,
                             comp1=1,
                             comp2=2,
                             legend=TRUE,
                             label=TRUE,
                             colors=NA,
                             annotation_file=NA){

  # data transformations
  md <- prep(t(data_in), scale = "none", centre = FALSE)
  pc <- pca(md, method="svd", center=FALSE, nPcs=ncol(data_in))
  var_3 <- R2cum(pc)[3] # accumulated variance
  pc_1 <- round(pc@R2[comp1]*100, 2)
  pc_2 <- round(pc@R2[comp2]*100, 2)

  pc_scores <- as.data.frame(scores(pc))
  pc_scores <- data.frame(pc_scores, "group"=annotation_file$group,
                          "labels"=annotation_file$labels)

  # initialise plot:
  # -2 is remove the group and file name variable.
  graphics::plot(1, type="n", xlim=c(min(pc_scores[comp1])-5,
                           max(pc_scores[comp1])+5),
                    ylim=c(min(pc_scores[comp2])-5,
                           max(pc_scores[comp2])+5),
       axes=TRUE,
       xlab=paste("PC", comp1, " - ", pc_1, "%", sep=""),
       ylab=paste("PC", comp2, " - ", pc_2, "%", sep=""),
       main=paste("PCA - ", title, sep="")
       )

  # add labels:
  graphics::abline(h = 0, v = 0, col = "gray", lty = 2)
  # add legend:
  if(legend==TRUE){
    legend("bottomright", c(as.character(unique(pc_scores$group))),
                            col=colors[annotation_file$group],
                          pch = c(rep(19, length(unique(pc_scores$group)))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
  if(label==TRUE){
    graphics::text(pc_scores[,comp1], pc_scores[,comp2], pc_scores$labels,
                   cex=0.5, pos=3, col="black")
  }
  graphics::points(pc_scores[,comp1], pc_scores[,comp2], cex = 1,
                   col = colors[annotation_file$group], pch=19)
}

