#' @title Computes differential analysis statistics for cell clusters
#'
#' @description This function aims to identify differences of abundance between populations (i.e., clusters or metaclusters) in two sets of samples. 
#'
#' The statistical test used for the comparisons can be defined by users, and may be paired or unpaired. Corrections can be computed for the p-value. 
#' For each population, the p-value, log2 fold-change and effect size relative to the reference comparisonSpl are computed.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the population to use. By default, all the population are used
#' @param level a character value indicating the type of population plotted. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param comparisonSpl a character value providing the names of the samples to be compared to the reference samples
#' @param referenceSpl a character value providing the names of the reference samples
#' @param comparisonTitle a character value providing the name of comparison. Will be used for storage within the CYTdata object and for plotting with plotVolcano
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcoxon.test' (default) or 't.test'
#' @param p.adjust a character value providing the type of p-value adjustment to use.
#' Possible values are: 'holm' (défault), 'none', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'. See p.adjust function for details. 
#' @param paired a boolean value indicating if a paired test should be performed. Defaults to FALSE.
#' @param verbose a boolean value. If TRUE, the function will print additional information. 
#' @param NFSValues a named numeric with blood count (numération formule sanguine) per sample. Optional. If provided, the abundance data will be multiplied by the blood count before comparison.  
#' @param checkOverwrite a boolean value. If TRUE, the function will check if a differential analysis with the same comparisonTitle has already been performed before overwriting. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#' @import rstatix
#'

runDiffAbundanceAnalysis <- function(CYTdata,
                                     population = NULL,
                                     level = c("clusters", "metaclusters"),
                                     comparisonSpl,
                                     referenceSpl,
                                     comparisonTitle = "cmp",
                                     test.statistics = c("wilcox.test", "t.test"),
                                     p.adjust = c("holm", "none", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                                     paired = FALSE,
                                     verbose = FALSE,
                                     NFSValues = NULL,
                                     checkOverwrite = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  
  checkmate::qassert(comparisonSpl, "S+")
  comparisonSpl = checkorderSamples(CYTdata, comparisonSpl, order = FALSE, checkDuplicates = TRUE)
  checkmate::qassert(referenceSpl, "S+")
  referenceSpl = checkorderSamples(CYTdata, referenceSpl, order = FALSE, checkDuplicates = TRUE)
  
  checkmate::qassert(NFSValues, c(0, "N*"))
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(test.statistics, "S1")
  p.adjust = match.arg(p.adjust)
  checkmate::qassert(p.adjust, "S1")
  checkmate::qassert(paired, "B1")
  checkmate::qassert(verbose, "B1")
  
  if (length(comparisonSpl) < 4) warning("Less that 4 samples in the comparison samples.")
  if (length(referenceSpl) < 4) warning("Less that 4 samples in the reference samples.")
  
  if (verbose) { cat("Differential abundance analysis : Compuation of the ", test.statistics,
                     " for: ", paste(comparisonSpl, collapse = ", "),
                     " vs. ", paste(referenceSpl, collapse = ", "), "\n\n") }
  
  checkmate::qassert(checkOverwrite, "B1")
  if (checkOverwrite && comparisonTitle %in% unique(CYTdata@DiffAbundanceAnalysis$Title)) {
    reply <- readline(prompt="Comparison already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    CYTdata@DiffAbundanceAnalysis = subset(CYTdata@DiffAbundanceAnalysis, Title != comparisonTitle)
  }
  
  if (level == "clusters") { abundance = CYTdata@Clustering@abundance }
  else { abundance = CYTdata@Metaclustering@abundance }
  
  abundance = abundance[population, c(comparisonSpl, referenceSpl)]
  popZero = rownames(abundance)[rowSums(abundance)==0]
  if (length(popZero)>0) {
    message("Warning : Following", level, " (", paste0(popZero, collapse=", "),
            ") are not present in sample's subset defined by 'comparisonSpl' and 'referenceSpl' arguments.")
  }
  abundance = abundance[rowSums(abundance)!=0,]
  
  if (!is.null(NFSValues)) {
    if (is.null(names(NFSValues))) {
      stop("Error : NFSValues argument is not a named vector, please see documentation.")
    }
    NFSValues = NFSValues[colnames(abundance)]
    if (sum(is.na(NFSValues))>0) {
      stop("Error : NFSValues does not contain some samples (",
           paste0(names(NFSValues)[is.na(NFSValues)], collapse=", "),
           ")as named elements.")
    }
    abundance = as.matrix(abundance) %*% diag(as.vector(NFSValues))
    abundance = as.data.frame(abundance)
    colnames(abundance) = names(NFSValues)
  }
  
  abundance = reshape2::melt(data.matrix(abundance))
  colnames(abundance) = c("popId", "samples", "value")
  abundance$group = ifelse(abundance$samples %in% referenceSpl, "referenceSpl", "comparisonSpl")
  
  # return(abundance)
  
  switch(test.statistics,
         wilcox.test = { test_fct <- rstatix::wilcox_test },
         t.test = { test_fct <- rstatix::t_test })
  
  statsTest = abundance %>%
    dplyr::group_by(popId) %>%
    test_fct(value ~ group, paired = paired) %>%
    rstatix::adjust_pvalue(method = p.adjust) %>%
    rstatix::add_significance("p")
  
  foldchange = abundance %>%
    dplyr::group_by(popId) %>%
    dplyr::summarise(log2FoldChange = log(mean(value[group == "comparisonSpl"])/
                                            mean(value[group == "referenceSpl"]), 2))
  
  stats = merge(statsTest, foldchange, by="popId")
  
  stats$Title = rep(comparisonTitle, nrow(stats))
  stats$Level = rep(level, nrow(stats))
  if (is.null(stats$p.adj)) { stats$p.adj = rep(NA, nrow(stats)) }
  
  if (nrow(CYTdata@DiffAbundanceAnalysis) == 0) {
    if(verbose) cat("Creating DiffAbundanceAnalysis data.frame with ", comparisonTitle, " statistical results..")
    CYTdata@DiffAbundanceAnalysis = stats
  } else {
    if(verbose) cat("Updating DiffAbundanceAnalysis data.frame with ", comparisonTitle, " statistical results..")
    CYTdata@DiffAbundanceAnalysis = rbind.data.frame(CYTdata@DiffAbundanceAnalysis, stats)
  }
  
  validObject(CYTdata)
  return(CYTdata)
}



#' @title Plots a volcano plot of a differential abundant analysis
#'
#' @description This function aims to visualize the results of a differential abundant analysis using a Volcano plot
#' In such representation, each in dot corresponds to a phenotypic family (population; clusters or metaclusters). Dots are positioned in a two dimensional space where the X-axis represents the log2(fold-change) and the Y-axis represents the -log10 of the p-value.
#' One horizontal line shows the p-value threshold and two vertical lines show the fold-change threshold.
#' The populations are plotted and colored according to their color in the CYTdata. 
#'
#' @param CYTdata a CYTdata object
#' @param comparison a character value containing the name of the comparison to display
#' @param thPvalue a numeric value containing the p-value threshold to use. Defaults to 0.05
#' @param thFoldchange a numeric value containing the fold-change threshold to use for considering changes in population abundance between reference and comparison as significant. Defaults to 2.
#' @param signifOnly a boolean. If TRUE, populations that do not show a significant change of abundance (in terms of pValue and foldchange) will be displayed in grey, with no name. 
#' @param displayAdjust a boolean. If TRUE, adjusted p-values will be used. 
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotVolcano <- function(CYTdata,
                        comparisonTitle,
                        thPvalue = 0.05,
                        thFoldchange = 2,
                        displayAdjust = TRUE,
                        signifOnly = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  stats = CYTdata@DiffAbundanceAnalysis
  
  checkmate::qassert(comparisonTitle, "S1")
  stats = subset(stats, Title == comparisonTitle)
  if (nrow(stats)==0){
    stop("Error : No differential abundant analysis has been performed with ", comparisonTitle, " title. Nothing can be vizualized")
  }
  
  level = unique(stats$Level)
  if (!all(level %in% c("clusters", "metaclusters"))) {
    stop("Error : ", comparisonTitle, " data must  differential abundance analysis on either clusters or metaclustersa
    and 'Level' column must contain either 'clusters' ot 'metaclusters'. Here, ", paste0(level, collapse=", "), " are givenas Level.")
  }
  if (length(level)!=1) {
    stop("Error : ", comparisonTitle, " data contains differential abundance analysis on both clusters and metaclusters.
         Each differential abundance analysis must be performed on either clusters or metaclusters.")
  }
  
  checkmate::qassert(thPvalue, "N1")
  if (!(thPvalue>=0 && thPvalue<=1)) {
    stop("Error : 'thPvalue' argument is a p-value threshold and must be a positive numeric between 0 and 1.")
  }
  checkmate::qassert(thFoldchange, "N1")
  if (thFoldchange<=1) { stop("Error : 'thFoldchange' argument must be a numerical value greater than 1.") }
  
  if (displayAdjust) {
    if (any(is.na(stats$p.adj))) {
      stop("Error : 'displayAdjust' argument set to TRUE but no adjusted p-value computed for ", comparisonTitle, ".")
    }
    else { stats$log10Pvalue = -log10(stats$p.adj) }
  }
  else { stats$log10Pvalue = -log10(stats$p) }
  
  stats$dir = "ns"
  stats$dir[stats$log10Pvalue > -log10(thPvalue) & stats$log2FoldChange > log2(thFoldchange)] = "up" #green
  stats$dir[stats$log10Pvalue > -log10(thPvalue) & stats$log2FoldChange < -log2(thFoldchange)] = "down" #red
  
  print(stats)
  
  if (level == "clusters") {
    palette = CYTdata@Clustering@palette
    stats$popId = factor(stats$popId, levels = levels(CYTdata@Clustering@clusters))
  }
  else {
    palette = CYTdata@Metaclustering@palette
    stats$popId = factor(stats$popId, levels = levels(CYTdata@Metaclustering@metaclusters))
  }
  
  if (signifOnly) {
    signifStats = subset(stats, dir != "ns")
    nosignifStats = subset(stats, dir == "ns")
    plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = signifStats,
                          ggplot2::aes_string(x="log2FoldChange",
                                              y="log10Pvalue",
                                              fill="popId"),
                          shape = 21, stroke = 0, size=8) +
      ggplot2::scale_fill_manual(values = palette) +
      ggrepel::geom_text_repel(data = signifStats,
                               ggplot2::aes_string(x = "log2FoldChange",
                                                   y = "log10Pvalue",
                                                   label = "popId"),
                               color = "black", size = 6, force = 4) +
      ggplot2::geom_point(data = nosignifStats,
                          ggplot2::aes_string(x="log2FoldChange",
                                              y="log10Pvalue"),
                          shape = 21, stroke = 0, fill="grey", size=5)
    
  }
  else {
    stats = dplyr::mutate_at(stats, "popId", as.character)
    plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = stats,
                          ggplot2::aes_string(x="log2FoldChange",
                                              y="log10Pvalue",
                                              fill="popId"),
                          shape = 21, stroke = 0, size=8) +
      ggplot2::scale_fill_manual(values = palette) +
      ggrepel::geom_text_repel(data = stats,
                               ggplot2::aes_string(x = "log2FoldChange",
                                                   y = "log10Pvalue",
                                                   label = "popId"),
                               color = "black", size = 6, force = 4)
  }
  
  max.fc = max(abs(stats$log2FoldChange[is.finite(stats$log2FoldChange)]))
  plot <- plot +
    ggplot2::geom_hline(yintercept = -log10(thPvalue), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(log2(thFoldchange), -log2(thFoldchange)), linetype = "dashed") +
    ggplot2::labs(title = paste(comparisonTitle, "(Volcano plot representation)", sep=" ")) +
    ggplot2::xlab("log2(fold-change)") + ggplot2::ylab("-log10(p-value)") +
    ggplot2::scale_x_continuous(limits = c(-max.fc, max.fc), breaks = seq(-10, 10, 1)) +
    ggplot2::scale_y_continuous(breaks = c(seq(0, 10, 1), -log10(thPvalue))) +
    ggplot2::guides(color=guide_legend(title=level)) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   legend.position = "none")
  return(plot)
}


