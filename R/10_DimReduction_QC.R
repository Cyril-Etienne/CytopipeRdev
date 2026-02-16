#' @title QC for dimensionality reduction: correlation between pairwise distances in the high-dimensional space and in the reduced space
#'
#' @description Compares distances between cells in the original, high-dimensional space and in the reduced space. This shows the preservation of the global structure in the dimension reduced space. 
#' Theoretically, a "perfect" dimension reduction should give a linear relationships between pairwise distances in both data spaces. 
#' 
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param subsetSize a numeric being the sample size used for computation. Defaults to 1000. Higher values will lead to memory saturation.
#' @param method a character being the method use to compute correlation. Possible values are : "pearson", "kendall", "spearman". Default value is "pearson"
#' @param markers a character vector of the markers to consider for the knn. Defaults to the markers used to compute the dimensionnality reduction (recommended).
#' 
#' @return a list containing correlation coefficient and boxplot of pairwise distances
#'
#' @export

computeDimRed_CorrelationQC <- function(CYTdata,
                                        subsetSize = 1000,
                                        method = c("pearson", "kendall", "spearman"), 
                                        markers = NULL){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (nrow(CYTdata@DimReduction@coordinates)==0) {
    stop("Error : Dimensionality reduction quality control (QC) require dimensionality step to be performed
  but 'DimReduction' slot from 'CYTdata' argument is empty")
  }
  
  checkmate::qassert(subsetSize, c(0,"N1"))
  method = match.arg(method)
  checkmate::qassert(method, "S1")
  
  if(!is.null(subsetSize)) {
    if (subsetSize<=0) {
      stop("Error : 'subsetSize'  argument must be positive integer ")
    }
    if (subsetSize>nrow(CYTdata@DimReduction@coordinates)) {
      stop("Error : 'subsetSize'  argument must be smaller than number of data point.")
    }
    ## Sample randomly subsetSize points to compute pairwise distances
    subIdx = sample(1:nrow(CYTdata@DimReduction@coordinates), subsetSize)
  }
  else { subIdx = 1:nrow(CYTdata@DimReduction@coordinates) }
  
  if(is.null(markers)) {markers = CYTdata@DimReduction@optional_parameters$markers}
  
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE) 
  
  
  ## Compute pairwise distances on the DR and original space for this subsample
  dist.DR = CYTdata@DimReduction@coordinates[subIdx,] %>% dist() %>% as.vector()
  dist.origin = CYTdata@matrix.expression[subIdx,markers] %>% dist() %>% as.vector()
  data.boxplot = data.frame("DR" = dist.DR, "OR" = cut(dist.origin, 50))
  
  ## Plot distance on the original space vs distance on the DR plot for each pair of point as boxplots
  plot <- ggplot2::ggplot(data.boxplot, ggplot2::aes(x=OR, y=DR)) +
    ggplot2::geom_boxplot(fill="slateblue", alpha=0.2) +
    ggplot2::xlab("Distance cuts in original space") + ggplot2::ylab("Distance in reduced space") +
    ggplot2::ggtitle("Pairwise distances in original space vs reduced space") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust=0.5, size = 20, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 20))
  
  ## Report correlation coefficient
  res = stats::cor(dist.DR, dist.origin, method=method)
  
  cat("\n\n - Pairwise distance", method, "correlation between reduced dimension
      and high dimension (Computed across pairwise distance from", subsetSize, "points) : ",
      res, "\n")
  return(list("boxplot" = plot,
              "cor" = res))
}


#' @title QC for dimensionality reduction : Computes proportion of nearest neighbors preservation in reduced dimension in comparison to high dimension
#'
#' @description Computes the fraction of k-nearest neighbors in the original high-dimensional data that are preserved as k-nearest neighbors in the embedding
#' The average is computed across all n points. KNN quantifies preservation of the local, or microscopic structure.
#' Please note that results will vary greatly depending on parameters knn and subsetSize. Please refer to Wang et al., 2023 (https://www.nature.com/articles/s41467-023-37478-w) for quality control of dimension reduction method in cytometry. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param knn a numeric being the number of nearest neighbors to compute. Defaults to 100. 
#' @param subsetSize a numeric being the sample size used for computation. Defaults to 10000, in order to avoid risks of memory shortage,
#' @param markers a character vector of the markers to consider for the knn. Defaults to the markers used to compute the dimensionnality reduction (recommended).
#'
#' @return a numeric proportion
#'
#' @export

computeDimRed_KNNQC <- function(CYTdata, 
                                knn = 100, 
                                subsetSize = 10000, 
                                markers = NULL){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (nrow(CYTdata@DimReduction@coordinates)==0) {
    stop("Error : Dimensionality reduction quality control (QC) require dimensionality step to be performed
  but 'DimReduction' slot from 'CYTdata' argument is empty")
  }
  
  checkmate::qassert(subsetSize, c(0,"N1"))
  checkmate::qassert(knn, "N1")
  checkmate::qassert(markers, c("S*", 0))
  
  if(!is.null(subsetSize)) {
    if (subsetSize<=0) {
      stop("Error : 'subsetSize'  argument must be positive integer ")
    }
    if (subsetSize>nrow(CYTdata@DimReduction@coordinates)) {
      stop("Error : 'subsetSize'  argument must be smaller than number of data point.")
    }
    ## Sample randomly subsetSize points to compute pairwise distances
    subIdx = sample(1:nrow(CYTdata@DimReduction@coordinates), subsetSize)
  }
  else { subIdx = 1:nrow(CYTdata@DimReduction@coordinates) }
  
  if(is.null(markers)) {markers = CYTdata@DimReduction@optional_parameters$markers}
  
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE) 
  
  data.or = CYTdata@matrix.expression[subIdx, markers]
  data.dr = CYTdata@DimReduction@coordinates[subIdx,]
  
  # finding KNN for each element of data.or and data.dr -> Matrix : nrow(data.or) rows and knn columns
  neighborMatrix.or = RANN::nn2(data.or, data.or, knn+1, searchtype = "standard")[[1]][,-1]
  neighborMatrix.dr = RANN::nn2(data.dr, data.dr, knn+1, searchtype = "standard")[[1]][,-1]
  
  
  proportion.vector = sapply(1:nrow(data.or),
                             function(i){
                               same.nn = intersect(neighborMatrix.or[i,], neighborMatrix.dr[i,])
                               return(length(same.nn))
                             })
  res = mean(proportion.vector)/knn
  cat("\n\n - Average proportion of nearest neighbours preserved in reduced dimension
      in comparison to high dimension nearest (Across a", knn, "points neighbourhood and
      with a computation based on a subset of ", subsetSize, "data points ) : ", res)
  return(res)
}

#' @title QC for dimensionality reduction : Computes proportion of nearest clusters preservation in reduced dimension in comparison to high dimension
#'
#' @description Computes the fraction of k-nearest class (clusters) means in the original data that are preserved as k-nearest class means in the reduced data. 
#' Computes the class mean only and averages across all classes. KNC quantifies preservation of the mesoscopic structure.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param knc a numeric being the number of nearest classes to compute. Defaults to 5. 
#' @param subsetSize a numeric being the sample size used for computation.
#' @param markers a character vector of the markers to consider for the knn. Defaults to the markers used to compute the dimensionnality reduction (recommended).
#' 
#' @return a numeric proportion
#'
#' @export
#'

computeDimRed_KNCQC <- function(CYTdata, 
                                knc = 5, 
                                subsetSize = NULL, 
                                markers = NULL ){
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  if (nrow(CYTdata@DimReduction@coordinates)==0) {
    stop("Error : Dimensionality reduction quality control (QC) require dimensionality step to be performed
  but 'DimReduction' slot from 'CYTdata' argument is empty")
  }
  if (length(CYTdata@Clustering@clusters)==0) {
    stop("Error : 'KNclass' Dimensionality reduction quality control (QC) require clustering step to be performed
  but 'Clustering' slot from 'CYTdata' argument is empty")
  }
  
  checkmate::qassert(subsetSize, c(0,"N1"))
  checkmate::qassert(knc, "N1")
  checkmate::qassert(markers, c("S*", 0))
  
  if(!is.null(subsetSize)) {
    if (subsetSize<=0) {
      stop("Error : 'subsetSize'  argument must be positive integer ")
    }
    if (subsetSize>nrow(CYTdata@DimReduction@coordinates)) {
      stop("Error : 'subsetSize'  argument must be smaller than number of data point.")
    }
    ## Sample randomly subsetSize points to compute pairwise distances
    subIdx = sample(1:nrow(CYTdata@DimReduction@coordinates), subsetSize)
  }
  else { subIdx = 1:nrow(CYTdata@DimReduction@coordinates) }
  
  if(is.null(markers)) {markers = CYTdata@DimReduction@optional_parameters$markers}
  
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE) 
  
  
  data.or = cbind.data.frame("clusters" = CYTdata@Clustering@clusters[subIdx],
                             CYTdata@matrix.expression[subIdx,markers])
  data.dr = cbind.data.frame("clusters" = CYTdata@Clustering@clusters[subIdx],
                             CYTdata@DimReduction@coordinates[subIdx,])
  
  class.means.or = plyr::ddply(data.or,
                               "clusters",
                               function(x){ return(colMeans(x[,-1])) })
  class.means.dr = plyr::ddply(data.dr,
                               "clusters",
                               function(x){ return(colMeans(x[,-1])) })
  
  neighborMatrix.or = RANN::nn2(class.means.or, class.means.or, knc+1, searchtype = "standard")[[1]][,-1]
  neighborMatrix.dr = RANN::nn2(class.means.dr, class.means.dr, knc+1, searchtype = "standard")[[1]][,-1]
  
  proportion.vector = sapply(1:nrow(class.means.or),
                             function(i){
                               same.nc = intersect(neighborMatrix.or[i,], neighborMatrix.dr[i,])
                               return(length(same.nc))
                             })
  
  res = mean(proportion.vector)/knc
  cat("\n\n - Average proportion of nearest clusters preserved in reduced dimension
      in comparison to high dimension nearest (Across a", knc, "clusters neighbourhood and
      with a computation based on a subset of ", subsetSize, "data points ) : ", res)
  return(res)
}




