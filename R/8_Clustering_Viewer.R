#' @title Plots the numbers of cells of each population (cluster or metacluster)
#'
#' @description This function aims to visualize the number of cells associated to each population. The populations can be clusters or metaclusters
#'
#' This representation displays the population in the X-axis and the total number of associated cells in the Y-axis.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the population to use. By default, all the population are used
#' @param level a character value indicating the type of population plotted. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param sort a boolean value indicating if population must be sorted by the number associated cluster
#' @param color.metadata a character value specifying the metadata used to color the barplot. By default, color.metadata is set to NULL (the barplot color is uniform)
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotPopulationCounts <- function(CYTdata,
                                 population = NULL,
                                 level = c("clusters", "metaclusters"),
                                 sort = TRUE,
                                 color.metadata = NULL) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.")}
  checkmate::qassert(population, c("0", "S*"))
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  checkmate::qassert(sort, "B1")
  checkmate::qassert(color.metadata, c("0", "S1"))
  
  if (level == "clusters"){
    if (length(CYTdata@Clustering@clusters)==0) { stop("Error : Clustering step required before vizualization") }
    cellcount = CYTdata@Clustering@cellcount
  }
  if (level == "metaclusters"){
    if (length(CYTdata@Metaclustering@metaclusters)==0) { stop("Error : Metaclustering step required before vizualization") }
    cellcount = CYTdata@Metaclustering@cellcount
  }
  
  if (is.null(population)) { population = rownames(cellcount) }
  cellcount = cellcount[rownames(cellcount) %in% population, ]
  
  pop.effectif = apply(cellcount, 1, sum)
  
  if (!is.null(color.metadata)){
    
    if (!color.metadata %in% colnames(CYTdata@metadata)){ stop("Error : 'color.metadata' argument
    is not a biological condition present in metadata (", paste0(colnames(CYTdata@metadata), collapse=","), ")") }
    
    matrix.count = merge(CYTdata@metadata, t(cellcount), by = "row.names")
    matrix.count = plyr::ddply(matrix.count, color.metadata, function(x){
      return(colSums(x[,rownames(cellcount)]))
    })
    matrix.count = reshape2::melt(matrix.count, id = c(color.metadata))
    colnames(matrix.count) = c("condition", "population", "count")
    
    if (sort) {
      matrix.count$population = factor(matrix.count$population,
                                       levels = rownames(cellcount)[order(pop.effectif, decreasing=TRUE)])
    }
    
    plot <- ggplot2::ggplot(data = matrix.count,
                            ggplot2::aes(x = population,
                                         y = count,
                                         fill = condition)) +
      ggplot2::geom_bar(stat="identity", position = "stack", color="black") +
      viridis::scale_fill_viridis(option = "turbo", discrete = TRUE)
    
  } else {
    
    matrix.count = data.frame("count" = pop.effectif, "population" = rownames(cellcount))
    
    if (sort) {
      matrix.count$population = factor(matrix.count$population,
                                       levels = rownames(cellcount)[order(pop.effectif, decreasing=TRUE)])
    }
    
    plot <- ggplot2::ggplot(data = matrix.count,
                            aes(x = population,
                                y = count)) +
      ggplot2::geom_bar(stat="identity")
  }
  
  if (is.null(color.metadata)) {titre = paste("Number of cells per", level)
      
  } else {titre = paste("Number of cells per", level, "with", color.metadata, "proportion")}
  
  plot <- plot +
    ggplot2::ggtitle(titre) +
    ggplot2::xlab(level) + ggplot2::ylab("Number of cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.key = ggplot2::element_blank())
  
  return(plot)
}







