#' @title Check and order population names according to clusters/metaclusters
#'
#' @description This functions checks that names in a character vector are all present in a CYTdata's clustering or metaclustering object
#' If order = TRUE, it also orders the character vector according to the order of the levels of the CYTdata's clusters/metaclusters.
#' If checkDuplicates = TRUE, it also checks for duplicates in the character vector 
#' The function is useful, for example, to check that no mistakes were made when writing a list of cluster names.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the population(s) (clusters or metaclusters) to check and potentially order
#' @param level a character value indicating whether to compare population to the clusters or the metaclusters. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param order logical. If TRUE, population will be ordered according to the order of the levels of the clusters/metaclusters
#' @param checkDuplicates logical. If TRUE, will check for duplicates inside the given population character vector. 
#'
#' @return a character vector
#'
#'@export
#'

checkorderPopulation <- function(CYTdata, population, level, order = TRUE, checkDuplicates = TRUE) {
  
  checkmate::qassert(population, c("0","S*"))
  checkmate::qassert(order, "B1")
  checkmate::qassert(level, "S1")
  level = match.arg(level)
  
  if (!is.null(population) && length(population)==0) {
    stop("Error : 'population' argument is an empty vector (length=0, but non NULL).")
  }
  
  checkmate::qassert(checkDuplicates, "B1")
  if (checkDuplicates) {
    populationdup = population[duplicated(population)]
    if (length(populationdup)>0) {
      stop("Error : 'population' argument contains duplicated values (",
           paste0(populationdup, collapse = ", "),
           "). It must be a vector of unique population levels.")
    }
  }
  
  if (level == "clusters") {
    popId = CYTdata@Clustering@clusters
    if (length(popId)==0) {
      stop("Error : 'level' argument requires clustering step to be performed (Clustering@clusters slot is empty).")
    }
  }
  else {
    popId = CYTdata@Metaclustering@metaclusters
    if (length(popId)==0) {
      stop("Error : 'level' argument requires metaclustering step to be performed (Metaclustering@metaclusters slot is empty).")
    }
  }
  
  if (is.null(population)) { population = levels(popId) }
  else {
    popErr = setdiff(population, levels(popId))
    if (length(popErr)>0) {
      stop("Error : 'population' argument providing identifiers not present in ",
           level, " vector (", paste0(popErr, collapse=", "), ").")
    }
    if (order) {
      # order population vector according to levels
      population = levels(popId)[levels(popId) %in% population]
    }
  }
  
  return(population)
}

#' @title Renames clusters or metaclusters within a CYTdata object.
#'
#' @description This function aims to rename clusters or metaclusters stored within a CYTdata object. 
#' It can rename only a subset of the clusters/metaclusters, or all of them. 
#' If merge = TRUE, clusters/metaclusters which present the same name after renaming will be merged.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param level a character value indicating whether to rename clusters or metaclusters. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param from a character vector providing the cluster/metacluster names to replace. By default, all of them are used, in the level order.
#' @param to a character vector providing the new cluster/metacluster names to use
#' @param to logical. If TRUE, clusters/metaclusters with duplicated names will be merged. If FALSE (default), duplicated names will result in an error. 
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

renamePopulations <-  function(CYTdata, level = c("clusters", "metaclusters"), from = NULL, to, merge = FALSE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  from = checkorderPopulation(CYTdata=CYTdata, population=from, level=level, order=FALSE, checkDuplicates=TRUE)
  
  checkmate::qassert(to, "S*")
  checkmate::qassert(merge, "B1")
  
  to_cat = to
  if (length(from)!=length(to)) {
    if (length(to) == 1) { to = rep(to, length(from)) }
    else {
      stop("Error : Length of argument 'from' (", length(from),
           ") and argument 'to' (", length(to), ") must be equal")
    }
  }
  
  if (level=="clusters") {
    popId = CYTdata@Clustering@clusters
    exPalette = CYTdata@Clustering@palette
  }
  else {
    popId = CYTdata@Metaclustering@metaclusters
    exPalette = CYTdata@Metaclustering@palette
  }
  
  cat("\n\nCurrent ", level, " names are (in the order) :", paste0(levels(popId), collapse = ", "), "\n")
  cat("\n\nThe following ", level, " :")
  cat("\n - ", paste0(from, collapse=", "))
  cat("\n\nwill be renamed, in the same order, by :")
  cat("\n - ", paste0(to_cat, collapse=", "))
  
  exLevels = levels(popId)
  newLevels = levels(popId)
  newLevels[match(from, newLevels)] = to
  
  duplicatednewLevels = unique(newLevels[duplicated(newLevels)])
  if (length(duplicatednewLevels)>0){
    if (!merge) {
      stop("Error : After renaming, several ", level, " have the same name (",  paste0(duplicatednewLevels, collapse=", "),
           " ), either by renaming to an already existing and unchanged ", level, " name,
           or by duplicate in the 'to' argument, or both. CYTdata was unchanged (as 'merge' is set to FALSE).")
    }
    else {
      message("Warning : After renaming, several ", level, " have the same name (",  paste0(duplicatednewLevels, collapse=", "),
              " ), either by renaming to an already existing and unchanged ", level, " name,
           or by duplicate in the 'to' argument, or both. These ", level, " will be merged according to their names (as 'merge' is set to TRUE).")
    }
  }
  
  newpopId = plyr::mapvalues(popId, from = from, to = to)
  
  # message("After renaming,  ", level, " levels contained duplicated name(s) and were merged. As it concerns color palette : \n
  #             - If the name of merged ", level, " was already associated to an existing ", level, " name before renaming,
  #             the color of the existing ", level, " is the new color of newly renamed ", level, ".\n
  #             - If the name of merged ", level, " was a new one (argument 'to' contain this name several times).
  #             The color of the resulting  ", level, " is the one of the ", level, " contained in 'from' argument
  #             and which was ordered first (in 'from' argument) among all the ", level, " renamed to this name.")
  
  newpalette = c()
  for (pop in levels(newpopId)) {
    if (sum(newLevels == pop)>1){ # if duplicated
      if (pop %in% exLevels[newLevels == pop]) { v = exPalette[pop] }
      else { pal = exPalette[from[match(TRUE, to==pop)[1]]] }
    }
    else {
      if (pop %in% to) { pal = exPalette[from[to==pop]] }
      else { pal = exPalette[pop] }
    }
    newpalette = c(newpalette, pal)
  }
  names(newpalette) = levels(newpopId)
  
  if (level=="clusters") {
    CYTdata@Clustering@clusters = newpopId
    cat("\n - Computing new cell", level, "count and abundance matrix...")
    CYTdata@Clustering@cellcount = compute.cellcount(CYTdata@Clustering@clusters, CYTdata@samples)
    CYTdata@Clustering@abundance = compute.abundance(CYTdata@Clustering@cellcount)
    CYTdata@Clustering@palette = newpalette
  }
  else {
    CYTdata@Metaclustering@metaclusters = newpopId
    cat("\n - Computing new cell", level, "count and abundance matrix...")
    CYTdata@Metaclustering@cellcount = compute.cellcount(CYTdata@Metaclustering@metaclusters, CYTdata@samples)
    CYTdata@Metaclustering@abundance = compute.abundance(CYTdata@Metaclustering@cellcount)
    CYTdata@Metaclustering@palette = newpalette
  }
  
  validObject(CYTdata) # check if Clustering, Metaclustering slots are ok
  return(CYTdata)
}

#' @title Replaces a CYTdata's clustering slot with its metaclustering slot, or vice-versa
#' 
#' @description This functions copy-pastes the clustering slot into the metaclustering slot or the metaclustering slot into the clustering slot. It does not empty the copy-pasted slot. 
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param direction a string specifying the direction of the change, "clustersToMetaclusters" or "metaclustersToClusters"
#' @param checkOverwrite a boolean value indicating whether to check if a clustering/metaclustering has already been performed on the CYTdata object. 
#' 
#' @return a CYTdata object
#'
#' @export
#'

changeLevel <- function(CYTdata,
                        direction = c("clustersToMetaclusters", "metaclustersToClusters"),
                        checkOverwrite = TRUE) {
  
  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
  
  direction = match.arg(direction)
  checkmate::qassert(direction, "S1")
  checkmate::qassert(checkOverwrite, "B1")
  
  if (direction == "clustersToMetaclusters") {
    if (length(CYTdata@Clustering@clusters)==0) {
      stop("Error : changeLevel with direction argument equal to 'clustersToMetaclusters' can not be performed on object with empty Clustering slot.
           Please perform Clusetring step")
    }
    if (checkOverwrite && length(CYTdata@Metaclustering@metaclusters)!=0){
      reply <- readline(prompt="Metaclustering slot is not empty, do you still want to continue and overwrite (yes or no): ")
      while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
      if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    }
    CYTdata@Metaclustering = methods::new("Metaclustering",
                                          metaclusters = CYTdata@Clustering@clusters,
                                          cellcount = CYTdata@Clustering@cellcount,
                                          abundance = CYTdata@Clustering@abundance,
                                          palette = CYTdata@Clustering@palette,
                                          optional_parameters = append(CYTdata@Clustering@optional_parameters, list("changeLevel" = direction)))
  }
  else {
    if (length(CYTdata@Metaclustering@metaclusters)==0) {
      stop("Error : changeLevel with direction argument equal to 'metaclustersToClusters' can not be performed on object with empty Metaclustering slot.
           Please perform Metaclusetring step")
    }
    if (checkOverwrite && length(CYTdata@Clustering@clusters)!=0){
      reply <- readline(prompt="Clustering slot is not empty, do you still want to continue and overwrite (yes or no): ")
      while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
      if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    }
    CYTdata@Clustering = methods::new("Clustering",
                                      clusters = CYTdata@Metaclustering@metaclusters,
                                      cellcount = CYTdata@Metaclustering@cellcount,
                                      abundance = CYTdata@Metaclustering@abundance,
                                      palette = CYTdata@Metaclustering@palette,
                                      optional_parameters = append(CYTdata@Metaclustering@optional_parameters, list("changeLevel" = direction)))
  }
  validObject(CYTdata)
  return(CYTdata)
}

















