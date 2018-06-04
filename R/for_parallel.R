#' (Internal) Replace missing values with other values,internally used for
#' parallelization
#'
#' @param X An input vector
#' @param missing A charactors representing a missing value
#' @param rep A vector of new values to replace missing values
#'
#' @return A vector with replaced values
#'
replace.missing <- function(X,missing=NA,rep){
  if (is.na(missing)){
    idx = which(is.na(X))
    X[idx] = rep[idx]
  }else{
    idx = which(X == missing)
    X[idx] = rep[idx]
  }
  return(X)
}
