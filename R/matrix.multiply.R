#' Calculate matrix multipication between a matrix and its
#' transpose for large data.
#'
#' @description Calculate matrix multiplication using "divide and
#' conquer technique", which accelerates the computation to be faster.
#'
#' @param X An input matrix to be processed.
#' @param window.size The window size of matrices to be devided. The default
#' value is 5.
#'
#' @return The multiplication matrix of \code{X.t(X)}.
#'
#' @export
#'
#' @examples
#'
#'
#' #Use the example files embedded in the package.
#' X <-matrix(runif(100), ncol=20)
#' R1  <- xxt(X)
#'
#' #Show the result (R1)
#' print(R1)
#' R2 <- X %*% t(X)
#'
#' #Show the result (R2)
#' print(R2)
#'
#'
#'

xxt <- function(X, window.size = 5){
  n.col = dim(X)[2]
  window.size = as.integer(window.size)
  if (window.size < 1){
    cat(paste0("In matrices.split.by.col(), window.size must be possitive integer\n"))
    return(NULL)
  }
  if (is.null(n.col)){
    cat(paste0("In matrices.split.by.col(), n.col is null.\n"))
    return(NULL)
  }
  if (n.col < window.size){
    cat(paste0("In matrices.split.by.col(), n.col is than window.size\n"))
    return(NULL)
  }

  no.split = n.col %/% window.size

  split.matrices = lapply(1:no.split, function(col){
    X[,((col-1)*window.size+1):(col*window.size)]})

  if ((no.split * window.size) < n.col){
    split.matrices[[length(split.matrices)+1]] = as.matrix(X[,(no.split*window.size+1):n.col])
  }

  #Can be parallel

  partial.mul = NULL
  for (i in 1:length(split.matrices)){
    partial.mul[[i]] = split.matrices[[i]] %*% t(split.matrices[[i]])
  }

  mul = NULL
  for (i in 2:length(partial.mul)){
    if (is.null(mul)){
      mul = partial.mul[[1]] + partial.mul[[i]]
    }else{
      mul = mul + partial.mul[[i]]
    }
  }

  return(mul)
}





