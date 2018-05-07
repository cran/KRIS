#' Calculate linear principal component analysis (PCA) from numeric data and
#' Single-nucleotide polymorphism (SNP) dataset
#'
#' @description Available for two types of data; numeric data and
#' Single-nucleotide polymorphism (SNP) dataset in additive coding (0, 1, and 2).
#'
#' @param X A data matrix which rows represent samples and columns represent
#' features.
#' @param PCscore To specify whether scaled PCs will be calculated or not. If
#' FALSE, eigenvectors are returned instead. Default = TRUE.
#' @param no.pc A number of PCs to be calculated. If no.pc is set, PCs are
#' patially calculated. Otherwise all PCs are obtained after calculation.
#' Default = NA.
#' @param data.type To specify a type of data matrix X. It can be set to
#' "linear" and "snp". Default = "linear".
#' @param XXT To specify how pricipal components (PCs) are calculated. If TRUE,
#' PCs are calculated from X.t(X), otherwise X is used directly. XXT is useful
#' option especially an input matrix X contains many columns. Enabling this
#' option, it helps to reduce computation complexity. Regardless the option XXT
#' is enable or not, optained PCs are the same. Default = TRUE.
#'
#' @return The returned value is a list with 2 objects, \code{$PC},
#' \code{$evalue},
#' \itemize{
#' \item \code{$PC} is a PC matrix which rows represent samples and columns
#' represent PCs.
#' \item \code{$evalue} is a vector of eigen values.
#' }
#'
#' @export
#' @importFrom stats median sd
#'
#' @import rARPACK
#'
#' @examples
#'
#' #Load simulated dataset
#' \donttest{
#' data(example_SNP)
#'
#' #Using default parameters
#' PCs <- cal.pc.linear(simsnp$snp)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $evalue
#' print(PCs$evalue[1:3])
#'
#' plot3views(PCs$PC[,1:3], sample_labels)
#'
#' #Calculate PCs without PC scores
#'
#' PCs <- cal.pc.linear(simsnp$snp, PCscore = FALSE)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $evalue
#' print(PCs$evalue[1:3])
#'
#' plot3views(PCs$PC[,1:3], sample_labels)
#'
#' #Calculate the top 3 PCs
#' PCs <- cal.pc.linear(simsnp$snp, no.pc = 3)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $evalue
#' print(PCs$evalue[1:3])
#'
#' plot3views(PCs$PC[,1:3], sample_labels)
#' }

cal.pc.linear <- function(X, PCscore = TRUE, no.pc = NA, data.type = "linear",
                          XXT = TRUE){

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


  #Resolve missing value by median
  #cat(paste0("Missing values are set as medians\n"))
  X.median = apply(X,2,median,na.rm=TRUE)
  missing.char=NA
  X = t(apply(X,1,replace.missing,missing=missing.char,rep=X.median))

  if (data.type == "snp"){
    XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
    XSD = (XPi * (1-XPi))^0.5
  }else{
    #Assume linear
    XSD = apply(X,2,sd)
  }

  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]

  if (is.na(no.pc)){
    no.pc = dim(A)[1]
  }

  if (XXT == TRUE){

    #To handle large matrix
    no.col = dim(A)[2]
    if (no.col <= 1500){
      AA=A %*% t(A)
    }else{
      AA = xxt(A,window.size = 1000)
    }

    UDV = svds(AA, k=no.pc)
    eigen.value = UDV$d
    PCs = UDV$u
  }else{
    UDV = svds(A, k=no.pc)
    eigen.value = UDV$d
    PCs = UDV$u
  }


  if (PCscore == TRUE){
    #Calculate the PCscore and the scaled-up PCs
    B.T=diag(1/sqrt(eigen.value)) %*% t(PCs) %*% A
    PCs=A %*% t(B.T)
  }

  return(list("PC"=PCs,"evalue"=eigen.value))
}


