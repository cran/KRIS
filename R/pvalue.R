#' Calculate the average fixation index (Fst) between two groups of individuals
#' from Single-nucleotide polymorphism (SNP)
#'
#' @description Fixation index (Fst) calculation was implemented using Hudson
#' method as in Bhatia (2013) and Hudson (1992).
#'
#' @param X A matrix contains the number 0, 1, and 2 representing SNP in
#' additive coding. Rows represent individuals and columns represent SNP.
#' @param idx.p1 An integer vector contains the row indices of first population
#' in the matrix X.
#' @param idx.p2 An integer vector contains the row indices of second population
#' in the matrix X.
#'
#' @return The function returns an average Fst value between 2 specified groups.
#'
#' @references
#' Bhatia, G., Patterson, N., Sankararaman, S., and Price, A.L. (2013).
#' Estimating and interpreting FST: The impact of rare variants. Genome Res. 23,
#' 1514-1521.
#'
#' Hudson, R.R., Slatkin, M., and Maddison, W.P. (1992). Estimation of levels of
#' gene flow from DNA sequence data. Genetics 132, 583-589.

#'
#' @export
#'
#' @seealso \code{\link{fst.each.snp.hudson}}
#'
#' @examples
#'
#'
#' #Load simulated dataset
#' data(example_SNP)
#'
#' idx1 <- which(sample_labels == 'pop1')
#' idx2 <- which(sample_labels == 'pop2')
#' fst <- fst.hudson(simsnp$snp, idx1, idx2)
#'
#' #Print out the Fst value between 'pop1' and 'pop2'
#' print(fst)
#'

fst.hudson <-function(X, idx.p1, idx.p2){
  prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){

    #Pop 1
    G = alleles[idx.p1]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n1=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n1
    #p.B=(no.BB*2 + no.AB)/n1
    #p1 = min(p.A,p.B,na.rm = T)
    p1 = p.A

    #Pop 2
    G = alleles[idx.p2]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n2=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n2
    #p.B=(no.BB*2 + no.AB)/n2
    #p2 = min(p.A,p.B,na.rm = T)
    p2 = p.A


    N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    D = p1*(1-p2) + p2*(1-p1)

    return(c(N,D))
  }

  set.fst = apply(X,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)
  #fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
  fst = mean(set.fst[1,]/set.fst[2,],na.rm = T)

  return(fst)
}

#' Calculate the fixation index (Fst) for all SNPs between two groups of
#' individuals from Single-nucleotide polymorphism (SNP)
#'
#' @description Fixation index (Fst) calculation was implemented using Hudson
#' method as in Bhatia (2013) and Hudson (1992).
#'
#' @param X A matrix contains the number 0, 1, and 2 representing SNP in
#' additive coding. Rows represent individuals and columns represent SNP.
#' @param idx.p1 An integer vector contains the row indices of first population
#' in the matrix X.
#' @param idx.p2 An integer vector contains the row indices of second population
#' in the matrix X.
#'
#' @return The function returns a matrix of pairwise Fst values for all SNPs
#' between 2 specified groups.
#'
#' @references
#' Bhatia, G., Patterson, N., Sankararaman, S., and Price, A.L. (2013).
#' Estimating and interpreting FST: The impact of rare variants. Genome Res. 23,
#' 1514-1521.
#'
#' Hudson, R.R., Slatkin, M., and Maddison, W.P. (1992). Estimation of levels of
#' gene flow from DNA sequence data. Genetics 132, 583-589.
#'
#' @export
#'
#' @seealso \code{\link{fst.hudson}}
#'
#' @examples
#'
#'
#' #Load simulated dataset
#' data(example_SNP)
#'
#' idx1 <- which(sample_labels == 'pop1')
#' idx2 <- which(sample_labels == 'pop2')
#' fst.pairwise <- fst.each.snp.hudson(simsnp$snp, idx1, idx2)
#'
#' #Print out the Fst values of the first three SNPs between 'pop1' and 'pop2'
#' print(fst.pairwise[1:3])
#'

fst.each.snp.hudson <-function(X, idx.p1, idx.p2){
  prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){

    #Pop 1
    G = alleles[idx.p1]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n1=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n1
    #p.B=(no.BB*2 + no.AB)/n1
    #p1 = min(p.A,p.B,na.rm = T)
    p1 = p.A

    #Pop 2
    G = alleles[idx.p2]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n2=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n2
    #p.B=(no.BB*2 + no.AB)/n2
    #p2 = min(p.A,p.B,na.rm = T)
    p2 = p.A


    N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    D = p1*(1-p2) + p2*(1-p1)

    return(c(N,D))
  }

  set.fst = apply(X,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)
  #fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
  vec.fst = set.fst[1,]/set.fst[2,]

  return(vec.fst)
}

