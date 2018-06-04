#' Synthetic dataset containing single nucleotide polymorphisms (SNP)
#'
#' The simsnp is the simulated dataset which consists of 3,000 independent SNPs
#' and 753 individuals belonging to one of three populations (250 individuals
#' each) and 3 outlying individuals. The pairwise genetic distance between
#' populations was set to Fst=0.01 as in Balding (1995).
#'
#' \describe{
#'   \item{ind.info}{A character matrix of 753 rows and 6 columns representing
#'   individuals and individual information respectively. The columns of
#'   ind.info represents sample_ID, family_ID, father_ID, mother_ID, gender,
#'   and phenotype respectively.}
#'   \item{snp.info}{A character matrix of 3,000 rows and 6 columns representing
#'   SNPs and SNP information respectively. The columns of snp.info represents
#'   SNP_CHR (chromosome), SNP_ID, centimorgan, position, allele1, and allele2
#'   respectively.}
#'   \item{snp}{A numeric matrix of 753 rows and 3,000 columns representing
#'   individuals and SNPs respectively. The matrix snp contains the number 0, 1,
#'   and 2 representing SNP in additive coding.}
#' }
#'
#' @name simsnp
#' @docType data
#' @format A list with 3 objects
#' @seealso \code{\link{sample_labels}}
#' @usage data(example_SNP)
#' @keywords simsnp
#' @references
#' Balding, D.J., and Nichols, R.A. (1995). A method for quantifying
#' differentiation between populations at multi-allelic loci and its
#' implications for investigating identity and paternity. Genetica 96, 3-12.
"simsnp"

#' Synthetic dataset containing population labels for the dataset simsnp.
#'
#' A dataset contains a character vector of 753 elements containing labels or
#' populations of 753 individuals which they belong. Three populations and
#' outliers were labeled as "pop1", "pop2", "pop3", and "outlier".
#'
#' @name sample_labels
#' @docType data
#' @format A vector with 753 elements.
#' @seealso \code{\link{simsnp}}
#' @usage data(example_SNP)
#' @keywords sample_labels
"sample_labels"

NULL
