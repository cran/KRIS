#' Read the binary PLINK format (BED, BIM, and FAM)
#'
#' @description Require the complete set of 3 files in the binary
#' PLINK format. It includes BED file, BIM file and BAM file. For more
#' information about the binary PLINK format, please check in the manual of
#' PLINK.
#'
#' @param bed A path of BED file
#' @param bim A path of BIM file
#' @param fam A path of FAM file
#' @param only.snp If TRUE, the function to read only SNP matrix, otherwise all
#' files are loaded. The default value is FALSE.
#'
#' @return The list containing the matrices of \code{$snp}, \code{$snp.info},
#' and \code{$ind.info}.
#' \itemize{
#' \item \code{$snp} is a SNP matrix from BED file.
#' \item \code{$snp.info} is a data.frame of SNP information from BIM file.
#' \item \code{$ind.info} is a data.frame of individual information from FAM file.
#' }
#'
#' @details For more details about the binary PLINK format, please check
#' \url{http://zzz.bwh.harvard.edu/plink/binary.shtml}
#'
#' @export
#' @md
#'
#' @importFrom utils read.table
#'
#' @seealso \code{\link{write.bed}}
#'
#' @examples
#'
#' \donttest{
#' #Use the example files embedded in the package.
#' bed <- system.file("extdata", "example_SNP.bed", package="KRIS")
#' bim <- system.file("extdata", "example_SNP.bim", package="KRIS")
#' fam <- system.file("extdata", "example_SNP.fam", package="KRIS")
#' snp <- read.bed(bed, bim, fam )
#'
#' #Check the objects inside 'snp'
#' ls(snp)
#'
#' #Preview $snp
#' print(snp$snp[1:10, 1:10])
#'
#' #Preview $snp.info
#' head(snp$snp.info)
#'
#' #Preview $ind.info
#' head(snp$ind.info)
#' }
read.bed <- function(bed, bim, fam, only.snp = FALSE){

  ret = NA

  if (!file.exists(bed)){
    cat(paste0("Error: BED file doesn't exist: ",bed,"\n"))
    return(ret)
  }
  if (!file.exists(bim)){
    cat(paste0("Error: BIM file doesn't exist: ",bim,"\n"))
    return(ret)
  }
  if (!file.exists(fam)){
    cat(paste0("Error: FAM file doesn't exist: ",fam,"\n"))
    return(ret)
  }
  #Read BIM file
  snp.info <- read.table(bim,header=FALSE)
  colnames(snp.info) = c("chr","ID","GD","position","allele1","allele2")

  #Read FAM file
  ind.info = read.table(fam,header=FALSE)
  colnames(ind.info) = c("FamID","IndID","PatID","MatID","sex","phenotype")

  no.ind = dim(ind.info)[1]
  no.snp = dim(snp.info)[1]

  if (only.snp == TRUE){
    snp.info = NA
    ind.info = NA
  }

  #Read BED file
  fh = file(bed, 'rb')
  #Read the first three bytes to check file format
  buff = readBin(fh, what="raw",n=3)
  if (sum(buff[1:2] == c('6c','1b')) != 2){
    cat(paste0("Error: BED file is not in a correct format: ",bed,"\n"))
    return(ret)
  }
  no.byte.to.read = NA
  no.loop = NA
  if (buff[3] == '01'){
    no.byte.to.read = ceiling(no.ind/4.0)
    no.loop = no.snp
  }else{
    no.byte.to.read = ceiling(no.snp/4.0)
    no.loop = no.ind
  }

  snp = readBin(fh, what="raw",n=(no.byte.to.read*no.loop))
  close(fh)

  length.snp = length(snp)
  no.part=10
  sub.length.snp = round(length.snp / no.part)

  tmp.snp = c()
  for (i in 1:no.part){
    if (i != no.part){
      idx1 = (i-1)*sub.length.snp + 1
      idx2 = i*sub.length.snp
      sub.snp =snp[idx1:idx2]
    }else{
      idx1 = (i-1)*sub.length.snp + 1
      idx2 = length.snp
      sub.snp =snp[idx1:idx2]
    }

    sub.snp.bits = matrix(as.integer(rawToBits(sub.snp)),ncol=2,byrow=T)

    tmp.sub.snp.bits = sub.snp.bits[,1]+sub.snp.bits[,2]
    tmp.sub.snp.bits[which((sub.snp.bits[,2]-sub.snp.bits[,1])<0)] = NA
    tmp.snp = c(tmp.snp,tmp.sub.snp.bits)
  }

  snp = matrix(tmp.snp,ncol=no.snp,byrow=F)[1:no.ind,]

  if (only.snp == FALSE){
    return(list("snp"=snp,"snp.info"=snp.info,"ind.info"=ind.info))
  }else{
    return(list("snp"=snp))
  }
}


#' Write an list of SNP object to the binary PLINK format (BED, BIM, and FAM)
#'
#' @description Write a SNP object to the files in the binary
#' PLINK format. For more information about the binary PLINK format, please
#' check in the manual of PLINK.
#'
#' @param object An object of SNP is a list consisting of 3 matrices, see the
#' \emph{Details} section for more details.
#' @param file A prefix of output files for BED, BIM and FAM to be saved.
#'
#' @return \code{NULL}.
#'
#' @details The \code{object} should contain:
#' \itemize{
#' \item \code{$snp} is a SNP matrix from BED file.
#' \item \code{$snp.info} is a data.frame of SNP information from BIM file.
#' \item \code{$ind.info} is a data.frame of individual information from FAM file.
#' }
#'
#' For more details about the binary PLINK format, please check
#' \url{http://zzz.bwh.harvard.edu/plink/binary.shtml}
#'
#' @export
#'
#' @importFrom utils write.table
#'
#' @seealso \code{\link{read.bed}}
#'
#' @examples
#'
#'
#' #Load example data
#' data(example_SNP)
#'
#' #Save 'simsnp' to the file as defined in 'save.file'
#' save.file <- file.path(tempdir(),"new_SNP")
#' write.bed(simsnp , save.file)

write.bed <- function(object, file){

  bed = paste0(file,'.bed')
  bim = paste0(file,'.bim')
  fam = paste0(file,'.fam')
  if (file.exists(bed)){
    cat(paste0("Overwrite the existed bed file\n"))
  }
  if (file.exists(bim)){
    cat(paste0("Overwrite the existed bim file\n"))
  }
  if (file.exists(fam)){
    cat(paste0("Overwrite the existed fam file\n"))
  }

  #Write BIM file
  if (is.null(object$snp.info)){
    cat(paste0("Couldn't find object$snp.info: a bim file was not created\n"))
  }else{
    write.table(object$snp.info,file=bim,quote=F,row.names=F,col.names=F)
  }

  #Write FAM file
  if (is.null(object$ind.info)){
    cat(paste0("Couldn't find object$ind.info: a fam file was not created\n"))
  }else{
    write.table(object$ind.info,file=fam,quote=F,row.names=F,col.names=F)
  }

  #Header of BED file
  header = c(0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0)

  no.ind = dim(object$ind.info)[1]
  no.snp = dim(object$snp.info)[1]

  fill.up = 4-(no.ind %% 4)
  tmp = NA
  if (fill.up != 4){
    tmp = matrix(rep(0,no.snp*fill.up),ncol=no.snp)
  }

  vec = object$snp
  if (!anyNA(tmp)){
    vec = rbind(vec,tmp)
  }

  vec = as.vector(vec)
  buff = matrix(rep(0,length(vec)*2),ncol=2)
  buff[which(vec == 1),1] = 0
  buff[which(vec == 1),2] = 1
  buff[which(vec == 2),1] = 1
  buff[which(vec == 2),2] = 1
  buff[which(is.na(vec)),1] = 1
  buff[which(is.na(vec)),2] = 0
  vec = NULL
  buff = c(header,as.vector(t(buff)))
  buff=packBits(as.raw(buff))

  fh = file(bed, 'wb')
  writeBin(buff,fh)
  close(fh)

  invisible(NULL)
}




