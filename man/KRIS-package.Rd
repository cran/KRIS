% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kris.R
\docType{package}
\name{KRIS-package}
\alias{KRIS}
\alias{KRIS-package}
\title{Keen and Reliable Interface Subroutines for bioinformatic analysis}
\description{
Provides useful functions which are needed for bioinformatic analysis such as
calculating linear principal components from numeric data and
Single-nucleotide polymorphism (SNP) dataset, calculating fixation index
(Fst) using Hudson method, creating scatter plots in 3 views, handling with
PLINK binary file format, detecting rough structures and outliers using
unsupervised clustering, and calculating matrix multiplication in the faster
way for big data.
}
\details{
The R package \pkg{KRIS} requires the package \pkg{rARPACK}.

Here is the list of functions in the R package \pkg{KRIS}:
\itemize{
\item \code{\link{cal.pc.linear}}
\item \code{\link{cal.pc.projection}}
\item \code{\link{fst.each.snp.hudson}}
\item \code{\link{fst.hudson}}
\item \code{\link{plot3views}}
\item \code{\link{read.bed}}
\item \code{\link{rubikclust}}
\item \code{\link{write.bed}}
\item \code{\link{xxt}}
}

Moreover, here is the list of example datasets in the R package \pkg{KRIS}:
\itemize{
\item \code{\link{simsnp}}
\item \code{\link{sample_labels}}
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://gitlab.com/kris.ccp/kris}
  \item Report bugs at \url{https://gitlab.com/kris.ccp/kris/-/issues}
}

}
\author{
\strong{Maintainer}: Kridsadakorn Chaichoompu \email{kridsadakorn@biostatgen.org}

Authors:
\itemize{
  \item Kristel Van Steen
  \item Fentaw Abegaz
  \item Sissades Tongsima
  \item Philip James Shaw
  \item Anavaj Sakuntabhai
  \item Luisa Pereira
}

}
\keyword{internal}
