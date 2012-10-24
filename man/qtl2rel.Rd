\name{qtl2rel}
\alias{qtl2rel}

\title{Convert data from R/qtl to QTLrel format}

\description{
  Convert the data for a QTL mapping experiment from the R/qtl format
  (see \url{http://www.rqtl.org}) to that used by QTLrel.
}

\usage{
qtl2rel(cross)
}

\arguments{
 \item{cross}{An object of class \code{"cross"}, as defined by the R/qtl package}
}

\value{
  A list with four components: \code{"ped"} (pedigree information),
  \code{"gdat"} (genotype data), \code{"pdat"} (phenotype data), and
  \code{"gmap"} (genetic map), in the formats used by QTLrel.
}

\details{
  The input cross must by an intercross (class \code{"f2"}).  X
  chromosome data is omitted.

  Simple pedigree information is created, assuming that it's a standard
  intercross.  
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
library(qtl)
data(listeria)
listeria <- listeria[as.character(1:19),]
reldat <- qtl2rel(listeria)
}

\seealso{ \code{\link{rel2qtl}} }

\keyword{manip}
