\name{summary.BLS.ret}
\alias{summary.BLS.ret}

\title{Summarize the result obtained from the BLS model }

\description{ 
 Summarize the result object obtained from the BLS method.
} 

\usage{ 
\method{summary}{BLS.ret}(object, \dots) 
}

\arguments{
\item{object}{an object obtained by \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}}
\item{\dots}{Other parameters}
}

\details{
The \code{summary} command summarizes some tables of the result object \code{r.bls}, and the command exports all of this values to R console.
}

\author{Zhong Wang and Nating Wang }

\references{
 Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.
}

\seealso{See \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}} to generate an object }

\examples{

##r.bls is a result object obtained from bls.simpe, bls.plink, bls.snpmat
##e.g. r.bls <- bls.simple(...);
#
# summary(r.bls);

}

\keyword{summary}
\keyword{BLS model}
