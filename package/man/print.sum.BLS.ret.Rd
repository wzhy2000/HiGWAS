\name{print.sum.BLS.ret}
\alias{print.sum.BLS.ret}

\title{Print the results obtained from the BLS model}

\description{ 
Print a result object obtained from the BLS method.
} 

\usage{ 
\method{print}{sum.BLS.ret}(x, \dots) 
}

\arguments{
\item{x}{An object obtained from \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}}
\item{\dots}{Other parameters}
}

\details{
The print command will print out the result object in the standard format. 
}

\author{Zhong Wang and Nating Wang }

\references{
 Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.
}

\seealso{See \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}} to generate an object }

\keyword{print}
\keyword{BLS model}
