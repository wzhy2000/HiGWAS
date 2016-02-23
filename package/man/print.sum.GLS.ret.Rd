\name{print.sum.GLS.ret}
\alias{print.sum.GLS.ret}


\title{Print the result obtained from the GLS model}

\description{ 
Print a result object obtained from the GLS method.
} 

\usage{ 
\method{print}{sum.GLS.ret}(x, \dots) 
}

\arguments{
\item{x}{An object obtained from \code{\link{gls.simple}} or \code{\link{gls.plink}} or \code{\link{gls.snpmat}}}
\item{...}{Other parameters}
}

\details{

The print command will print out the result object in the standard format. 
}

\author{Zhong Wang and Nating Wang }

\references{
Li, J., Wang, Z., Li, R., & Wu, R. (2015). Bayesian group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics, 9(2), 640-664.\cr
}

\seealso{See \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}} to generate an object }


\keyword{print}
\keyword{GLS model}