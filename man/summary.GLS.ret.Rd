\name{summary.GLS.ret}
\alias{summary.GLS.ret}

\title{Summarize the result obtained from the GLS model }

\description{ 
Summarize the result object obtained from the GLS method.
} 


\usage{ 
\method{summary}{GLS.ret}(object, \dots) 
}

\arguments{
\item{object}{an object obtained from \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
\item{\dots}{Other parameters}
}

\details{
The \code{summary} command summarizes some tables of the result object \code{r.gls}, and the command exports all of this values to R console.
}

\author{Zhong Wang and Nating Wang }

\references{
Li, J., Wang, Z., Li, R., & Wu, R. (2015). Bayesian group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics, 9(2), 640-664.\cr
}

\examples{
##r.gls is a result object obtained from gls.simpe, gls.plink, gls.snpmat
##e.g. r.gls <- gls.simple(...);
#
# summary(r.gls);
}

\keyword{summary}
\keyword{GLS model}
