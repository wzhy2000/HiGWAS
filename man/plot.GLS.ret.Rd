\name{plot.GLS.ret}
\alias{plot.GLS.ret}

\title{Plot the results of the GLS model}

\description{ 
  Plot a result object obtained by the GLS method.
} 

\usage{ 
\method{plot}{GLS.ret}(x, y=NULL, ..., fig.prefix=NULL) 
}

\arguments{
\item{x}{An object obtained from \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
\item{y}{An null parameter, not used.}
\item{fig.prefix}{String, the prefix of output file name}
\item{...}{Other parameters}
}

\details{
The plot function can output three types of PDF figure, including: \cr

1) The Manhattan figure exported by the fGWAS method (*.fgwas.pdf). \cr

2) The genetic effects of all SNPs estimated by the variable selection procedure (*.varsel.pdf). \cr

3) The genetic effects of significant SNPs estimated by the refit procedure (*.refit.pdf). \cr

The Manhattan figure gives -log10(p-value) for each SNP. The variable selection only selects the SNPs with -log10(p-value) 
greater than the threshold value specified in the control parameters.\cr

In the GLS model, the figures of genetic effects will output the time-varying additive and dominant curves for each significant SNP.
}

\author{Zhong Wang and Nating Wang }

\references{
Li, J., Wang, Z., Li, R., & Wu, R. (2015). Bayesian group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics, 9(2), 640-664.\cr
}

\seealso{See \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}} to generate an object }

\examples{

##r.gls is a result object obtained from gls.simpe, gls.plink, gls.snpmat
##e.g. r.gls <- gls.simple(...);
#
# plot(r.gls, fig.prefix="r.gls.pdf");

}

\keyword{plot}
\keyword{GLS model}

