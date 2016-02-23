\name{gls.simple}
\alias{gls.simple}

\title{ Scan Simple format dataset by the Group lasso model}

\description{ 
Scan Simple format dataset with longitudinal measurements and detect the joint additive and dominant effects by the Group lasso model (GLS model).
} 

\usage{
gls.simple(file.phe, 
      file.snp, 
      Y.prefix, 
      Z.prefix, 
      covar.names,
      refit = TRUE, 
      add.used = TRUE, 
      dom.used = TRUE, 
      fgwas.filter = FALSE, 
      options = NULL)}

\arguments{
\item{file.phe}{Input CSV file containing longitudinal measurements and covariate data.}
\item{file.snp}{Input CSV file containing chromosome, position and genotypes of SNPs.}
\item{Y.prefix}{String, the prefix of column names of response variables.}
\item{Z.prefix}{String, the prefix of column names of measure time.}
\item{covar.names}{String or vector of string indicating the column names of covariates. }
\item{refit}{Logical, if TRUE, the refit procedure should be processed. \cr Default: TRUE.}
\item{add.used}{Logical, if TRUE, additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, If TRUE, the SNPs will be filtered by the fGWAS model before the lasso method is applied. \cr Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}. See the details in \code{\link{bls.plink}}.
}
}

\value{
An object of class \code{GLS.ret}, see "Details" in the \code{\link{gls.plink}}.
}

\details{

The simple format is introduced in the section of \code{\link{bls.simple}}. Other details can be found in the section of \code{\link{gls.plink}}

}

\author{Zhong Wang and Nating Wang }

\references{
(1)Group lasso model\cr
Li, J., Wang, Z., Li, R., & Wu, R. (2015). Bayesian group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics, 9(2), 640-664.\cr
(2)fGWAS model\cr
Li, J., Wang, Z., Li, Y., & Wu, R. (2010). Functional genome-wide association studies of longitudinal traits. Handbook of Adaptive Designs in Pharmaceutical and Clinical Development.\cr
}

\examples{

# Create simulation dataset using the default parameters
gls.phe.out <- "gls.test.simple.phe"
gls.snp.out <- "gls.test.simple.snp"

r.sim <- gls.simulate(gls.phe.out, gls.snp.out, simu_n=600, simu_p=10);

# Call the Grpup lasso model to estimate the effects.
r.gls <- gls.simple(gls.phe.out, gls.snp.out, 
                Y.prefix="Y", 
                Z.prefix="Z", 
                covar.names="X", 
                fgwas.filter = FALSE, 
                options=list(nPiecewise.ratio=0) );	

# Show the significant SNPs and effects,
summary(r.gls);

# Plot the results
plot(r.gls);

}

\keyword{functions}
\keyword{GLS model}
\keyword{simple}


