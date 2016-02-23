\name{gls.snpmat}
\alias{gls.snpmat}

\title{ scan matrix dataset by the Group lasso model }

\description{ 
  Scan matrix dataset with longitudinal measurements and detect the joint additive and dominant effects by the Group lasso model (GLS model).
} 

\usage{
gls.snpmat(phe.mat,
      snp.mat,
      Y.prefix,
      Z.prefix,
      covar.names, 
      refit = TRUE,
      add.used = TRUE,
      dom.used = TRUE,
      fgwas.filter = FALSE,
      options = NULL)}

\arguments{
\item{phe.mat}{Matrix, phenotypic data containing longitudinal measurements and covariate data.}
\item{snp.mat}{Matrix, genotypic data containing chromosome, position and genotypes of SNPs.}
\item{Y.prefix}{String, prefix of column names of response variables.}
\item{Z.prefix}{String, prefix of column names of measurement times.}
\item{covar.names}{String or vector of string indicating column names of covariates. }
\item{refit}{ Logical, if TRUE, the refit procedure should be processed. \cr Default: TRUE.}
\item{add.used}{Logical, if TRUE, additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, If TRUE, the SNPs will be filtered by the fGWAS model before the lasso method is applied. \cr Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}.\cr
see details.}
}

\value{
An object of class \code{GLS.ret}, see "Details" in the \code{\link{gls.plink}}.
}

\details{
(1) Matrix dataset\cr
The phenotypic data and genotypic data are stored in the matrix format.\cr

(2) fGWAS model\cr
A SNP filter strategy named fGWAS is proposed in this package. This method is used to detect the
marginal genetic effects. If the parameter \code{fgwas.filter} is true, the fGWAS method will be
employed to remove the nonsignificant SNPs. Once SNPs are filtered by the fGWAS method, significant 
SNPs are remained to do LASSO process.\cr

(3) Varible selection and refit\cr
A two-stage procedure based on Bayesian Lasso or Group Lasso is proposed, including variable selection 
and refit process. It is necessary to refit the lasso model after variable selection, where only selected SNPs 
are included in the final process. The parameter \code{refit} indicates whether to carry out the refit function.\cr

(4) Control parameters\cr
The control parameters include the following items.\cr
\itemize{
\item{\code{nParallel.cpu}}    {Default:0, the number of the CPU used to calculate.}
\item{\code{nPiecewise.ratio}} {Default:2, the ratio of grouping SNPs.}
\item{\code{nMcmcIter}}        {Default:2000, the number of iteration of the MCMC algorithm execution. }
\item{\code{fBurnInRound}}     {Default:0.3, a ratio to cut off unstable iteration in order to make the result precise.}
\item{\code{fRhoTuning}}       {Default:0.095, a ratio used for narrowing the scope of parameter rho to reduce parameter oscillation.}
\item{\code{fQval.add}}        {Default:0.05, }
\item{\code{fQval.dom}}        {Default:0.05, the Q-values of additive and dominant genetic effects respectively.}
\item{\code{fgwas.cutoff}}     {Default:0.05, a p-value cutoff for fGWAS model to select the significant marginal SNPs.}
\item{\code{debug}}            {Default:FALSE, logical value, indicating whether or not debug information is outputted.}
}
}

\author{Zhong Wang and Nating Wang }

\references{
(1)Group lasso model\cr
Li, J., Wang, Z., Li, R., & Wu, R. (2015). Bayesian group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics, 9(2), 640-664.\cr
(2)fGWAS model\cr
Li, J., Wang, Z., Li, Y., & Wu, R. (2010). Functional genome-wide association studies of longitudinal traits. Handbook of Adaptive Designs in Pharmaceutical and Clinical Development.\cr
}

\examples{

gls.phe.out <- "gls.test.simple.phe";
gls.snp.out <- "gls.test.simple.snp";
r.sim <- gls.simulate( gls.phe.out, gls.snp.out, simu_n=400, simu_p=10)

tb.phe <- read.csv(gls.phe.out, header=TRUE);
# Set IDs as the row name
rownames(tb.phe) <- tb.phe[,1];
tb.phe <- tb.phe[,-1];

tb.snp <- read.csv(gls.snp.out);

r.gls <- gls.snpmat(tb.phe, tb.snp, 
                Y.prefix="Y", 
                Z.prefix="Z", 
                covar.names="X", 
                fgwas.filter = FALSE );
	
summary(r.gls);

plot(r.gls);

}

\keyword{functions}
\keyword{GLS model}
