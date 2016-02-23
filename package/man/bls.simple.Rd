\name{bls.simple}
\alias{bls.simple}

\title{Scan simple format dataset by the Bayesian lasso model}

\description{ 
Scan simple format dataset with single measurement and detect the joint additive and dominant effects by the Bayesian lasso model (BLS model).
} 

\usage{ 
bls.simple(file.phe, 
      file.snp, 
      Y.name, 
      covar.names, 
      refit=TRUE,
      add.used=TRUE, 
      dom.used=TRUE, 
      fgwas.filter=FALSE, 
      options=NULL)}

\arguments{
\item{file.phe}{Input CSV file containing one single measurement and covariate data.}
\item{file.snp}{Input CSV file containing chromosome, position and genotypes of SNPs.}
\item{Y.name}{String indicating the column name of response variable.}
\item{covar.names}{String or vector of string indicating the column names of covariates.}
\item{refit}{ Logical, if TRUE, the refit procedure should be processed. \cr Default: TRUE.}
\item{add.used}{Logical, if TRUE, the additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, the dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, If TRUE, the SNPs will be filtered by fGWAS model before the lasso method is applied. \cr Default: FALSE.}
\item{options}{A list containing the control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}. See "Details" in the \code{\link{bls.plink}}.}
}

\value{
An object of class \code{BLS.ret}, see "Details" in the \code{\link{bls.plink}}.
}

\details{
\emph{Simple format dataset}\cr
The simple format is a user-defined format which is designed to store small amount of SNPs 
for the users who do not use PLINK.\cr
The genotypic data are stored in the CSV format, where each line describes a single SNP and must 
start with 2 columns of chromosome information (chromosome number and SNP position).Three genotypes 
(aa=0,Aa=1,AA=2) and missing data (coded as -1 or NA ) are valid SNP values.\cr

}

\author{Zhong Wang and Nating Wang }

\examples{


# Set the filenames for simple format dataset.
bls.phe.out <- "bls.test.simple.phe"; 
bls.snp.out <- "bls.test.simple.snp";
# Makea simulation data using the defalut parameters.
bls.simulate(bls.phe.out, bls.snp.out,simu_n= 500, simu_p=50,
            simu_add_pos=c( 10, 20, 30), simu_dom_pos=c( 30, 40, 50));
# Call Beyasian lasso model to estimate the joint effects.
r.bls <- bls.simple(bls.phe.out,bls.snp.out,
      Y.name="Y",
      covar.names=c(),
      refit=FALSE,
      options=list(nPiecewise.ratio=0));
      
# List the significant information detcted by the BLS model
summary(r.bls);
# Plot the data
plot(r.bls);

}
 
\references{
(1)Beyasian lasso model\cr
Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.\cr
(2)fGWAS model\cr
Li, J., Wang, Z., Li, Y., & Wu, R. (2010). Functional genome-wide association studies of longitudinal traits. Handbook of Adaptive Designs in Pharmaceutical and Clinical Development.\cr
}

\keyword{functions}
\keyword{BLS model}
\keyword{simple}



