\name{bls.snpmat}
\alias{bls.snpmat}

\title{ Scan matrix dataset by the Bayesian lasso model}

\description{ 
Scan matrix dataset with single measurement and detect the joint additive and dominant effects by the Bayesian lasso model (BLS model).
} 

\usage{
bls.snpmat(phe.mat, 
      snp.mat, 
      Y.name, 
      covar.names, 
      refit = TRUE, 
      add.used = TRUE, 
      dom.used = TRUE, 
      fgwas.filter = FALSE, 
      options = NULL)}

\arguments{
\item{phe.mat}{Matrix, phenotypic data containing one single measurement and covariate data.}
\item{snp.mat}{Matrix, genotypic data containing chromosome, position and genotypes of SNPs.}
\item{Y.name}{String, the column name of response variable.}
\item{covar.names}{String or vector of string indicating the column names of covariates.}
\item{refit}{ Logical, if TRUE, the refit procedure should be processed. \cr Default: TRUE.}
\item{add.used}{ Logical, if TRUE, the additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{ Logical, if TRUE, the dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{ Logical, If TRUE, the SNPs will be filtered by fGWAS model before the lasso method is applied. \cr Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}. See details in \code{\link{bls.plink}}.
}
}

\value{
An object of class \code{BLS.ret}, see "Details" in the \code{\link{bls.plink}}.
}

\details{

The phenotypic data and genotypic data are stored in the matrix. The first two columns in the genotypic matrix must be chromosome and position information and the other columns are SNP data encode by 0,1,2 and NA.

The phenotypic matrix inludes the following columns:\cr
1) Individual ID\cr
2) Response varible(Y)\cr
3) One or more covariates(X, X_1,...\cr

}

\author{Zhong Wang and Nating Wang }


\references{
(1)Beyasian lasso model\cr
Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.\cr
(2)fGWAS model\cr
Li, J., Wang, Z., Li, Y., & Wu, R. (2010). Functional genome-wide association studies of longitudinal traits. Handbook of Adaptive Designs in Pharmaceutical and Clinical Development.\cr
}

\examples{

# Generate the simulation data set.
bls.phe.out <- "bls.test.simple.phe";
bls.snp.out <- "bls.test.simple.snp";
r.sim <- bls.simulate( bls.phe.out, bls.snp.out,simu_n= 500, simu_p=100, 
            simu_add_pos=c( 10, 20, 30), simu_dom_pos=c( 30, 40, 50));

# Load the phenotypic data into a matrix
tb.phe<-read.csv(bls.phe.out, header=TRUE);
# Set IDs as the row name
rownames(tb.phe) <- tb.phe[,1];
tb.phe <- tb.phe[,-1];

# Load the genotypic data into a matrix
tb.snp<-read.csv(bls.snp.out);

# Call Beyasian lasso model to estimate the joint effects.
r.bls <- bls.snpmat(tb.phe, tb.snp, 
      Y.name="Y", 
      covar.names=c("X_1","X_2"), 
      fgwas.filter = FALSE );

# Show the signiciant SNPs and effects.
summary(r.bls);

# Plot the data
plot(r.bls);

}

\keyword{functions}
\keyword{BLS model}

