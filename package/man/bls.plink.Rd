\name{bls.plink}
\alias{bls.plink}

\title{ Scan PLINK dataset by the Bayesian lasso model}

\description{ 
Scan PLINK dataset with single measurement and detect the joint additive and dominant effects by the Bayesian lasso Model(BLS Model).
} 

\usage{
bls.plink(file.phe, 
      file.plink.bed, file.plink.bim, file.plink.fam,
      Y.name, 
      covar.names, 
      refit = TRUE, 
      add.used = TRUE, 
      dom.used = TRUE, 
      fgwas.filter = FALSE, 
      options = NULL,
      force.split=TRUE,
      plink.command=NULL)}

\arguments{
\item{file.phe}{Input CSV file containing one single measurement and covariate data.}
\item{file.plink.bed}{Input PLINK data file, a binary file containing genotype information.}
\item{file.plink.bim}{Input PLINK data file, an extended map file.}
\item{file.plink.fam}{Input PLINK data file, family information.}
\item{Y.name}{String indicating column name of response variable.}
\item{covar.names}{String or vector of string indicating column names of covariates. }
\item{refit}{ Logical, if TRUE, the refit procedure should be processed. \cr Default: TRUE.}
\item{add.used}{Logical, if TRUE, additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, if TRUE, the SNPs will be filtered by fGWAS model before the LASSO method is applied. \cr Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}.See blow details. }
\item{force.split}{Logical, if TRUE, the PLINK command will be applied to extract the data for each single chromosome and then the analysis will be performed on these single chromosomes separately. \cr Default: TRUE.}
\item{plink.command}{String, the full path of PLINK command, if not assigned, data extraction maybe failed due to PLINK command can not be found from defulat path. }
}

\value{
An object of class \code{BLS.ret} including the following elements. 
\item{varsel}{Matrix, additive and dominant effects estimated by the procedure of variable selection for each SNP. 
Eleven columns are available, 
\itemize{
\item{(1) chromosome group}{}
\item{(2) position}{}
\item{(3) significant sign of additive effect, 1: siginificant, 0:nonsignificant.}{}
\item{(4) mean value of additive effect}{}
\item{(5) minimum value of addtive effect}{}
\item{(6) maximum value of addtive effect}{}
\item{(7) significant sign of dominant effect, 1: siginificant, 0:nonsignificant.}{}
\item{(8) mean value of dominant effect}{}
\item{(9) minimum value of dominant effect}{}
\item{(10) maximum value of dominant effect}{}
\item{(11) Heritability rate}{}
}
}
\item{varsel_cov}{Matrix, covariate effects estimated by the variable selection procedure, Four columns as follows:
\itemize{
\item{(1) significant sign of covariate effect, 1: siginificant, 0:nonsignificant.}{}
\item{(2) overall mean value of covariate effect}{}
\item{(3) minimum value of covariate effect}{}
\item{(4) maximum value of covariate effect}{}
}
}
\item{refit}{Matrix, additive and dominant effects estimated by the refit procedure for the SNPs selected by the variable selection, Eleven columns are same as the matrix of varible selection(\code{varsel}).}
\item{refit_cov}{Matrix, covariate effects estimated by the refit procedure, 4 columns are same as the matrix of varible selection(\code{varsel_cov}).}
\item{fgwas.filter}{Matrix, the p-values estimated by fGWAS method for each SNP, 5 columns as follows:  
\itemize{
\item{(1) SNP ID}{}
\item{(2) chromosome group}{}
\item{(3) position}{}
\item{(4) likelihood ratio}{}
\item{(5) p-value}{}
}
}
\item{options}{ List, data file names and all calling parameters, including the control parameters defined in \code{options} }
}

\details{
(1) PLINK binary dataset\cr
The PLINK binary dataset includes three files, \emph{bed},\emph{bim} and \emph{fam} file.\cr
The \emph{bed} file is a compressed binary file containing genotype information.\cr
The \emph{bim} file is an extended map file where each line of this file defines a SNP information 
containing exactly 4 columns: chromosome, SNP identifier, genetic distance and base-pair 
position (bp units). \cr
The \emph{fam} file is a family file containing the individuals and phenotypic information including 
six columns: Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; 
other=unknown) and phenotype.\cr
This package does \strong{not} use the phenotype information stored in the fam file. The Phenotype 
values and covariates are stored into a CSV file which should be filled with numerical values and 
missing values encoded as NA.\cr

(2) fGWAS model\cr
A SNP filter strategy named fGWAS is proposed in this package. This method is used to detect the
marginal genetic effects. If the parameter \code{fgwas.filter} is TRUE, the fGWAS method will be
employed to remove the non-significant SNPs. Once SNPs are filtered by the fGWAS method, significant 
SNPs are remained to do lasso process.\cr

(3) Variable selection and refit\cr
A two-stage procedure based on Bayesian lasso or Group lasso is proposed, including variable selection 
and refit process. It is necessary to refit the lasso model after variable selection, where only selected SNPs 
are included in the final process. The parameter \code{refit} indicates whether to carry out the refit function.\cr

(4) Control parameters\cr
The control parameters include the following items.\cr
\itemize{
\item{\code{nParallel.cpu}}    {Default:1, the number of the CPU used to calculate.}
\item{\code{nPiecewise.ratio}} {Default:2, the ratio of grouping SNPs.}
\item{\code{nMcmcIter}}        {Default:2000, the number of iteration of the MCMC algorithm. }
\item{\code{fBurnInRound}}     {Default:0.3, a ratio to cut off unstable iterations in order to obtain the reasonable results.}
\item{\code{fRhoTuning}}       {Default:0.095, a ratio used for narrowing the scope of parameter rho to reduce parameter oscillation.}
\item{\code{fQval.add}}        {Default:0.05, }
\item{\code{fQval.dom}}        {Default:0.05, the Q-values of additive and dominant genetic effects respectively.}
\item{\code{fgwas.cutoff}}     {Default:0.05, a p-value cutoff for fGWAS model to select the significant marginal SNPs.}
\item{\code{debug}}            {Default:FALSE, logical value, indicating whether or not debug information is outputted.}
}

}

\author{Zhong Wang and Nating Wang }

\references{
(1)Beyasian lasso model\cr
Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.\cr
(2)fGWAS model\cr
Li, J., Wang, Z., Li, Y., & Wu, R. (2010). Functional genome-wide association studies of longitudinal traits. Handbook of Adaptive Designs in Pharmaceutical and Clinical Development.\cr
}

\examples{

bls.phe.out <- "bls.test.plink.phe"  
bls.snp.out <- "bls.test.plink.snp"

sigsnp <- c(1:5)*5;

r.sim <- bls.simulate( bls.phe.out, bls.snp.out, 
                simu_grp=1, simu_n= 600, simu_p=50, 
                simu_snp_rho = 0.1, 
                simu_rho     = 0.4, 
                simu_sigma2  = 9, 
                simu_mu      = 24, 
                simu_cov_range=c( 0, 1),
                simu_cov_effect = c( 0, 2 ), 
                simu_add_pos   = c( sigsnp[1], sigsnp[2], sigsnp[3]), 
                simu_add_effect= c( 2.2, -2.5, 2.0 ),  
                simu_dom_pos   = c( sigsnp[3], sigsnp[4], sigsnp[5]), 
                simu_dom_effect= c( 2.8, 2.0, -2.5 ),
                simu_t_range = c(-1, 1), 
                plink.format = TRUE,
                debug = FALSE );

r.bls <- bls.plink( bls.phe.out, 
                r.sim$file.plink.bed,
                r.sim$file.plink.bim,
                r.sim$file.plink.fam,
                Y.name="Y", 
                covar.names=c("X_1", "X_2"), 
                fgwas.filter=FALSE );

summary(r.bls);

plot(r.bls);

}

\keyword{functions}
\keyword{BLS model}
\keyword{PLINK}

