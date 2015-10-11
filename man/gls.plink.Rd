\name{gls.plink}
\alias{gls.plink}

\title{ Scan PLINK dataset by the Group lasso model }

\description{ 
Scan PLINK dataset with longitudinal measurements and detect the joint additive and dominant effects by the Group lasso model (GLS model).
} 

\usage{
gls.plink(file.phe,
      file.plink.bed,file.plink.bim,file.plink.fam,
      Y.prefix,
      Z.prefix,
      covar.names,
      refit = TRUE,
      add.used = TRUE,
      dom.used = TRUE,
      fgwas.filter = FALSE, 
      options = NULL,
      force.split=TRUE,
      plink.command=NULL)
}
\arguments{
\item{file.phe}{Input CSV file containing longitudinal measurements and covariate data.}
\item{file.plink.bed}{Input PLINK data file, a binary file containing genotype information.}
\item{file.plink.bim}{Input PLINK data file, an extended map file.}
\item{file.plink.fam}{Input PLINK data file, family information.}
\item{Y.prefix}{String, the prefix of column names of response variables.}
\item{Z.prefix}{String, the prefix of column names of measurement times.}
\item{covar.names}{String or vector of string indicating column names of covariates.}
\item{refit}{ Logical, if TRUE, the refit procedure should be processed. \cr Default: TRUE.}
\item{add.used}{Logical, if TRUE, the additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, the dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, If TRUE, the SNPs will be filtered by fGWAS model before the lasso method is applied. \cr Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}. See the details in \code{\link{bls.plink}}.}
\item{force.split}{Logical, if TRUE, the PLINK command will be applied to extract the data for each single chromosome and then the analysis will be performed on these single chromosomes separately. \cr Default: TRUE.}
\item{plink.command}{String, the full path of PLINK command, if not assigned, data extraction maybe failed due to PLINK command can not be found from defulat path. }

}
\value{
An object of class \code{GLS.ret} including the following elements. 
\item{varsel_add}{Matrix, \strong{additive} effects estimated by the procedure of variable selection for the SNPs selected by fGWAS filter. Twenty-one columns are available, 
\itemize{
\item{(1) chromosome group}{}
\item{(2) position}{}
\item{(3) the sign of \strong{1st} Legendre polynomial of \strong{additive} effcts}{}
\item{(4) the sign of \strong{2nd} Legendre polynomial of \strong{additive} effcts}{}
\item{(5) the sign of \strong{3rd} Legendre polynomial of \strong{additive} effcts}{}
\item{(6) the sign of \strong{4th} Legendre polynomial of \strong{additive} effcts}{}
\item{(7) the norm of mean values of \strong{additive} effects, i,e, L=sqrt(add_1^2 + add_2^2 + add_3^2 + add_4^2)}{}
\item{(8)  the mean value of \strong{1st} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(9)  the mean value of \strong{2nd} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(10) the mean value of \strong{3rd} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(11) the mean value of \strong{4th} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(12) the norm of minimum values of \strong{additive} effects}{}
\item{(13) the minimum value of \strong{1st} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(14) the minimum value of \strong{2nd} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(15) the minimum value of \strong{3rd} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(16) the minimum value of \strong{4th} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(17) the norm of maximum values of \strong{additive} effects}{}
\item{(18) the maximum value of \strong{1st} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(19) the maximum value of \strong{2nd} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(20) the maximum value of \strong{3rd} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(21) the maximum value of \strong{4th} Legendre coefficient of \strong{additive} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
}
}
\item{varsel_dom}{Matrix, \strong{dominant} effects estimated by the procedure of variable selection for the SNPs selected by fGWAS filter. 
21 columns are same as \code{varsel_add}.
}
\item{varsel_cov}{Matrix, covariate effects estimated by the variable selection procedure, 19 columns as follows:
\itemize{
\item{(1) the sign of \strong{1st} Legendre polynomial of \strong{covariate} effcts}{}
\item{(2) the sign of \strong{2nd} Legendre polynomial of \strong{covariate} effcts}{}
\item{(3) the sign of \strong{3rd} Legendre polynomial of \strong{covariate} effcts}{}
\item{(4) the sign of \strong{4th} Legendre polynomial of \strong{covariate} effcts}{}
\item{(5) the norm of mean values of \strong{covariate} effects, i,e, L=sqrt(add_1^2 + add_2^2 + add_3^2 + add_4^2)}{}
\item{(6)  the mean value of \strong{1st} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(7)  the mean value of \strong{2nd} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(8) the mean value of \strong{3rd} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(9) the mean value of \strong{4th} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(10) the norm of minimum values of \strong{covariate} effects}{}
\item{(11) the minimum value of \strong{1st} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(12) the minimum value of \strong{2nd} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(13) the minimum value of \strong{3rd} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(14) the minimum value of \strong{4th} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(15) the norm of maximum values of \strong{covariate} effects}{}
\item{(16) the maximum value of \strong{1st} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(17) the maximum value of \strong{2nd} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(18) the maximum value of \strong{3rd} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
\item{(19) the maximum value of \strong{4th} Legendre coefficient of \strong{covariate} effects in the range of (fQval.add, 1 - fQval.add) of sorted steable iteration.}{}
}
}

\item{refit_add}{Matrix, \strong{additive} effects estimated by the refit procedure for the SNPs selected by variable selection, 21 columns are same as the matrix of varible selection. See also (\code{varsel_add}).}
\item{refit_dom}{Matrix, \strong{dominant} effects estimated by the refit procedure for the SNPs selected by variable selection, 21 columns are same as the matrix of varible selection. See also (\code{varsel_dom}).}
\item{refit_cov}{Matrix, \strong{covariate} effects estimated by the refit procedure for the SNPs selected by variable selection, 19 columns are same as the matrix of varible selection. See also (\code{varsel_cov}).}
\item{fgwas.filter}{Matrix, the p-values estimated by fGWAS method for each SNP, 5 columns as follows:  
\itemize{
\item{(1) SNP ID}{}
\item{(2) chromosome group}{}
\item{(3) position}{}
\item{(4) likelihood ratio}{}
\item{(5) p-value}{}
}
}
\item{varsel_PSRF}{Matrix, the potential scale reduction factors in the variable selection procedure. }
\item{refit_PSRF}{Matrix, the potential scale reduction factors in the refit procedure.}
\item{options}{ List, data file names and all parameters, including the control parameters defined in \code{options} }
}
\details{
(1) PLINK dataset\cr
The PLINK dataset includes three files, \emph{bed},\emph{bim} and \emph{fam} file.\cr
The \emph{bed} file is a compressed binary file containing genotype information.\cr
The \emph{bim} file is an extended map file where each line of this file defines a SNP information 
containing exactly 4 columns: chromosome, SNP identifier, genetic distance and base-pair 
position (bp units). \cr
The \emph{fam} file is a family file containing the individuals and phenotypic information including 
six columns: Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; 
other=unknown) and phenotype.\cr
This package does \strong{not} use the phenotype information stored in the fam file. The phenotype 
values and covariates are stored into a CSV file which should be filled with numerical values and 
missing values encoded as NA.\cr

(2) fGWAS model\cr
A SNP filter strategy named fGWAS is proposed in this package. This method is used to detect the
marginal genetic effects. If the parameter \code{fgwas.filter} is true, the fGWAS method will be
employed to remove the nonsignificant SNPs. Once SNPs are filtered by the fGWAS method, significant 
SNPs are remained to do LASSO process.\cr

(3) Variable selection and refit\cr
A two-stage procedure based on Bayesian Lasso or Group Lasso is proposed, including variable selection 
and refit process. It is necessary to refit the lasso model after variable selection, where only selected SNPs 
are included in the refit process. The parameter \code{refit} indicates whether to carry out the refit function.\cr

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

# Set the filename of simulation data set
gls.phe.out <- "gls.test.simple.phe"  
gls.snp.out <- "gls.test.simple.snp"

# Set the additive effects for 3 SNPs
a_effect <- array(0, dim=c(3,4));
a_effect[1,]<-c( 1.04, 0.885, -2.055, 0.545);
a_effect[2,]<-c( 1.17, -0.20, 0.74, -4.715);
a_effect[3,]<-c( 1.40, -2.25, 1.00,  0.00);

# Set the dominant effects for 3 SNPs
d_effect <- array(0, dim=c(3,4));
d_effect[1,]<-c( 1.49, -2.135, 4.82, 1.425);
d_effect[2,]<-c( 1.045, 1.320, 1.905,  1.535);
d_effect[3,]<-c( 1.265, -1.225, 2.710, -1.96);

# Set the coefficient for 2 covariates
cov_effect <- array(0, dim=c(2,4));
cov_effect[1,]<-c( 2.49, -1.135, 0.82, 0.425);
cov_effect[2,]<-c( -1.045, 2.320, 0.905,  0.535);

#Set the significant positions
sigsnp <- c(1, 3, 5,7,9);

#Create simulation dataset.
r.sim <- gls.simulate( gls.phe.out, gls.snp.out, 
                simu_n= 400, simu_grp=1, simu_p=10, 
                simu_snp_rho=0.4, simu_rho=0.1, simu_sigma2= 4, 
                simu_mu= c(13.395, -3.08, 1.875, -3.195),  
                simu_cov_effect = cov_effect, 
                simu_cov_range  = c(-1,1),
                simu_add_pos    = c( sigsnp[1], sigsnp[2], sigsnp[3] ),
                simu_add_effect = a_effect,  
                simu_dom_pos    = c( sigsnp[3], sigsnp[4], sigsnp[5] ),
                simu_dom_effect = d_effect, 
                simu_z_range    = c(30,60), simu_z_count = c(5,12), 
                plink.format=TRUE, 
                debug=FALSE);

# Call theGroup model to estimate the genetic effects 
r.gls <- gls.plink(gls.phe.out, 
                r.sim$file.plink.bed,
                r.sim$file.plink.bim,
                r.sim$file.plink.fam,
                Y.prefix="Y", 
                Z.prefix="Z", 
                covar.names=c("X_1","X_2"), 
                fgwas.filter = FALSE );	
		

# Show the significant SNPs and effects		
summary(r.gls);

# Plot the results
plot(r.gls);

}

\keyword{functions}
\keyword{GLS model}
\keyword{PLINK}
