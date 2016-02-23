\name{bls.simulate}
\alias{bls.simulate}

\title{ Generate simulation data of the Beyasian lasso model}

\description{ 
Generate the simulation dataset of the Beyasian lasso model(BLS model).
} 

\usage{
bls.simulate(file.phe.out,
      file.snp.out,
      simu_grp = 1,
      simu_n = 500,
      simu_p =1000, 
      simu_snp_rho = 0.1, 
      simu_snp_missing = 0.002, 
      simu_rho = 0.4, 
      simu_sigma2 = 3, 
      simu_mu = 26, 
      simu_cov_range = c(0,1),
      simu_cov_effect= c( 0, 2 ), 
      simu_add_pos = c( 100, 200, 300), 
      simu_add_effect= c(2.2,-2.5,2.0), 
      simu_dom_pos = c(300, 500, 700), 
      simu_dom_effect= c(2.8,2.0,-2.5), 
      simu_t_range = c(-1,1), 
      plink.format = FALSE,
      debug = FALSE)
}

\arguments{
\item{file.phe.out}{ String, the name of output phenotypic data file. }
\item{file.snp.out}{ String, the name of output genotypic data file. }
\item{simu_grp}{ Integer, the number of groups that SNPs divided. \cr Default: 1.}
\item{simu_n}{ Integer, sample size. \cr Default: 500. }
\item{simu_p}{ Integer, number of SNPs. \cr Default: 1000. }
\item{simu_snp_rho}{ Float, the correlation coefficient between two adjacent SNPs. \cr Default: 0.1.}
\item{simu_snp_missing}{Float, the ratio of missing SNPs. \cr Default: 0.002.}
\item{simu_rho}{ Float, the correlation coefficient between two adjacent time points. \cr Default: 0.4.}
\item{simu_sigma2}{ Float, the variance of AR(1) convariance matrix. \cr Default: 3.}
\item{simu_mu}{ Float, overall mean of phenotypic data. \cr Default: 26. }
\item{simu_cov_range}{ Vector, range of covariate values. \cr Default: c(0,1).}
\item{simu_cov_effect}{ Vector indicating coefficient of the covariates. \cr Default: c(0,2), means two covariates which the effects are 0 and 2.}
\item{simu_add_pos}{ Vector, indicating positions of the significant SNPs with additive effects. \cr Default: c(100, 200, 300).}
\item{simu_add_effect}{ Vector, indicating additive effects of significant SNPs. \cr Default: c(2.2, -2.5, 2).}
\item{simu_dom_pos}{ Vector, indicating positions of the significant SNPs with dominant effects. \cr Default: c(300, 500,  700).}
\item{simu_dom_effect}{ Vector, dominant effects of significant SNPs. \cr Default: c(2.8, 2, -2.5). }
\item{simu_t_range}{ Vector, range of time points. \cr Default: c(-1, 1) }
\item{plink.format}{ Logical, if TRUE, save the SNP data as PLINK format, otherwise, save as simple format. \cr Default: FALSE. }
\item{debug}{ Logical, if TRUE, run the command in debug model which inputs more messages for debugging. \cr Default: FALSE. }
}

\value{
If PLINK format is specified, 4 filenames, including phenotype file,\emph{bed} file, \emph{bim} file and \emph{fam} file, are returned in a list object. 
Otherwise, two simple format filenames, the phenotypic data file and the genotypic data file, are returned.
}

\details{
The function uses the pre-defined parameters to create a simulation dataset. Two data formats are supported in this function, one is PLINK format, which functions are provided by snpStats package, another is simple-format, which is created by this package.\cr

For PLINK format, 4 files are generated. The \emph{bed} file, \emph{bim} file and \emph{fam} file can be accessed by PLINK command or the snpStats package.\cr

For simple format, 2 files are generated. The genotypic data file is encoded 2, 1, 0, NA to indicate genotypes QQ, Qq, qq and missing data respectively. The phenotypic file and genotypic file can be viewed by any text editor.\cr

It is assumed that the trait is controlled by some SNPs, and the positions of these SNPs are specified by the parameter \code{simu_a_pos} and \code{simu_d_pos}. 
The genetic effects of overall mean and causal SNPs are specified by the parameter \code{simu_mu}, \code{simu_a_effect} and \code{simu_d_effect}.\cr

Except genetic effects, covariates can be simulated in this function. The covariate effects are specified by parameter \code{simu_covar_effect}.\cr

}

\author{Zhong Wang and Nating Wang }

\references{
}

\examples{

# Create the Simple format dataset by defulat parameter
bls.simulate("bls.simple.phe", "bls.simple.snp");

# Create the PLINK dataset by defulat parameter
bls.simulate("bls.simple.phe", "bls.simple.snp", plink.format=TRUE);

}

\keyword{functions}
\keyword{BLS model}
\keyword{simulate}


