\name{gls.plink.tped}
\alias{gls.plink.tped}

\title{ The GLS model }

\description{ 
   Identifying significant SNPs and estimating their additive and dominant genetic effects.
} 

\usage{gls.plink.tped(file.phe,file.plink.tped,file.plink.tfam,
Y.prefix,Z.prefix,covar.names,refit=TRUE,add.used=TRUE,dom.used=TRUE,
options=NULL )}

\arguments{
\item{file.phe}{Input phenotypic data file containing longitudinal phenotypic data and covariate data.}
\item{file.plink.tped}{Input file, a file listing genotypes.}
\item{file.plink.tfam}{Input file, a file containing genotype information.}
\item{Y.prefix}{String, the prefix of column names of response variables.}
\item{Z.prefix}{String, the prefix of measurement names.}
\item{covar.names}{String or vector with strings indicating the names of covariates.}
\item{refit}{ Logical, if TRUE, refit should be processed. Default: TURE.}
\item{add.used}{Logical, if TRUE, the program will estimate additive effects of SNPs. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, the program will estimate dominant effects of SNPs. \cr Default: TRUE.}
\item{options}{List, containing nParallel.cpu, nPiecewise.ratio, nMcmcIter, fBurnInRound, fRhoTuning, fQval.add, fQval.dom, fgwas.cutoff, and debug.\cr
 Defaults to "0", "2", "2000", "0.3", "0.095", "0.05", "0.05", "0.05", "F", respectively.}
}

\value{
An object of class GLS.ret. Details to see \code{\link{bls.plink}}.
}

\details{
A two-stage procedure for multi-SNP modeling and analysis in GWASs is proposed, by first filtering the less significant SNPs if the number of SNPs is greatly large, and then formulating Bayesian lasso to select a subset of significant SNPs. Moreover, to ameliorate the bias of the parameter estimates, it is necessary to refit the fGWAS model after variable selection, where only selected SNPs are included in the final model. The parameter \strong{refit} indicates whether to carry out the refit function.\cr
The phenotype file should be filled with numerical values and missing values recorded as NA. In this file one row represents one subject.\cr
Input genotypic files contains two text files: one (TPED) containing SNP and genotype information where one row is a SNP; one (TFAM) containing individual and family information, where one row is an individual. The first 4 columns of a TPED file are chromosome, snp identifier, Genetic distance (morgans) and Base-pair position (bp units).Then all genotypes are listed for all individuals for each particular SNP on each line. 
The TFAM file is just the first six columns of a standard PED file. In otherwords, we have just taken the standard PED/MAP file format, but swapped all the genotype information between files, after rotating it 90 degrees.\cr
}

\author{Zhong Wang and Nating Wang }


\references{
 Li, J., Wang, Z., Li, R., & Wu, R. (2015). Bayesian group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics, 9(2), 640-664.\cr

}

\keyword{internal}



