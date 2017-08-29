# gwas.lasso

implementing LASSO methods to identify significant SNPs and estimate joint genetic effects for non-longitudinal traits or longitudinal traits in GWAS.

# Reference

[1]. Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. *Bioinformatics*, 27(4), 516-523.

[2]. Li. J., Wang, Z., Li, R., Wu, R.(2015).Bayesian Group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. *The Annals of Applied Statistics*. 9(1).

## Abstract:

The GWAS Lasso package is developed to identify significant SNPs that control phenotypic variation and estimate their additive and dominant genetic effects based on the Bayesian Lasso or Group Lasso model. The package provides two statistical models, **BLS** can detect the association using one single time measured phenotypc and **GLS** solve the association based on the longitudinal phenotype.

## Document

> 1) Vignette (https://github.com/wzhy2000/gwas.lasso/blob/master/gwaslasso.vignette.pdf)

> 2) Manual (https://github.com/wzhy2000/gwas.lasso/blob/master/gwaslasso.manual.pdf)

## Installation Instructions:

### Required software and packages
    
> 1. R (http://www.r-project.org/)
    
> 2. Package [snpStats](http://bioconductor.org/packages/release/bioc/html/snpStats.html), nlme, snowfall.

Please install the required R packages before you install the fGWAS package. After the  installation of the dependencies, please install the **gwas.lasso** as following steps.

### Install package on LINUX or Mac OSX

```
git clone https://github.com/wzhy2000/gwas.lasso.git

cd gwas.lasso

R CMD INSTALL package

```

### Install package on Windows

1) Please download windows package from (https://github.com/wzhy2000/gwas.lasso/raw/master/windows/gwas.lasso.zip)

2) Install the package in R GUI by selecting the menu "Packages|Install package(s) from local zip files..."

## Usage Instructions

GWAS lasso is an R package which provides:

> 1) Two Lasso models to analyze the joint genetice effects accumulated by the multiple significant SNPs.

> 2) **GLS** model which is used to associate SNPs with the longitudinal phenotype data(traits).

> 3) **BLS** model which is used to associate SNPs with the single measured phenotype.

> 4) Data analysis pipeline starting from PLINK genotype data, or Simple format genotype data, or SNP matrix. 

> 5) Detecting the significant SNPs and export the results.

> 6) Drawing the genetic effects for each significant SNP.


The following codes show how to call above steps in R.

We don't attach any data set in the package, so here we use the simulation to generate the phenotype taits andgenotype data. The simulation function returns a list containing one phenotype object and one genotype object.

```
library(gwas.lasso);
## generate for BLS model
bls.simulate(“bls.phe.csv”, “bls.gen.csv”);
## generate the longitudinal traits for GLS model
gls.simulate (“gls.phe.csv”, “gls.gen.csv”);;
```

Call SNP scaning using **BLS** model. 

```
r.bls <- bls.simple(“bls.phe.csv”, “bls.gen.csv”, Y.name="Y", covar.names=c("X_1","X_2"));
r.bls;

```

Call SNP scaning using **GLS** model. 

```
r.gls <- gls.simple(“gls.phe.csv”, “gls.gen.csv”,Y.prefix="Y",Z.prefix="Z", covar.names=c("X_1","X_2"));
r.gls;

```

Plot genetic effects for all SNPs in a PDF file.

```
plot (r.bls, fig.prefix="bls-ret");
plot (r.gls, fig.prefix="gls-ret");
```

Show the significant SNPs and effects

```
summary(r.bls);
summary(r.gls);
```

All functions and examples in the GWAS.lasso package are available in the manual (https://github.com/wzhy2000/gwas.lasso/blob/master/gwaslasso.manual.pdf).
