# gwas.lasso

implementing LASSO methods to identify significant SNPs and estimate joint genetic effects for non-longitudinal traits or longitudinal traits in GWAS.

# Reference

[1]. Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. *Bioinformatics*, 27(4), 516-523.

[2]. Li. J., Wang, Z., Li, R., Wu, R.(2015).Bayesian Group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. *The Annals of Applied Statistics*. 9(1).

## Abstract:

The GWAS Lasso package is developed to identify significant SNPs that control phenotypic variation and estimate their additive and dominant genetic effects based on the Bayesian Lasso or Group Lasso model. The package provides two statistical models, **BLS** can detect the association using one single time measured phenotypc and **GLS** solve the association based on the longitudinal phenotype.

## Document

> 1) Vignette (https://github.com/wzhy2000/gwas.lasso/blob/master/gwaslasso-vignette.pdf)

> 2) Manual (https://github.com/wzhy2000/gwas.lasso/blob/master/gwaslasso-manual.pdf)

## Installation Instructions:

### Required software and packages
    
> 1. R (http://www.r-project.org/)
    
> 2. Package [snpStats](http://bioconductor.org/packages/release/bioc/html/snpStats.html), nlme, snowfall.

Please install the required R packages before you install the fGWAS package. After the  installation of the dependencies, please install the **gwas.lasso** as following steps.

### Install fGWAS on LINUX or Mac OSX

```
git clone https://github.com/wzhy2000/gwas.lasso.git

cd gwas.lasso

R CMD INSTALL package

```

### Install fGWAS on Windows

1) Please download windows package from (https://github.com/wzhy2000/fGWAS/raw/master/windows/gwas.lasso.zip)

2) Install the package in R GUI by selecting the menu "Packages|Install package(s) from local zip files..."

## Usage Instructions

fGWAS is an R package which provides:

> 1) Loading the genotype data(SNP) from PLINK data files or simple SNP data table.

> 2) Loading the longitudinal phenotype data(traits) from CSV file with the covariate file or the measure time file.

> 3) Scaning SNP data set to estimate log-likelihood ratio and the genetic effetcs of each genotype.

> 4) Detecting the significant SNPs and export the results.

> 5) Drawing the genetic effects for each significant SNP.


The following codes show how to call above steps in R.

We don't attach any data set in the package, so here we use the simulation to generate the phenotype taits andgenotype data. The simulation function returns a list containing one phenotype object and one genotype object.

```
library(gwas.lasso);
r <- fg.simulate("Logistic", "AR1", 2000, 500, 1:7, sig.pos=250 );
```

Call SNP scaning in a short range (245:255) using 'fgwas' method. 

```
obj.scan <- fg.snpscan(r$obj.gen, r$obj.phe, method="fgwas", snp.sub=c(245:255) );
obj.scan;
```

Plot Manhattan figure for all SNPs in a PDF file.

```
plot(obj2.scan, file.pdf="temp.fwgas.obj2.scan.pdf");
```

Select significant SNPs and plot the varing genetic effects in PDF.

```
tb.sig <- fg.select.sigsnp(obj2.scan, sig.level=0.001, pv.adjust = "bonferroni")
plot.fgwas.curve( obj2.scan, tb.sig$INDEX, file.pdf="temp.fwgas.obj2.curve.pdf");
```

All functions and examples in the fGWAS are available in the manual (https://github.com/wzhy2000/gwas.lasso/blob/master/gwaslasso-manual.pdf).
