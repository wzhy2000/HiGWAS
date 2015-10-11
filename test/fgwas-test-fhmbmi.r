library(gwas.lasso)

#Yale BulldongN
file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1.bed"  
file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1.bim"  
file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1.fam"  
file.phe.nonlong  <- "/home/zw224/f/bmi/bmi-pheno-age-mean.csv";
file.phe.long    <- "/home/zw224/f/bmi/bmi-phenos-sex-longtime.csv"

#louise
#file.plink.bed <- "/home2/zw224/w/bmi/FHS-bmi-v1.bed"
#file.plink.bim <- "/home2/zw224/w/bmi/FHS-bmi-v1.bim"
#file.plink.fam <- "/home2/zw224/w/bmi/FHS-bmi-v1.fam"
#file.phe.long  <- "/home2/zw224/w/bmi/bmi2-phenos-mean.csv"

# stempede
#file.phe.long  <- "/work/03350/tg826494/test/bmidata/bmi-pheno-age-mean.csv";

if(1)
{
	plink.fhm <- load_plink_binary(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.nonlong, verbose=TRUE );

	snp.range<- c(1:500)+10000;

	snp.mat <- as( plink.fhm$snp.mat$genotypes[, snp.range, drop=F ], "numeric");

	ret1 <- fgwas.scan( plink.fhm$phe.mat, t(snp.mat), Y.name="Y", covar.names=c("X"), options=list(nParallel.cpu=1), longitudinal = FALSE);

	save( ret1, file="fgwas-test-fhmbmi1.rdata");

	summary(ret1)

	plot(ret1, fig.prefix="fgwas-test-fhmbmi1");
}

if(1)
{
	plink.fhm <- load_plink_binary(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, verbose=TRUE );

	snp.range<- c(1:500)+10000;

	snp.mat <- as( plink.fhm$snp.mat$genotypes[, snp.range, drop=F ], "numeric");
	
	phe.mat <- plink.fhm$phe.mat;
	Z <- phe.mat[,c(2:27)];
	phe.mat[,c(2:27)] <- 2*(Z-min(Z,na.rm=T))/(max(Z,na.rm=T)-min(Z,na.rm=T)) -1
	
	ret2 <- fgwas.scan( phe.mat, t(snp.mat), Y.name="Y", Z.name="Z", covar.names=c("X"), options=list(nParallel.cpu=1), longitudinal = TRUE);

	save( ret2, file="fgwas-test-fhmbmi2.rdata");

	summary(ret2)

	plot(ret2, fig.prefix="fgwas-test-fhmbmi2");
}
