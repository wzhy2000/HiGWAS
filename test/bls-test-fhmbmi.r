library(gwas.lasso)

#Yale BulldongN
file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1.bed"  
file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1.bim"  
file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1.fam"  
file.phe.long  <- "/home/zw224/f/bmi/bmi-pheno-mean.csv"  
#file.phe.long  <- "/home/zw224/f/bmi/bmi-pheno-age-mean.csv";

#louise
#file.plink.bed <- "/home2/zw224/w/bmi/FHS-bmi-v1.bed"
#file.plink.bim <- "/home2/zw224/w/bmi/FHS-bmi-v1.bim"
#file.plink.fam <- "/home2/zw224/w/bmi/FHS-bmi-v1.fam"
#file.phe.long  <- "/home2/zw224/w/bmi/bmi2-phenos-mean.csv"

# stempede
#file.phe.long  <- "/work/03350/tg826494/test/bmidata/bmi-pheno-age-mean.csv";

ret1 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="mean", covar.names=c(), fgwas.filter = T,  refit = T , options=list(nParallel.cpu=15, fgwas.cutoff=0.15));

save( ret1, file="bls-test-fhmbmi.rdata");

summary(ret1)

plot(ret1, fig.prefix="bls-test-fhmbmi");

