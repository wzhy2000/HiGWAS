library(gwas.lasso)

#Yale BulldongN
file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1.bed"  
file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1.bim"  
file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1.fam"  
file.phe.long  <- "/home/zw224/f/bmi/bmi-phenos-sex-longtime.csv"  

#TACC Stampede
#file.plink.bed <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.bed"  
#file.plink.bim <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.bim"  
#file.plink.fam <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.fam"  
#file.phe.long  <- "/work/03350/tg826494/test/bmidata/bmi-phenos-sex-longtime.csv"  

ret1 <- gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=c("X"), fgwas.filter = T , options=list(nParallel.cpu=15) );	

save( ret1, file="gls-test-fhmbmi.rdata");

summary(ret1);
plot(ret1, fig.prefix="gls-test-fhmbmi");


