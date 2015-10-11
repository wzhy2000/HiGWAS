library(gwas.lasso)

#file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1.bed"  
#file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1.bim"  
#file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1.fam"  

#Yale BulldongN
#file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1-chr2.bed"  
#file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1-chr2.bim"  
#file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1-chr2.fam"  

#TACC Stampede
file.plink.bed <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.bed"  
file.plink.bim <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.bim"  
file.plink.fam <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.fam"  

#file.phe.long  <- "/home/zw224/f/bmi/bmi-pheno-mean.csv"  
#file.phe.long  <- "/home/zw224/f/bmi/bmi-pheno-age-mean.csv";
file.phe.long  <- "/work/03350/tg826494/test/bmidata/bmi-pheno-age-mean.csv";


ret1<-ret2<-ret3<-ret4<-ret5<-ret6<-c();

ret1 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c(),    fgwas.filter = T,  refit = F , options=list(nParallel.cpu=7));

save( ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-plink.rdata");
summary(ret1)
plot(ret1, fig.prefix="bls-test-plink-ret1");

ret2 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c("X"), fgwas.filter = T, add.used=T, dom.used=F, options=list(nParallel.cpu=7) );

save( ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-plink.rdata");
summary(ret2)
plot( ret2, fig.prefix="bls-test-plink-ret2");

ret3 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c("X"), fgwas.filter = F , options=list(nParallel.cpu=7));	

save( ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-plink.rdata");
summary(ret3)
plot( ret3, fig.prefix="bls-test-plink-ret3");

ret4 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c("X"), fgwas.filter = T , add.used=T, dom.used=F);	

save( ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-plink.rdata");
summary(ret4)
plot( ret4, fig.prefix="bls-test-plink-ret4");

ret5 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c("X"), fgwas.filter = T, refit = F , options=list(nParallel.cpu=7));	

save( ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-plink.rdata");
summary(ret5)
plot( ret5, fig.prefix="bls-test-plink-ret5");

ret6 <- bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c("X"), fgwas.filter = T, refit=F, add.used=F, dom.used=T , options=list(nParallel.cpu=7));

save( ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-plink.rdata");
summary(ret6)
plot( ret6, fig.prefix="bls-test-plink-ret6");

