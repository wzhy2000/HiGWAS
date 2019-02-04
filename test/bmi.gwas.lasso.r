#run on cbsudanko
library(fGWAS);

file.plink.bed = "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bed"
file.plink.bim = "/home/zw355/proj/gwas2/bmi-c1c2-qc2.bim"
file.plink.fam = "/home/zw355/proj/gwas2/bmi-c1c2-qc2.fam"

obj.gen <- fg.load.plink( file.plink.bed, file.plink.bim, file.plink.fam, plink.command=NULL, chr=NULL, options=list())
obj.gen;
tb.pca <- fg.get.pca(obj.gen, "~/lib/plink1.9/plink")

tb <- read.csv("/home/zw355/proj/gwas2/phe-cov-time-64.csv");
all(rownames(tb$ID)==tb.pca[,2])

tb[, c("X_2","X_3","X_4","X_5","X_6")] <- tb.pca[,c(3:7)]
file.phe.long <- "gwas.lasso.phe.csv"
write.csv(tb, file=file.phe.long, quote=F, row.names=F);

library(Rgtsvm);
selectGPUdevice(1)

#GPU node
file.plink.bed <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.2000.bed"  
file.plink.bim <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.2000.bim"  
file.plink.fam <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.2000.fam"  
file.phe.long <- "gwas.lasso.phe.csv"
library(HiGWAS)

ret1 <- c();
ret1 <- gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=paste("X", 1:6, sep="_"), fgwas.filter = F , gpu.used=T, options=list(nLegendre=2, nMcmcIter=3000) );	

save(ret1, file="bmi.gwas.lasso.gpu.mcmc3000.2000.rdata");



#GPU node
file.plink.bed <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.3700.bed"  
file.plink.bim <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.3700.bim"  
file.plink.fam <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.3700.fam"  
file.phe.long <- "gwas.lasso.phe.csv"
library(HiGWAS)

ret1 <- c();
ret1 <- gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=paste("X", 1:6, sep="_"), fgwas.filter = F , gpu.used=T, options=list(nLegendre=2) );	

save(ret1, file="bmi.gwas.lasso.gpu.3700.rdata");


#GPU node
file.plink.bed <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.6091.bed"  
file.plink.bim <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.6091.bim"  
file.plink.fam <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.6091.fam"  
file.phe.long <- "gwas.lasso.phe.csv"
library(HiGWAS)

ret1 <- c();
ret1 <- gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=paste("X", 1:6, sep="_"), fgwas.filter = F , gpu.used=T, options=list(nLegendre=2) );	

save(ret1, file="bmi.gwas.lasso.gpu.6091.rdata");

library(Rgtsvm);
selectGPUdevice(1)
#GPU node
file.plink.bed <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.29343.bed"  
file.plink.bim <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.29343.bim"  
file.plink.fam <- "/home/zw355/src/fGWAS/test/bmi/gwas.lasso.bmi.29343.fam"  
file.phe.long <- "gwas.lasso.phe.csv"
library(HiGWAS)

ret1 <- c();
ret1 <- gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=paste("X", 1:6, sep="_"), fgwas.filter = F , gpu.used=T, options=list(nLegendre=2, fQval.add  = 0.05, fQval.dom  = 0.05, nMcmcIter=3000) );	

save(ret1, file="bmi.gwas.lasso.gpu.fQval.change.29343.rdata");
