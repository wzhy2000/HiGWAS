library(gwas.lasso)

phe.out <- "bls.test.simple.10k.phe"  
snp.out <- "bls.test.simple.10k.snp"


sigsnp <- c(1, 1050, 1060, 5050, 5060);
show(sigsnp);

ret1 <- ret2 <- ret3 <- ret4 <- ret5 <- ret6 <- ret7 <- ret8 <- ret9 <- ret10 <- c();

r<-bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 1000, simu_p=10000, 
		simu_snp_rho = 0.1, 
		simu_rho     = 0.4, 
		simu_sigma2  = 9, 
		simu_mu      = 24, 
		simu_cov_range=c( 0, 1),
		simu_cov_effect = c( 0, 2 ), 
		simu_add_pos   = c( sigsnp[1], sigsnp[2], sigsnp[3]), 
		simu_add_effect= c( 2.2, -2.5, 2.0 ),  
		simu_dom_pos   = c( sigsnp[3], sigsnp[4], sigsnp[5]), 
		simu_dom_effect= c( 2.8, 2.0, -2.5 ),
		simu_t_range = c(-1, 1), 
		plink.format=T,
		debug=F );


snp.out.bed <- r$file.plink.bed
snp.out.bim <- r$file.plink.bim
snp.out.fam <- r$file.plink.fam

ret1 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=1) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret1);

ret2 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=2) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret2);

ret3 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=3) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret3);

ret4 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=4) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret4);

ret5 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=5) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret5);

ret6 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=6) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret6);

ret7 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=7) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret7);

ret8 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=8) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret8);

ret9 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=9) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret9);

ret10 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c("X_1","X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=7, nPiecewise.ratio=10) );	
save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret9, ret10, file="bls-parallel-demo.rdata");
summary(ret10);
