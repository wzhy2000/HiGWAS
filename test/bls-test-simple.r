library(gwas.lasso)

phe.out <- "bls.test.simple.phe"  
snp.out <- "bls.test.simple.snp"
snp.out.bed <- "bls.test.simple.snp.bed"
snp.out.bim <- "bls.test.simple.snp.bim"
snp.out.fam <- "bls.test.simple.snp.fam"


sigsnp <- sample(1:100)[1:5];
show(sigsnp);

ret0 <-ret1<-ret2<-ret3<-ret4<-ret5<-ret6<-c();

bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 600, simu_p=100, 
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


ret0 <- bls.plink( phe.out, snp.out.bed, snp.out.bim, snp.out.fam, Y.name="Y", covar.names=c(), refit=T, fgwas.filter=F );	
summary(ret0);

bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 600, simu_p=2000, 
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
		plink.format=F,
		debug=F );

ret1 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c(), refit=F, options=list(nPiecewise.ratio=0) );	

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret1)
plot(ret1, fig.prefix="bls-test-simple-ret1");

ret2 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c("X_2"), refit=F );	

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret2)
plot(ret2, fig.prefix="bls-test-simple-ret2");

ret3 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c("X_1","X_2"), refit=F , options=list(nParallel.cpu=7) );	

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret3)
plot(ret3, fig.prefix="bls-test-simple-ret3");

ret4 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c("X_1","X_2"), refit=T, add.used=T, dom.used=F, options=list(nParallel.cpu=7) );

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret4)
plot(ret4, fig.prefix="bls-test-simple-ret4");

ret5 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c("X_1","X_2"), refit=T, add.used=T, dom.used=F, fgwas.filter=T, options=list(nParallel.cpu=7) );

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret5)
plot(ret5, fig.prefix="bls-test-simple-ret5");

ret6 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c(), refit=F, add.used=T, dom.used=F, fgwas.filter=F, options=list(nParallel.cpu=7) );

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret6)
plot(ret6, fig.prefix="bls-test-simple-ret6");

