library(gwas.lasso)

phe.out <- "bls.test.snpmat.phe"  
snp.out <- "bls.test.snpmat.snp"

sigsnp <- c(11, 22, 33, 444, 55);

bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 400, simu_p=3000, 
		simu_snp_rho = 0.1, 
		simu_rho     = 0.4, 
		simu_sigma2  = 9, 
		simu_mu      = 24, 
		simu_cov_effect = c( 0, 2 ), 
		simu_add_pos   = c( sigsnp[1], sigsnp[2], sigsnp[3]), 
		simu_add_effect= c( 2.2, -2.5, 2.0 ),  
		simu_dom_pos   = c( sigsnp[3], sigsnp[4], sigsnp[5]), 
		simu_dom_effect= c( 2.8, 2.0, -2.5 ),
		simu_cov_range=c( 0, 1),
		simu_t_range = c(-1, 1), 
		debug=F );


tb.phe<-read.csv(phe.out);
tb.snp<-read.csv(snp.out);

rownames(tb.phe) <- tb.phe[,1];
tb.phe <- tb.phe[,-1];

ret1<-ret2<-ret3<-ret4<-ret5<-ret6<-c();

ret1 <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1","X_2"), fgwas.filter = T , options=list(nParallel.cpu=7));	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-snpmat.rdata");
summary(ret1)
plot(ret1, fig.prefix="bls-test-snpmat-ret1");

ret2 <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1","X_2"), fgwas.filter = F );	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-snpmat.rdata");
summary(ret2)
plot(ret2, fig.prefix="bls-test-snpmat-ret2");

ret3 <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1"), fgwas.filter = T, options=list(nParallel.cpu=7) );	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-snpmat.rdata");
summary(ret3)
plot(ret3, fig.prefix="bls-test-snpmat-ret3");

ret4 <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1"), fgwas.filter = T,  refit = F, options=list(nParallel.cpu=7) );	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-snpmat.rdata");
summary(ret4)
plot(ret4, fig.prefix="bls-test-snpmat-ret4");

ret5 <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c(), fgwas.filter = T,  refit = F, options=list(nParallel.cpu=7) );	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-snpmat.rdata");
summary(ret5)
plot(ret5, fig.prefix="bls-test-snpmat-ret5");

ret6 <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1","X_2"), fgwas.filter=T,  refit=T, add.used=T, dom.used=F, options=list(nParallel.cpu=7) );

save(ret1, ret2, ret3, ret4, ret5, ret6, file="bls-test-snpmat.rdata");
summary(ret6)
plot(ret6, fig.prefix="bls-test-snpmat-ret6");

