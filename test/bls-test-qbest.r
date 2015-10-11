library(gwas.lasso)

phe.out <- "bls.test.qbest.phe"  
snp.out <- "bls.test.qbest.snp"
snp.out.bed <- "bls.test.qbest.snp.bed"
snp.out.bim <- "bls.test.qbest.snp.bim"
snp.out.fam <- "bls.test.qbest.snp.fam"


sigsnp <- sample(1:100)[1:5];
show(sigsnp);

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

ret1 <- bls.simple( phe.out, snp.out, Y.name="Y", covar.names=c("X_1", "X_2"), fgwas.filter=F, options=list(nPiecewise.ratio=0) );	

save(ret1, file="bls-test-qbest.rdata");

summary(ret1);
r.qval <- bls.best.qval(ret1, paste("G0-SNP", sigsnp, sep="") );
show(r.qval);

r.qval <- bls.best.qval(ret1, paste("G0-SNP", c(sigsnp, 20, 30, 10), sep="") );
show(r.qval);

r.qval <- bls.best.qval(ret1, paste("G0-SNP", c(sigsnp, 2000, 3000, 1000), sep="") );
show(r.qval);

r.qval <- bls.best.qval(ret1, paste("G0-SNP", c(3000), sep="") );
show(r.qval);

snp.names <- bls.qval.cutoff(ret1, 0.05, 0.15, refit.select=FALSE);
