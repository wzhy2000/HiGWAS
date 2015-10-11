library(gwas.lasso)

phe.out <- "gls.test.qbest.phe"  
snp.out <- "gls.test.qbest.snp"
snp.out.bed <- "gls.qbest.simple.snp.bed"
snp.out.bim <- "gls.qbest.simple.snp.bim"
snp.out.fam <- "gls.qbest.simple.snp.fam"

a_effect <- array(0, dim=c(3,4));
d_effect <- array(0, dim=c(3,4));
cov_effect <- array(0, dim=c(2,4));

sigsnp <- sample(1:30)[1:5];
show(sigsnp);

a_effect[1,]<-c( 1.04, 0.885, -2.055, 0.545);
a_effect[2,]<-c( 1.17, -0.20, 0.74, -4.715);
a_effect[3,]<-c( 1.40, -2.25, 1.00,  0.00);

d_effect[1,]<-c( 1.49, -2.135, 4.82, 1.425);
d_effect[2,]<-c( 1.045, 1.320, 1.905,  1.535);
d_effect[3,]<-c( 1.265, -1.225, 2.710, -1.96);

cov_effect[1,]<-c( 2.49, -1.135, 0.82, 0.425);
cov_effect[2,]<-c( -1.045, 2.320, 0.905,  0.535);


gls.simulate( phe.out, snp.out, simu_n= 600, simu_grp=1, simu_p=30, simu_snp_rho=0.4, simu_rho=0.1, simu_sigma2= 4, 
		 simu_mu= c(13.395, -3.08, 1.875, -3.195),  
		 simu_cov_effect = cov_effect, 
		 simu_cov_range  = c(-1,1),
		 simu_add_pos    = c( sigsnp[1], sigsnp[2], sigsnp[3] ),
		 simu_add_effect = a_effect,  
		 simu_dom_pos    = c( sigsnp[3], sigsnp[4], sigsnp[5] ),
		 simu_dom_effect = d_effect, 
		 simu_z_range    = c(30,60), simu_z_count = c(5,12), debug=F);

ret1 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter = F,  options=list(nPiecewise.ratio=0)  );	

save(ret1, file="gls-test-qbest.rdata");

summary(ret1);
r.qval <- gls.best.qval(ret1, paste("G0-SNP", sigsnp, sep="") );
show(r.qval);

r.qval <- gls.best.qval(ret1, paste("G0-SNP", sigsnp, sep="") );
show(r.qval);

r.qval <- gls.best.qval(ret1, paste("G0-SNP", c(sigsnp, 20, 30, 10), sep="") );
show(r.qval);

r.qval <- gls.best.qval(ret1, paste("G0-SNP", c(sigsnp, 2000, 3000, 1000), sep="") );
show(r.qval);

r.qval <- gls.best.qval(ret1, paste("G0-SNP", c(3000), sep="") );
show(r.qval);

vs.names <- gls.qval.cutoff(ret1, 0.05, 0.10, refit.select=FALSE );
show(vs.names);

refit.names <- gls.qval.cutoff(ret1, 0.05, 0.10, refit.select=TRUE );
show(refit.names);
