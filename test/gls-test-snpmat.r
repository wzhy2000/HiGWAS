library(gwas.lasso)


phe.out <- "gls.test.snpmat.phe"  
snp.out <- "gls.test.snpmat.snp"

a_effect <- array(0, dim=c(3,4));
d_effect <- array(0, dim=c(3,4));
cov_effect <- array(0, dim=c(2,4));

sigsnp <- sample(1:4000)[1:5];
show(sigsnp);

a_effect[1,]<-c( 1.04, 0.885, -2.055, 0.545);
a_effect[2,]<-c( 1.17, -0.20, 0.74, -4.715);
a_effect[3,]<-c( 1.40, -2.25, 1.00,  0.00);

d_effect[1,]<-c( 1.49, -2.135, 4.82, 1.425);
d_effect[2,]<-c( 1.045, 1.320, 1.905,  1.535);
d_effect[3,]<-c( 1.265, -1.225, 2.710, -1.96);

cov_effect[1,]<-c( 2.49, -1.135, 0.82, 0.425);
cov_effect[2,]<-c( -1.045, 2.320, 0.905,  0.535);

gls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 800, simu_p=4000, simu_snp_rho=0.4, simu_rho=0.1, simu_sigma2= 4, 
		 simu_mu         = c(13.395, -3.08, 1.875, -3.195),  
		 simu_cov_effect = cov_effect, 
		 simu_cov_range  = c(-1,1),
		 simu_add_pos    = c( sigsnp[1], sigsnp[2], sigsnp[3] ),
		 simu_add_effect = a_effect,  
		 simu_dom_pos    = c( sigsnp[3], sigsnp[4], sigsnp[5] ),
		 simu_dom_effect = d_effect, 
		 simu_z_range    = c(30,60), simu_z_count = c(5,12), debug=F);


tb.phe<-read.csv(phe.out);
tb.snp<-read.csv(snp.out);

rownames(tb.phe) <- tb.phe[,1];
tb.phe <- tb.phe[,-1];


ret1<-ret2<-ret3<-ret4<-ret5<-ret6<-ret7<-ret8<-c();

ret1 <- gls.snpmat(tb.phe, tb.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter = T , options=list(nParallel.cpu=7));	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-test-snpmat.rdata");
summary(ret1);
plot(ret1, fig.prefix="gls-test-snpmat-ret1");

ret2 <- gls.snpmat(tb.phe, tb.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter = F );	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-test-snpmat.rdata");
summary(ret2);
plot(ret2, fig.prefix="gls-test-snpmat-ret2");

ret3 <- gls.snpmat(tb.phe, tb.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1"), fgwas.filter = T , options=list(nParallel.cpu=7));	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-test-snpmat.rdata");
summary(ret3);
plot(ret3, fig.prefix="gls-test-snpmat-ret3");

ret4 <- gls.snpmat(tb.phe, tb.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1"), fgwas.filter = T,  refit = F, options=list(nParallel.cpu=7) );	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-test-snpmat.rdata");
summary(ret4);
plot(ret4, fig.prefix="gls-test-snpmat-ret4");

ret5 <- gls.snpmat(tb.phe, tb.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c(), fgwas.filter = T,  refit = F , options=list(nParallel.cpu=7));	

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-test-snpmat.rdata");
summary(ret5);
plot(ret5, fig.prefix="gls-test-snpmat-ret5");

ret6 <- gls.snpmat(tb.phe, tb.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter=T,  refit=F, add.used=T, dom.used=F, options=list(nParallel.cpu=7) );

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-test-snpmat.rdata");
summary(ret6);
plot(ret6, fig.prefix="gls-test-snpmat-ret6");
