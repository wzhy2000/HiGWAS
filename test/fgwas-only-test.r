library(gwas.lasso)

phe.out <- "bls.test.fgwas.phe"  
snp.out <- "bls.test.fgwas.snp"

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

ret1 <- fgwas.scan(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1","X_2"), options=list(nParallel.cpu=1), longitudinal=FALSE);	
save(ret1, file="fgwas-bls-test.rdata");
summary(ret1)
plot(ret1, fig.prefix="fgwas-bls-test");


phe.out <- "gls.test.fgwas.phe"  
snp.out <- "gls.test.fgwas.snp"

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

ret2 <- fgwas.scan( tb.phe, tb.snp, Y.nam="Y", Z.name="Z", covar.names=c("X_1","X_2"), options=list(nParallel.cpu=0), longitudinal=TRUE);
save(ret2, file="fgwas-gls-test.rdata");
summary(ret2);
plot(ret2, fig.prefix="fgwas-gls-test");