library(gwas.lasso)

phe.out <- "gls.test.simple.phe"  
snp.out <- "gls.test.simple.snp"

a_effect <- array(0, dim=c(3,4));
d_effect <- array(0, dim=c(3,4));
cov_effect <- array(0, dim=c(2,4));

sigsnp <- sample(1:100)[1:5];
show(sigsnp);

a_effect[1,]<-c( 1.04, 0.885, -2.055, 0.545);
a_effect[2,]<-c( 1.17, -0.20, 0.74, -4.715);
a_effect[3,]<-c( 1.40, -2.25, 1.00,  0.00);

d_effect[1,]<-c( 1.49, -2.135, 4.82, 1.425);
d_effect[2,]<-c( 1.045, 1.320, 1.905,  1.535);
d_effect[3,]<-c( 1.265, -1.225, 2.710, -1.96);

cov_effect[1,]<-c( 2.49, -1.135, 0.82, 0.425);
cov_effect[2,]<-c( -1.045, 2.320, 0.905,  0.535);

simu.data <- gls.simulate( phe.out, snp.out, simu_n= 100, simu_grp=1, simu_p=200, simu_snp_rho=0.4, simu_rho=0.1, simu_sigma2= 4, 
		 simu_mu= c(13.395, -3.08, 1.875, -3.195),  
		 simu_cov_effect = cov_effect, 
		 simu_cov_range  = c(-1,1),
		 simu_add_pos    = c( sigsnp[1], sigsnp[2], sigsnp[3] ),
		 simu_add_effect = a_effect,  
		 simu_dom_pos    = c( sigsnp[3], sigsnp[4], sigsnp[5] ),
		 simu_dom_effect = d_effect, 
		 simu_z_range    = c(30,60), simu_z_count = c(5,12),plink.format=T, debug=F);
		 
if(simu.data$err==0)
{
	file.plink.bed <- simu.data$file.plink.bed  
	file.plink.bim <- simu.data$file.plink.bim 
	file.plink.fam <- simu.data$file.plink.fam 
	
}

ret1 <- gls.plink( phe.out, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = T , options=list(nParallel.cpu=1) );


