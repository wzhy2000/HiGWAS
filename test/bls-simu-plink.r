library(gwas.lasso)

phe.out <- "bls.test.simple.phe"  
snp.out <- "bls.test.simple.snp"

sigsnp <- sample(1:2000)[1:5];
show(sigsnp);

simu.data <- bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 1000, simu_p=2000, 
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
		simu_t_range = c(-1, 1), plink.format=T,
		debug=F );

if(simu.data$err==0)
{
	file.plink.bed <- simu.data$file.plink.bed;
	file.plink.bim <- simu.data$file.plink.bim;
	file.plink.fam <- simu.data$file.plink.fam;
	#file.phe.long  <- "phe.mean.csv";	
}

ret1<-ret2<-ret3<-ret4<-ret5<-ret6<-c();
ret1 <- bls.plink( phe.out, file.plink.bed, file.plink.bim, file.plink.fam, Y.name="Y", covar.names=c(),    fgwas.filter = T,  refit = F , options=list(nParallel.cpu=1));

save(ret1, ret2, ret3, ret4, ret5, ret6, sigsnp, file="bls-test-simple.rdata");
summary(ret1)
plot(ret1, fig.prefix="bls-test-simple-ret1");


