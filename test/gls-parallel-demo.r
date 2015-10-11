library(gwas.lasso)
 
# simu_id=as.integer(get_par("ID"))
# simu_p=as.integer(get_par("P"))
# simu_n=as.integer(get_par("N"))
# simu_sigma2=as.numeric(get_par("SIGMA2"))
# simu_fQval.add = as.numeric(get_par("FQVAL.ADD"))
# simu_fQval.dom = as.numeric(get_par("FQVAL.DOM"))
# simu_nPiecewise.ratio = as.numeric(get_par("NPIECEWISERATIAO"))

 simu_id = 1
 simu_p  = 1000
 simu_n  = 100
 simu_sigma2 = 9
 simu_fQval.add = 0.05
 simu_fQval.dom = 0.05
 simu_nPiecewise.ratio = 1
 
phe.out <- "gls.simple1.10k.phe"  
snp.out <- "gls.simple1.10k.snp"
 
a_effect <- array(0, dim=c(3,4));
d_effect <- array(0, dim=c(3,4));
cov_effect <- array(0, dim=c(2,4));
  
sigsnp <- c(1, 105, 106, 505, 506);

simu_add_pos <- sigsnp[1:3]
simu_dom_pos <- sigsnp[3:5]

a_effect[1,]<-c( 1.04, 0.885, -2.055, 0.545);
a_effect[2,]<-c( 1.17, -0.20, 0.74, -4.715);
a_effect[3,]<-c( 1.40, -2.25, 1.00,  0.00);

d_effect[1,]<-c( 1.49, -2.135, 4.82, 1.425);
d_effect[2,]<-c( 1.045, 1.320, 1.905,  1.535);
d_effect[3,]<-c( 1.265, -1.225, 2.710, -1.96);
 
cov_effect[1,]<-c( 2.49, -1.135, 0.82, 0.425);
cov_effect[2,]<-c( -1.045, 2.320, 0.905,  0.535);
 
start.tm <- proc.time();
gls.simulate( phe.out, snp.out, simu_n= simu_n, simu_grp=1, simu_p=simu_p, simu_snp_rho=0.4, simu_rho=0.1, simu_sigma2= simu_sigma2, 
 		 simu_mu= c(13.395, -3.08, 1.875, -3.195),  simu_cov_effect=cov_effect, simu_cov_range=c(-1,1),
 		 simu_add_effect=a_effect,  simu_dom_effect=d_effect,  simu_add_pos= simu_add_pos, simu_dom_pos=simu_dom_pos,
 		 simu_z_range = c(30,60), simu_z_count = c(5,12), debug=F);
 		 
ret1 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , refit=F, 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=1, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, file="gls-parallel-demo.rdata");

gls.best.qval(ret1, paste("G0-SNP", sigsnp, sep="") );

ret2 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=2, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, file="gls-parallel-demo.rdata");

gls.best.qval(ret2, paste("G0-SNP", sigsnp, sep="") );

ret3 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=3, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, ret3, file="gls-parallel-demo.rdata");

gls.best.qval(ret3, paste("G0-SNP", sigsnp, sep="") );


ret4 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=4, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, ret3, ret4, file="gls-parallel-demo.rdata");

gls.best.qval(ret4, paste("G0-SNP", sigsnp, sep="") );

ret5 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=5, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

gls.best.qval(ret5, paste("G0-SNP", sigsnp, sep="") );

save(ret1, ret2, ret3, ret4, ret5, file="gls-parallel-demo.rdata");


ret6 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=6, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, ret3, ret4, ret5, ret6, file="gls-parallel-demo.rdata");

gls.best.qval(ret6, paste("G0-SNP", sigsnp, sep="") );

ret7 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=7, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, file="gls-parallel-demo.rdata");

gls.best.qval(ret7, paste("G0-SNP", sigsnp, sep="") );

ret8 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=8, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, file="gls-parallel-demo.rdata");

gls.best.qval(ret8, paste("G0-SNP", sigsnp, sep="") );

ret10 <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1", "X_2"), fgwas.filter = F , 
 			options=list(nParallel.cpu=7, nPiecewise.ratio=10, fQval.add = simu_fQval.add, fQval.dom = simu_fQval.dom) );

save(ret1, ret2, ret3, ret4, ret5, ret6, ret7, ret8, ret10, file="gls-parallel-demo.rdata");

gls.best.qval(ret10, paste("G0-SNP", sigsnp, sep="") );
