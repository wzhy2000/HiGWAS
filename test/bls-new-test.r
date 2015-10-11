library(gwas.lasso)
library(parallel)

get_cmd_option<-function(arg.name, def)
{
	v.str <- "";
	args <- commandArgs();
 	for(i in 1:length(args))
 	{
 		if (substr(args[i], 1, nchar(arg.name)) == arg.name )
 		{
 			v.str <- substring(args[i], nchar(arg.name)+2);
 			return(v.str);
 		}
 	}

 	return(def)
}

remove_corralted_snp <-function(snp.mat, simu_sigsnp )
{
	a <- c(1:NROW(snp.mat));
	snp.unsig <- a [! a %in% simu_sigsnp]
	
	snp.rem <- c();
	for(i in simu_sigsnp)
	{
		snp.cors <- unlist( mclapply(snp.unsig, function(j){
			x <- cor( as.numeric(snp.mat[i, -c(1:2)], na.rm=T), as.numeric(snp.mat[j, -c(1:2)]) ); 
			return(x); 
		}, mc.cores=simu_cores ));

show(range(abs(snp.cors)));
			
		snp.rem <- c( snp.rem, snp.unsig[ which(abs(snp.cors)>0.8)] );
	}
	
	if(length(snp.rem)>0)
	{
		show(snp.rem);
		snp.mat <- snp.mat[-unique(snp.rem),];
	}
	
	return(snp.mat);	
}

create_bls_simu<-function()
{
	temp.file <- tempfile();
	phe.out <- paste(temp.file, "phe.out", sep="-");
	snp.out <- paste(temp.file, "snp.out", sep="-");
	
	simu_sigsnp <- sample(1:simu_p)[1:5];
	
	bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= simu_n, simu_p=simu_p,  
			simu_snp_rho   = 0.1, 
			simu_rho       = 0.4, 
			simu_sigma2    = simu_sigma2, 
			simu_snp_missing = 0,
			simu_mu        = 24, 
			simu_cov_range = c( 0, 1),
			simu_cov_effect= c( 0, 2 ), 
			simu_add_pos   = c( simu_sigsnp[1], simu_sigsnp[2], simu_sigsnp[3]), 
			simu_add_effect= c( simu_snp1a, simu_snp2a, simu_snp3a ),  
			simu_dom_pos   = c( simu_sigsnp[3], simu_sigsnp[4], simu_sigsnp[5]), 
			simu_dom_effect= c( simu_snp3d, simu_snp4d, simu_snp5d ),
			simu_t_range   = c( -1, 1), 
			plink.format   = F,
			debug=F );

	tb.snp <- read.csv(snp.out);
	#tb.snp <- remove_corralted_snp(tb.snp, simu_sigsnp );

	tb.phe <- read.csv(phe.out);
	rownames( tb.phe ) <- tb.phe[,1];
	tb.phe <- tb.phe[, -1];

	rownames(tb.snp)[simu_sigsnp] <- paste("Y-SNP", c(1:5), sep="");

	unlink(phe.out);
	unlink(snp.out);

	r.bls  <- bls.snpmat(tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1", "X_2"), refit=F, fgwas.filter=F, options=list(nParallel.cpu=1, nPiecewise.ratio=0) ); 

	sig.detect <- match( paste("Y-SNP", c(1:5), sep=""), rownames(r.bls$varsel_Qbest) );
	show(r.bls$varsel_Qbest[sig.detect,]);
	show(r.bls$varsel[sig.detect,]);

	r.qval <- bls.best.qval( r.bls, paste("Y-SNP", 1:5, sep=""));
	
	cat("N=", simu_n, "P=", simu_p, "SIGMA2=", simu_sigma2, "ADD=", simu_snp1a, simu_snp2a, simu_snp3a, "DOM=", simu_snp3d, simu_snp4d, simu_snp5d, "\n");
	show(r.qval);

	q.add <- max(r.qval[ 1:3 , 2] );
	q.dom <- max(r.qval[ 3:5 , 3] );
	
	cat("Q.val=", q.add, q.dom, "\n");

	r.sig <- bls.qval.cutoff( r.bls, q.add, q.dom );
	show(r.sig);

	return(list(phe.mat=tb.phe, snp.mat=tb.snp, r.bls=r.bls, q.add=q.add, q.dom=q.dom, r.qval=r.qval, r.sig=r.sig));
}

create_simu_snp<-function( maf, n.snp, chr, pos )
{
	allel0 <- ifelse(runif(n.snp)<=maf, 1, 0);
	allel1 <- ifelse(runif(n.snp)<=maf, 1, 0);
	return(c(chr, pos, c(allel0+allel1)) );
}

new_bls_snpmat<-function(phe.mat, snp.mat, grp.size, est.effect=c(2, 1, 0.75, 0.5, 0.1) )
{
	phe.mat0 <- phe.mat;
	
	snp.new <- array(NA ,dim=c( length(est.effect)*2, NCOL(snp.mat) ) );
	for(i in 1:length(est.effect))
	{
		snp.new[i,] <- create_simu_snp( runif(1, 0.3, 0.5), NCOL(snp.mat)-2, 99, i*2-1 );;
		snp.new[length(est.effect)+i,] <- create_simu_snp( runif(1, 0.3, 0.5), NCOL(snp.mat)-2, 100, i*2 );
	}
	
	rownames(snp.new) <- c( paste("Y-ADD", c(1:length(est.effect)), sep=""), paste("Y-DOM", c(1:length(est.effect)), sep="") );
	colnames(snp.new) <- colnames(snp.mat)
	est.effectX <- rep( est.effect , 2 );

	snp.add  <- snp.new[c(1:length(est.effect)),-c(1,2)]-1;
	snp.dom  <- 1- abs(snp.new[c(1:length(est.effect))+length(est.effect),-c(1,2)]-1);
	phe.delt <- colSums(rbind(snp.add, snp.dom)*est.effectX);
	phe.mat$Y <- phe.mat$Y + phe.delt;
	
	loop <- c( 1:(NROW(snp.mat)/grp.size));
	
	r.bls <- mclapply( loop, function(i){
		snp.mat0 <- rbind( snp.mat[c(((i-1)*grp.size+1):(i*grp.size)),], snp.new );
		r.bls0  <- bls.snpmat(phe.mat, snp.mat0, Y.name="Y", covar.names=c("X_1", "X_2"), refit=F, fgwas.filter=F, options=list(nParallel.cpu=1, nPiecewise.ratio=0) ); 
		return(r.bls0);
	}, mc.cores=simu_cores);
	
	qbest <- c();
	qbest.add <- c();
	qbest.dom <- c();
	for(i in loop)
	{
		i.add <- match(  paste("Y-ADD", c(1:length(est.effect)), sep=""), rownames(r.bls[[i]]$varsel_Qbest));
		i.dom <- match(  paste("Y-DOM", c(1:length(est.effect)), sep=""), rownames(r.bls[[i]]$varsel_Qbest));

		qbest.add <- rbind( qbest.add, r.bls[[i]]$varsel_Qbest[i.add, ]);
		qbest.dom <- rbind( qbest.dom, r.bls[[i]]$varsel_Qbest[i.dom, ]);
		qbest <- rbind( qbest, r.bls[[i]]$varsel_Qbest[-c(i.add, i.dom), ]);
	}	

	q.add <- c();
	q.dom <- c();
	for(i in 1:length(est.effect))
	{
		i.add <- which( paste("Y-ADD", i, sep="") == rownames(qbest.add));
		cat( "ADD=[", est.effect[i], "]", min(qbest.add[i.add,1]), mean(qbest.add[i.add,1]), max(qbest.add[i.add,1]), "\n" )
		i.dom <- which( paste("Y-DOM", i, sep="") == rownames(qbest.dom));
		cat( "DOM=[", est.effect[i], "]", min(qbest.dom[i.dom,5]), mean(qbest.dom[i.dom,5]), max(qbest.dom[i.dom,5]), "\n" )
		q.add0 <- min(qbest.add[i.add,1]);
		q.dom0 <- min(qbest.dom[i.dom,5]);
		
		sig.snp <-sig.add<-sig.dom<-c();
		if(q.add0 < 0.1 ) 
		{
			sig.add <- which( qbest[,1]<=q.add0)
			q.add <- c(q.add, q.add0)
		}
		else 
			q.add0 <- NA;

		if(q.dom0 < 0.15) 
		{
			sig.dom <- which( qbest[,5]<=q.dom0)
			q.dom <- q.dom0
		}
		else 
			q.dom0 <-NA;

		sig.snp <- qbest[unique(c( sig.add, sig.dom )),];
		
		sig.detect <- match( paste("Y-SNP", c(1:5), sep=""), rownames(sig.snp) );
		cat("Effect=[", est.effect[i], "] Q=", q.add0, q.dom0, "SIG=", row.names(sig.snp)[sig.detect], "in", NROW(sig.snp), " SNPs \n");
		
	}	
	
	q.add <- min(q.add);
	q.dom <- min(q.dom);
	
	sig.snp <- qbest[unique( which( qbest[,1]<=q.add | qbest[,5]<=q.dom) ),];
	show(rownames(sig.snp));
	sig.detect <- match( paste("Y-SNP", c(1:5), sep=""), rownames(sig.snp) );

	cat("!Final TEST: q.add=", q.add, "q.dom=", q.dom, "Sig.SNP=", row.names(sig.snp)[sig.detect] , "in", NROW(sig.snp), " SNPs \n");

	varsel.snpname <- rownames(sig.snp);

	snp.mat.sig   <- snp.mat[c(varsel.snpname),];
	snp.mat.unsig <- snp.mat[ - match(varsel.snpname, rownames(snp.mat) ),];
	snp.mat.sig.x <- sample( NROW(snp.mat.unsig) )[ 1: c(grp.size - NROW(snp.mat.sig) )]
	#snp.mat0 <- rbind( snp.mat.sig, snp.mat.unsig[ snp.mat.sig.x,] );
	snp.mat0 <- snp.mat.sig;
	
	r.bls  <- bls.snpmat( phe.mat0, snp.mat0, Y.name="Y", covar.names=c("X_1", "X_2"), refit=T, fgwas.filter=F, options=list(nParallel.cpu=1, nPiecewise.ratio=0) ); 
	
	print(summary(r.bls));
	
	sig.detect <- match( paste("Y-SNP", c(1:5), sep=""), rownames(r.bls$varsel_Qbest) );
	show(r.bls$varsel_Qbest[sig.detect,]);
	show(r.bls$varsel[sig.detect,]);
	sig.detect <- match( paste("Y-SNP", c(1:5), sep=""), rownames(r.bls$refit_Qbest) );
	show(r.bls$refit_Qbest[sig.detect,]);
	show(r.bls$refit[sig.detect,]);

	return(r.bls);
}

simu_n      <- as.numeric(get_cmd_option("N", 500));
simu_p      <- as.numeric(get_cmd_option("P", 2000));
simu_sigma2 <- as.numeric(get_cmd_option("SIGMA2",  10));
simu_snp1a  <- as.numeric(get_cmd_option("SNP1.A", 1));
simu_snp2a  <- as.numeric(get_cmd_option("SNP2.A", 1));
simu_snp3a  <- as.numeric(get_cmd_option("SNP3.A", 1));
simu_snp3d  <- as.numeric(get_cmd_option("SNP3.D", 1));
simu_snp4d  <- as.numeric(get_cmd_option("SNP4.D", 1));
simu_snp5d  <- as.numeric(get_cmd_option("SNP5.D", 1));

simu_rdata  <- as.numeric(get_cmd_option("RDATA", 1 ));
simu_loop   <- as.numeric(get_cmd_option("LOOP", 10 ));
simu_cores  <- as.numeric(get_cmd_option("CPU", 5 ));

simu_sigsnp  <- c( 1, 2, 3, 4, 5 );
simu_grpsize <- c(10, 50, 100, 125, 200, 250, 500, 1000, 1500, 2000, 2500, 4000, 5000);
#simu_grpsize <- c(10, 50, 100 );

r <- create_bls_simu();
#x <- new_bls_snpmat(r$phe.mat, r$snp.mat, 250, c(1, 0.5, 0.1) )


