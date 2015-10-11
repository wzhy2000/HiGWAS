plink_fgwas_bigdata <- function( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe , plink.command , 
	    						          Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS")
{
	cat( "Big PLINK dataset is splitted by PLINK command......\n");

	tb.bim <- read.table(file.plink.bim, header=F);
	tb.fam <- read.table(file.plink.fam, header=F);

	chrs <- unique(tb.bim[,1]);
	
	cat("*", length(chrs), "chromosomes/groups", NROW(tb.fam), "individuals and", NROW(tb.bim), "SNPs are stored in the PLINK dataset.\n")
	
	if(is.null(plink.command)) plink.command <- "plink";
	t <- try(system( paste( plink.command, "--noweb", sep=" "), ignore.stdout=TRUE, ignore.stderr=TRUE ));
	
	if(class(t)=="try-error")
	{
		cat("! No PLINK command can be found in your environment( plink.command=",plink.command, ")\n")
		return(list(error=T, err.info="No PLINK command."));
	}

	if(length(chrs)<=1)
	{
		cat("! No more chromosomes need to split.\n")
		return(list(error=T, err.info="No more chromosomes need to split."));
	}
	
	cat( "SNP Filtering by fGWAS method......\n");
	cat("* p-value Threshold :", fgwas.cutoff, "\n" );
	
	snp.mat <- c();
	phe.mat <- NULL;
	r.fgwas <- c();
	
	for(chr in chrs)
	{
		# "PLINK --chr 0" will output all SNPs, no extraction.
		if (as.character(chr) == "0" ) next;
		
		cat(" ", chr, "(th) chromosome/group is being proceeded.\n");	

		tmp <- tempfile();
		str.cmd <- paste( plink.command, "--noweb", 
					"--bed", file.plink.bed,
					"--bim", file.plink.bim,
					"--fam", file.plink.fam,
					"--chr", chr,
					"--make-bed",
					"--out", tmp, 
					sep=" ");
					
		t <- try(system( str.cmd, ignore.stdout=TRUE, ignore.stderr=TRUE) );
		if(class(t)=="try-error")			
		{
			cat("! Error in PLINK command.\n")
			return(list(error=T, err.info="Error in PLINK command."));
		}

		chr.pd <- load_plink_binary( paste(tmp, ".bed", sep=""),  paste(tmp, ".bim", sep=""), paste(tmp, ".fam", sep=""), file.phe );
		if( is.null(chr.pd) )
		{
			cat("! Package snpStats can not open PLINK data.\n")
			return(list(error=T, err.info="Package snpStats can not open PLINK data."));
		}

		chr.filter <- plink_fgwas_filter( chr.pd, Y.name, Z.name, covar.names, op.cpu=op.cpu, fgwas.cutoff=fgwas.cutoff, lasso_method=lasso_method, verbose=FALSE, only.sig.snp=TRUE );
		if( chr.filter$error )
			return(chr.filter);
			
		cat(" ", NROW(chr.filter$snp.mat), "SNPs are significant at the pvalue level in ", chr, "(th) chromosome/group.\n");	
		
		phe.mat <- chr.pd$phe.mat;
		snp.mat <- rbind(snp.mat, chr.filter$snp.mat);
		r.fgwas <- rbind(r.fgwas, chr.filter$r.fgwas);
		
		try(unlink(paste(tmp, ".bed", sep="")));
		try(unlink(paste(tmp, ".bim", sep="")));
		try(unlink(paste(tmp, ".fam", sep="")));
	}

	cat("*", NROW(snp.mat), " SNPs are detected by fGWAS method.\n")
	
	return(list(error=F, snp.mat=snp.mat, r.fgwas=r.fgwas, phe.mat=phe.mat) );
}	    						          

plink_fgwas_filter<-function( pd, Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS", verbose=TRUE, only.sig.snp=FALSE)
{
	n.snp <- NCOL( pd$snp.mat$genotypes );
	n.ind <- NROW( pd$snp.mat$genotypes );
		
	get_sub_snpmat<- function(pd.obj, idx.snp)
	{
		snp.sub <- get_plink_subsnp( pd.obj$snp.mat, idx.snp );

		# Append Chr and pos. information to SNP.MAT
		snp.mat <- cbind( snp.sub$info[,c(2,3)], snp.sub$snp )
		return(snp.mat);
	}
	
	r <- fgwas_filter( n.snp, n.ind, get_sub_snpmat, pd, pd$phe.mat, Y.name, Z.name, covar.names, op.cpu, fgwas.cutoff, lasso_method, verbose=verbose, only.sig.snp=only.sig.snp );
	return(r);
}

snpmat_fgwas_filter<-function( phe.mat, snp.mat, Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS", verbose=TRUE, only.sig.snp=FALSE)
{
	n.snp <- NROW( snp.mat );
	n.ind <- NCOL( snp.mat ) -2 ;
		
	get_sub_snpmat<- function(snpmat.big, idx.snp)
	{
		return(snpmat.big[idx.snp,,drop=F]);
	}
	
	r <- fgwas_filter( n.snp, n.ind, get_sub_snpmat, snp.mat, phe.mat, Y.name, Z.name, covar.names, op.cpu, fgwas.cutoff, lasso_method, verbose=verbose, only.sig.snp=only.sig.snp );
	return(r);
}

fgwas_filter<-function( n.snp, n.ind, f_get_snpmat, snp.obj, phe.mat, Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS", verbose=TRUE, only.sig.snp=FALSE, include.na.pvalue=TRUE )
{
	if(verbose)
	{
		cat( "SNP Filtering by fGWAS method......\n");
		cat("* SNP Count =", n.snp, "\n" );
		cat("* Sample Count =", n.ind, "\n" );
		cat("* p-value Threshold =", fgwas.cutoff, "\n" );
	}
	
	snp.sect0 <- seq( 1, n.snp, 10000 );
	snp.sect1 <- c( snp.sect0[-1]-1, n.snp );
	
	if( class(phe.mat)=="matrix") phe.mat <- as.data.frame(phe.mat);

	r.fgwas <- c();
	snp.mat <- c();
	
	for(i in 1:length(snp.sect0))
	{
		snp.mat0 <- f_get_snpmat( snp.obj, snp.sect0[i]:snp.sect1[i] );

		if(verbose)	cat("  Calculated SNP Range =", snp.sect0[i], snp.sect1[i], "\n" );

		r.fgwas0 <- list();
		if(lasso_method=="BLS")
			r.fgwas0 <- bls.fgwas( phe.mat, snp.mat0, Y.name, covar.names, op.cpu )
		else
			r.fgwas0 <- gls.fgwas( phe.mat, snp.mat0, Y.name, Z.name, covar.names, op.cpu )
		
		if( r.fgwas0$error )
			return(r.fgwas0);

		#adjust SNP.IDX field.
		r.fgwas0$r[,1] <- r.fgwas0$r[,1] + snp.sect0[i]-1;

		r.fgwas <- rbind( r.fgwas, r.fgwas0$r);
	}
	
	
	r.filter0 <- get_sigsnp_nomulti_correction( f_get_snpmat, snp.obj, r.fgwas, n.snp, n.ind, fgwas.cutoff, only.sig.snp=only.sig.snp, include.na.pvalue=include.na.pvalue );
	if( r.filter0$error )
		stop( r.filter0$err.info );

	snp.mat <- r.filter0$snp.mat;
	
	f.per <- round( NROW(snp.mat)/n.snp*100, digits=2);
	if(verbose) cat("*", NROW(snp.mat), "SNPs(", f.per ,"%) are detected by the fGWAS method.\n" );

	return(list(error=F, snp.mat=snp.mat, r.fgwas=r.fgwas) );
}


gls.fgwas <- function( phe.mat, snp.mat, Y.prefix, Z.prefix, covar.names=NULL, op.cpu=0)
{
	sample.ids <- intersect( rownames(phe.mat), colnames(snp.mat)[-c(1:2)] );
	if(length(sample.ids)==0)
	{
		cat("! No same population data in the phenotypic data and genotypic data.\n");
		return(list(error=T, err.info="No same population data in the phenotypic data and genotypic data."))
	}
	
	phe.mat <- phe.mat[ sample.ids, , drop=F]
	snp.mat <- cbind( snp.mat[,c(1,2)], snp.mat[ , sample.ids, drop=F]);
	
	y.p <- grep( Y.prefix, colnames(phe.mat) ); 
	z.p <- grep( Z.prefix, colnames(phe.mat) );
	x.p <- c(); 
	if(!is.null(covar.names))
		x.p <- match( covar.names, colnames(phe.mat) );

	Y.sub <- substring( colnames(phe.mat)[y.p], nchar(Y.prefix)+1);
	Z.names <- paste( Z.prefix, Y.sub, sep="");
	z.p <- match( Z.names, colnames(phe.mat) );
	if( length( which(is.na(z.p)) )>0 )
	{
		cat("! The Z columns are not matched with Y columns.\n");
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))
	}
	
	len.x <- length(x.p);
	len.y <- length(y.p);
	len.z <- length(z.p);
	if ( len.y != len.z)
	{
		cat("! The Z columns are not matched with Y columns.\n");
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))
	}
	
	if( length(x.p)>0 && length(which(is.na(x.p)))>0)
	{
		cat("! The covariate names are not matched with phenotypical data.\n");
		return(list(error=T, err.info="The covariate names are not matched with phenotypical data."))
	}
	
	phe.mat <- cbind(phe.mat, ID = c(1:NROW(phe.mat)) );
	id.p <- NCOL(phe.mat);

	phe.gls.mat <- array(0,dim=c(0, 1 + len.x + 2 ))
	for(i in 1:len.y)
	{
		phe.tmp <- phe.mat[, c( ID=id.p, Y=y.p[i], Z=z.p[i], x.p), drop=F ];
		
		phe.tmp.col <- colnames(phe.tmp);
		phe.tmp.col[ c(1:3) ] <- c("ID", "Y", "Z");
		colnames(phe.tmp) <- phe.tmp.col;
		
		phe.tmp.missing <- which( is.na(phe.tmp$Z) | is.na(phe.tmp$Y) );
		if(length(phe.tmp.missing)>0)
			phe.tmp <- phe.tmp[-phe.tmp.missing, , drop=F ];
		phe.gls.mat <- rbind( phe.gls.mat, phe.tmp );
	}

	reg.str0 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+")),  "Y ~ 1" );
	reg.str1 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+"), "+ as.factor(SNP)") ,  "Y ~ 1 + as.factor(SNP)" );

	cat("* H0 =", as.character(reg.str0), "\n" );
	cat("* H1 =", as.character(reg.str1), "\n" );

	r0 <- try( gls( as.formula(reg.str0), phe.gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML" ) );	
	if(class(r0)=="try-error")
	{
		cat("! Failed to call gls() method.\n");
		return(list( error=T, err.info="Failed to call gls() method.") )
	}
	
	cpu.fun<-function( sect )
	{
		if( (sect-1)*n.percpu+1 > NROW(snp.mat) )
			return(NULL);

		range.fr <- (sect-1)*n.percpu+1;
		range.to <- sect*n.percpu;
		range.to <- ifelse( range.to > NROW(snp.mat),  NROW(snp.mat), range.to );

		r.gls <- array( NA, dim = c( (range.to-range.fr+1), 7 ) );
		r.gls[,1]<-c(range.fr:range.to);
	
		reg.f0 <- as.formula(reg.str0);
		reg.f1 <- as.formula(reg.str1);

		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(3:NCOL(snp.mat)) ] ); 
			snp.gls <- snp[ phe.gls.mat$ID ];
			gls.mat <- cbind( phe.gls.mat, SNP=snp.gls );
	
			maf <- mean(snp.gls, na.rm=T)/2;
			if(maf>0.5) maf <- 1-maf;

			pv.max <- NA;
			if( maf <= 10^(-4)) pv.max <- 1.0;
	
			# remove NA row
			na.row <- unique( which(is.na(gls.mat))%% NROW(gls.mat) );
			if(length(na.row)>0) gls.mat <- gls.mat[-na.row,,drop=F];
	
			r0 <- try( do.call("gls", args = list(reg.f0, gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML") ) );
			r1 <- try( do.call("gls", args = list(reg.f1, gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML") ) );
			if(any(class(r1)=="try-error") || any(class(r0)=="try-error") )
			{
				r.gls[i-range.fr+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, NA, pv.max );
				next;	
			}

			r <- try( anova( r0,r1 ) );                   
			if(any(class(r)=="try-error"))
			{
				r.gls[i-range.fr+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf,NA, pv.max );
				next;	
			}

			r.gls[(i-range.fr+1),] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, r[2,8], r[2,9]);
		}
		
		return(r.gls);
	}
	

	r.gls <- c();
	if( op.cpu>1 && require("snowfall") )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		snowfall::sfInit(parallel = TRUE, cpus = op.cpu, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / op.cpu );
		snowfall::sfExport("n.percpu", "phe.gls.mat", "snp.mat", "covar.names", "reg.str0", "reg.str1" );
		
		gls.cluster <- snowfall::sfClusterApplyLB( 1:op.cpu, cpu.fun);
		snowfall::sfStop();

		cat("Stopping parallel computing......\n");
		
		r.gls <- do.call( rbind, gls.cluster );
	}		
	else
	{
		cat("Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);	
		r.gls <- cpu.fun(1);
	}

	colnames(r.gls) <- c("SNP.IDX", "CHR", "POS", "N.miss", "MAF", "L.Ratio", "pv");
	rownames(r.gls) <- rownames(snp.mat) 
	
	return(list(error=F, r=r.gls));
}

bls.fgwas <- function( phe.mat, snp.mat, Y.name, covar.names=NULL, op.cpu=0)
{
	sample.ids <- intersect( rownames(phe.mat), colnames(snp.mat)[-c(1:2)] );
	if(length(sample.ids)==0)
	{
		cat("! No same population data in the phenotypic data and genotypic data.\n")
		return(list(error=T, err.info="No same population data in the phenotypic data and genotypic data."))
	}
	
	phe.mat <- phe.mat[ sample.ids, c(Y.name, covar.names), drop=F]
	phe.mat <- cbind(ID = c(1:NROW(phe.mat)), phe.mat);
	snp.mat <- cbind( snp.mat[,c(1,2)], snp.mat[, sample.ids, drop=F]);
	
	len.x <- 0;
	x.p   <- c();
	if( !is.null( covar.names ) )
	{
		x.p <- match( covar.names, colnames(phe.mat) );
		len.x <- length(x.p);
	}
	
	phe.missing <- which( is.na(phe.mat[, dim(phe.mat)[2]]));
	if(length(phe.missing)>0)
		phe.mat <- phe.mat[-phe.missing, , drop=F ];

	str01 <- paste(Y.name, "~", paste(covar.names,collapse= "+") );
	str00 <- paste(Y.name,"~ 1");
	reg.str0 <- ifelse( len.x>0, str01, str00 );

	cat("* H0 =", as.character(reg.str0), "\n" );

	str01 <- paste(Y.name, "~", paste(covar.names,collapse= "+"), "+ as.factor(SNP)") ;
	str00 <- paste(Y.name, "~ 1 + as.factor(SNP)" );
	reg.str1 <- ifelse( len.x>0, str01, str00 );

	cat("* H1 =", as.character(reg.str1), "\n" );
	
	r0 <- try( do.call("gls", args = list(as.formula(reg.str0), phe.mat, method="ML" ) ) );
	if(class(r0)=="try-error")
	{
		cat("! Failed to call gls() method.\n")
		return(list( error=T, err.info="Failed to call gls() method.") )
	}
	
	cpu.fun<-function( sect )
	{
		if( (sect-1)*n.percpu+1 > NROW(snp.mat) )
			return(NULL);

		range.fr <- (sect-1)*n.percpu+1;
		range.to <- sect*n.percpu;
		range.to <- ifelse( range.to > NROW(snp.mat),  NROW(snp.mat), range.to );

		r.bls <- array( NA, dim = c( (range.to-range.fr+1), 7 ) );
		r.bls[,1]<-c(range.fr:range.to);
	
		reg.f0 <- as.formula(reg.str0);
		reg.f1 <- as.formula(reg.str1);
		
		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(3:NCOL(snp.mat)) ] ); 
			snp.bls <- snp[ phe.mat$ID ];
			bls.mat <- cbind( phe.mat, SNP=snp.bls );
			
			maf <- mean(snp.bls, na.rm=T)/2;
			if(maf>0.5) maf <- 1-maf;

			# if maf closes to 0, it will lead to an error: contrasts not defined for 0 degrees of freedom
			pv.max <- NA;
			if( maf <= 10^(-4)) pv.max <- 1.0;
			
			# remove NA row
			na.row <- unique( which(is.na(bls.mat))%% NROW(bls.mat) );
			if(length(na.row)>0) bls.mat <- bls.mat[-na.row,,drop=F];

			r0 <- try( do.call("gls", args = list( reg.f0, bls.mat, method="ML" ) ) );
			r1 <- try( do.call("gls", args = list( reg.f1, bls.mat, method="ML" ) ) );
			
			if(any(class(r0)=="try-error") || any(class(r1)=="try-error"))
			{
				r.bls[i-range.fr+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, NA, pv.max );
				next;	
			}

			r <- try( anova( r0,r1 ) );                   
			if(any(class(r)=="try-error"))
			{
				r.bls[i-range.fr+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, NA, pv.max );
				next;	
			}

			r.bls[(i-range.fr+1),] <- c( i, snp.mat[i,1], snp.mat[i,2], length(na.row), maf, r[2,8], r[2,9]);
		}
		
		return(r.bls);
	}
	
	r.bls <- c();
	if( op.cpu>1 && require("snowfall") )
	{
		cat("\n  Starting parallel computing, snowfall/snow......\n"); 
		snowfall::sfInit(parallel = TRUE, cpus = op.cpu, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / op.cpu );
		snowfall::sfExport("n.percpu", "phe.mat", "snp.mat", "covar.names", "reg.str0", "reg.str1" );
		
		bls.cluster <- snowfall::sfClusterApplyLB( 1:op.cpu, cpu.fun);
		snowfall::sfStop();

		cat("  Stopping parallel computing......\n\n");

		r.bls <- do.call( rbind, bls.cluster );

	}		
	else
	{
		cat("  Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);	
		r.bls <- cpu.fun(1);
	}
	
	colnames(r.bls) <- c("SNP.IDX", "CHR", "POS", "N.miss", "MAF", "L.Ratio", "pv");
	rownames(r.bls) <- rownames(snp.mat) 
	
	return(list(error=F, r=r.bls));
	
}

get_sigsnp_nomulti_correction<-function( f_get_snpmat, snp.obj, r.fgwas, n.snp, n.ind, fgwas.cutoff=0.05, only.sig.snp=FALSE, include.na.pvalue=TRUE )
{
	if( n.snp <= n.ind && !only.sig.snp)
	{
		snp.mat <- f_get_snpmat( snp.obj, 1:n.snp );
		return(list(error=F, snp.mat=snp.mat ));
	}	
	
	pv <- r.fgwas[,7];

	# in order to keep NA element in the vector , replace NA with Inf.
	pv [ is.na(pv) ] <- Inf;
	
	n.sig <- length( which( pv <= fgwas.cutoff ) );
	cat("  SNPs with p-value <=", fgwas.cutoff, ":", n.sig, "\n"); 
	
	#sorting pv field
	pv.sort  <- sort.int( pv, decreasing=F, index.return=T);

	sel.p <- c();	
	if ( only.sig.snp )
	{
		## Here is (n.snp > n.ind) && (only.sig.snp)
		if ( length( which( pv.sort$x <= fgwas.cutoff) ) >0)
		{
			max.idx <- max( which( pv.sort$x <= fgwas.cutoff) );
			sel.p <- pv.sort$ix[1:max.idx];
		}
	}
	else
	{
		## Here is (n.snp > n.ind) && (!only.sig.snp)
		if( n.sig < n.ind )
		{
			cat("  The count of significant SNPs is less than individual count.\n"); 
			sel.p <- pv.sort$ix[1:n.ind];
		}
		else
		{
			## Assert max.idx > 0
			max.idx <- max( which( pv.sort$x <= fgwas.cutoff) );
			sel.p <- pv.sort$ix[1:max.idx];
		}
	}
	
	if( include.na.pvalue )
		sel.p<- c( sel.p, which(is.na(r.fgwas[,7])) ) ;
	
	# if NO SNPs are slected.
	if( length(sel.p)==0 )
		return(list(error=F, snp.mat=NULL));
	
	sel.id <- r.fgwas[ sort(sel.p), 1 ];

	snp.mat <- f_get_snpmat( snp.obj, sel.id );

	return(list(error=F, snp.mat=snp.mat));
}


fgwas.scan<-function( phe.mat, snp.mat, Y.name, Z.name=NULL, covar.names, options=list(nParallel.cpu=1), longitudinal=FALSE)
{
	n.snp <- NROW(snp.mat);
	n.ind <- NROW(phe.mat);
	
	cat("fGWAS Scanning......\n");
	cat("* SNP Count =", n.snp, "\n" );
	cat("* Sample Count =", n.ind, "\n" );
	
	if( class(phe.mat)=="matrix") phe.mat <- as.data.frame(phe.mat);
	
	ind.snpmat <- colnames(snp.mat)[-c(1,2)];
	ind.phemat <- rownames(phe.mat);
	ind.sames  <- intersect( ind.phemat, ind.snpmat );
	
	if(length(ind.sames)==0)
	{
		cat("! No same individuals' names in SNP matrix and phenotypic matrix.\n")
		return(NULL);	
	}

	if( length(ind.sames) < length(ind.phemat) || length(ind.sames) < length(ind.snpmat) )
	{
		cat("! Different size of individuals in SNP matrix and phenotypic matrix. snp.mat=", length(ind.snpmat), " phe.mat=", length(ind.phemat), "\n")
		
		phe.mat <- phe.mat[ match( ind.sames, rownames(phe.mat)),,drop=F]
		snp.mat <- snp.mat[,c(1,2, match( ind.sames, colnames(snp.mat))),drop=F]
	}

	r.fgwas0 <- list();
	if ( longitudinal )
		r.fgwas0 <- gls.fgwas( phe.mat, snp.mat, Y.name, Z.name, covar.names, op.cpu=options$nParallel.cpu )
	else
		r.fgwas0 <- bls.fgwas( phe.mat, snp.mat, Y.name, covar.names, op.cpu=options$nParallel.cpu )
		
	if( r.fgwas0$error )
		stop(r.fgwas0$err.info)

	r <- r.fgwas0$r;
	class(r)<-"fgwas.ret0";
	return(r);
}


plot.fgwas.ret0 <- function( x, y=NULL, ..., fig.prefix=NULL )
{
	r.fgwas <- x;
	
	if( missing(fig.prefix)) fig.prefix <- "fgwas.plot";
	
	filter.man <- r.fgwas[, c(1,2,7), drop=F]
	draw_man_fgwas( filter.man, fig.prefix, "fgwas" );

	invisible(fig.prefix);
}

summary.fgwas.ret0 <- function(object, ..., fgwas.cutoff=0.05)
{
	r.fgwas <- object;

	fgwas.sig <- which( r.fgwas[,7] <= fgwas.cutoff );
	cat(length(fgwas.sig), "SNPs at the significant level:",  fgwas.cutoff, "\n")

	if(length(fgwas.sig)>0)
	{
		fgwas_sigs <- r.fgwas[ fgwas.sig, , drop=F];
		fgwas.sig.inc <- order( fgwas_sigs[,7] );
		
		cat("Top 20 SNPs:\n");
		show( head( fgwas_sigs[fgwas.sig.inc,],20 ) );
	}	
	
	invisible();
}

find_fgwas_pvalue<-function( fgwas.ret, snp.names)
{
	pv <- rep(NA, length(snp.names));

	idx <- match( snp.names, rownames(fgwas.ret));
	if( length(which(!is.na(idx)) ) > 0)
		pv[ which(!is.na(idx)) ] <-  fgwas.ret[ idx[!is.na(idx)], 7 ];
	
	return(pv);
}