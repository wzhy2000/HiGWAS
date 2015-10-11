bls.simulate<-function( file.phe.out, file.snp.out, simu_grp=1, simu_n= 500, simu_p=1000, 
		simu_snp_rho = 0.1, 
		simu_snp_missing = 0.002, 
		simu_rho     = 0.4, 
		simu_sigma2  = 3, 
		simu_mu      = 26, 
		simu_cov_range = c( 0, 1),
		simu_cov_effect= c( 0, 2 ), 
		simu_add_pos   = c( 100, 200, 300), 
		simu_add_effect= c( 2.2, -2.5, 2.0 ),  
		simu_dom_pos   = c( 300, 500, 700), 
		simu_dom_effect= c( 2.8, 2.0, -2.5 ),
		simu_t_range = c(-1, 1), 
		plink.format = FALSE,
		debug=FALSE )
{
	if( !missing(simu_grp) && length(simu_grp) > 1)
		stop("! The parameter of simu_grp is not a single valid value.");

	if( !missing(simu_n) && length(simu_n) > 1)
		stop("! The parameter of simu_n is not a single valid value.");

	if( !missing(simu_p) && length(simu_p) > 1)
		stop("! The parameter of simu_p is not a single valid value.");
	
	if( !missing(simu_snp_rho) && length(simu_snp_rho) > 1)
		stop("! The parameter of simu_snp_rho is not a single valid value.");

	if( !missing(simu_snp_missing) && length(simu_snp_missing) > 1)
		stop("! The parameter of simu_snp_missing is not a single valid value.");

	if( !missing(simu_rho) && length(simu_rho) > 1)
		stop("! The parameter of simu_rho is not a single valid value.");
	
	if( !missing(simu_sigma2) && length(simu_sigma2) > 1)
		stop("! The parameter of simu_sigma2 is not a single valid value.");

	if( !missing(simu_mu) && length(simu_mu) > 1)
		stop("! The parameter of simu_mu is not a single valid value.");

	if ( length(simu_add_pos) != length(simu_add_effect ) )
		stop("! The length of simu_add_pos is same as simu_add_effect.");

	if ( length(simu_dom_pos) != length(simu_dom_effect ) )
		stop("! The length of simu_dom_pos is same as simu_dom_effect.");

	if ( length(simu_add_pos)>0 && length(which(simu_add_pos<=0 | simu_add_pos>simu_p))>0  )
		stop("! The parameter of simu_add_pos should be in correct SNP range.");

	if ( length(simu_dom_pos)>0 && length(which(simu_dom_pos<=0 | simu_dom_pos>simu_p))>0  )
		stop("! The parameter of simu_dom_pos should be in correct SNP range.");

	if ( length(simu_add_effect)>0 && length(which(is.na(simu_add_effect)))>0  )
		stop("! The parameter of simu_add_pos has NA values.");

	if ( length(simu_dom_effect)>0 && length(which(is.na(simu_dom_effect)))>0  )
		stop("! The parameter of simu_dom_pos has NA values.");

	if ( length(simu_cov_effect)>0 && length(which(is.na(simu_cov_effect)))>0  )
		stop("! The parameter of simu_cov_effect has NA values.");

	if ( length(simu_cov_range)>0 && length(which(is.na(simu_cov_range)))>0  )
		stop("! The parameter of simu_cov_range has NA values.");

	if ( length(simu_t_range)>0 && length(which(is.na(simu_t_range)))>0  )
		stop("! The parameter of simu_t_range has NA values.");
	
	if ( length(simu_t_range)!=2)
		stop("! The parameter of simu_t_range should be a valid range.");

	if ( length(simu_cov_range)!=2)
		stop("! The parameter of simu_cov_range should be a valid range.");

	if (simu_cov_range[1]>simu_cov_range[2])
		simu_cov_range<-c(simu_cov_range[2], simu_cov_range[1]);
		
	if (simu_t_range[1]>simu_t_range[2])
		simu_t_range<-c(simu_t_range[2], simu_t_range[1]);
	
	sigp<-unique(c(simu_add_pos, simu_dom_pos))
	simu_sigp <- length(sigp);
	simu_add_len <- length(simu_add_pos);
	simu_dom_len <- length(simu_dom_pos);
	simu_cov_count <- length(simu_cov_effect);

	err <- 0;
	out <- .C("bls_simulate", 
		   as.character(file.phe.out),		     # char* szPhe_out
		   as.character(file.snp.out), 		     # char* szSnp_out
		   as.integer(simu_grp), 		         # int nSimu_grp
		   as.integer(simu_n), 			         # int nSimu_n
		   as.integer(simu_p), 			         # int nSimu_p
		   as.double(simu_snp_rho), 	         # double fSimu_snp_rho
		   as.double(simu_snp_missing), 	     # double fSimu_snp_missing
		   as.double(simu_rho), 		         # double fSimu_rho
		   as.double(simu_sigma2), 		         # double fSimu_sigma2
		   as.double(simu_mu),			         # double fSimu_mu
		   as.integer(simu_cov_count),
		   as.double(as.vector(simu_cov_effect)),# double* pfSimu_cov_coeff
		   as.integer(simu_sigp),		         # int nSimu_sig_p
		   as.integer(simu_add_len),
		   as.integer(as.vector(simu_add_pos)),	 # int* nSimu_a_pos
		   as.double(as.vector(simu_add_effect)),# double* pfSimu_a_effect
		   as.integer(simu_dom_len),
		   as.integer(as.vector(simu_dom_pos)),	 # int* nSimu_d_pos
		   as.double(as.vector(simu_dom_effect)),# double* pfSimu_d_effect
		   as.double(as.vector(simu_cov_range)), # double* pfSimu_cov_range
		   as.double(as.vector(simu_t_range)),	 # double* pfSimu_t_range
		   as.integer(debug),
		   as.integer(err) );
	
	if( plink.format )
	{
		tb.snp <- read.csv(file.snp.out, header=T);
		r <- convert_simpe_to_plink( tb.snp, file.snp.out );
		
		unlink(file.snp.out);
		return(list(err=err,  
			file.simple.phe = file.phe.out,
			file.plink.bed = r$file.plink.bed,
			file.plink.bim = r$file.plink.bim,
			file.plink.fam = r$file.plink.fam));

	}
	
	return(list(err=err,  
		file.simple.phe = file.phe.out,
		file.simple.snp = file.snp.out));
}

bls.simple<-function(file.phe, file.snp, Y.name, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, fgwas.filter=FALSE, options=NULL)
{
	cat( "[ BLASSO SIMPLE ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.snp) || missing(Y.name) || missing(covar.names) )
		stop("! file.phe, file.snp, Y.name and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.name) && length(Y.name)==1 ) )
		stop("! The parameter of Y.name should be assigned with a outcome column in the phenotypic data.");
	if ( !missing("covar.names") && length(covar.names)>0 && !is.character(covar.names) )
		stop("! The parameter of covar.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(refit) && length(refit)==1 ) )
		stop("! The parameter of refit should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(add.used) && length(add.used)==1 ) )
		stop("! The parameter of add.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(dom.used) && length(dom.used)==1 ) )
		stop("! The parameter of dom.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
		stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Data File = ",  file.phe, "\n");
	cat("* Simpe SNP File = ",  file.snp, "\n");

	show_bls_parameters( Y.name, covar.names, refit, add.used, dom.used, fgwas.filter ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
		options0[names(options)] <- options;
		options <- options0;
	}
	
	cat( "Checking the optional items......\n");
	show_options( options);
	
	options$params <- list( file.phe     = file.phe, 
				file.snp     = file.snp, 
				Y.name       = Y.name, 
				covar.names  = covar.names, 
				refit        = refit, 
				add.used     = add.used, 
				dom.used     = dom.used, 
				fgwas.filter = fgwas.filter);

	
	r.bls <- list();
	r.filter <- list();
	
	if( options$nPiecewise.ratio==0 && !fgwas.filter )
	{
		cat( "Genetic Effect Analysis by BLASSO method......\n");
		r.bls <- .Call("bls_simple", 
			file.phe,
			file.snp, 
			Y.name, 
			paste(covar.names, collapse=","), 
			refit,
			add.used,
			dom.used,
			options$nMcmcIter,
			options$fBurnInRound,
			options$fRhoTuning,
			options$fQval.add,
			options$fQval.dom,
			ifelse( options$debug, 3, 1) );
	}
	else
	{
		simple <- read_simple_bls_data( file.phe, file.snp);

		subset_op <- function(snpmat, sub.idx)
		{
			return( snpmat[sub.idx,,drop=F] );
		}
		
		if(fgwas.filter)
		{
			r.filter <- snpmat_fgwas_filter( simple$phe.mat, simple$snp.mat, Y.name, NULL, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "BLS");
			
			if( r.filter$error ) stop(r.filter$err.info);
			if( is.null(r.filter$snp.mat)) return( wrap_fgwas_ret( r.filter, options) ); 
		
			r.bls <- snpmat_parallel(
				NROW( r.filter$snp.mat ),
				subset_op,
				r.filter$snp.mat,
				simple$phe.mat,
				Y.name, 
				NULL,
				covar.names,
				refit,
				add.used,
				dom.used,
				options$nPiecewise.ratio,
				options$nMcmcIter,
				options$fBurnInRound,
				options$fRhoTuning,
				options$fQval.add,
				options$fQval.dom,
				options$debug,
				options$nParallel.cpu,
				"BLS");
		}
		else
		{
			r.bls <- snpmat_parallel(
				NROW( simple$snp.mat ),
				subset_op,
				simple$snp.mat,
				simple$phe.mat,
				Y.name, 
				NULL,
				covar.names,
				refit,
				add.used,
				dom.used,
				options$nPiecewise.ratio,
				options$nMcmcIter,
				options$fBurnInRound,
				options$fRhoTuning,
				options$fQval.add,
				options$fQval.dom,
				options$debug,
				options$nParallel.cpu,
				"BLS");
		}
	}
	
	if(!is.null(r.bls) && !is.na(r.bls))
	{
		r <- wrap_BLS_ret( r.bls, r.filter, options );
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
}

bls.plink<-function( file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.name, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, fgwas.filter=FALSE, options=NULL, force.split=FALSE, plink.command=NULL )      
{
	cat( "[ BLASSO PLINK ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.plink.bed) || missing(file.plink.bim) || missing(file.plink.fam) || 
		missing(Y.name) || missing(covar.names) )
		stop("! file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.name and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.name) && length(Y.name)==1 ) )
		stop("! The parameter of Y.name should be assigned with a outcome column in the phenotypic data.");
	if ( !missing("covar.names") && length(covar.names)>0 && !is.character(covar.names) )
		stop("! The parameter of covar.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(refit) && length(refit)==1 ) )
		stop("! The parameter of refit should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(add.used) && length(add.used)==1 ) )
		stop("! The parameter of add.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(dom.used) && length(dom.used)==1 ) )
		stop("! The parameter of dom.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
		stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(force.split) && length(force.split)==1 ) )
		stop("! The parameter of force.split should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Data File = ",  file.phe, "\n");
	cat("* PLINK BED File = ",  file.plink.bed, "\n");
	cat("* PLINK BIM File = ",  file.plink.bim, "\n");
	cat("* PLINK FAM File = ",  file.plink.fam, "\n")
	cat("* PLINK Command = ",   plink.command, "\n")
	cat("* Force Split by PLINK Command = ", force.split, "\n")

	show_bls_parameters( Y.name, covar.names, refit, add.used, dom.used, fgwas.filter ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
		options0[names(options)] <- options;
		options <- options0;
	}
	
	cat( "Checking the optional items......\n");
	show_options( options);

	options$params <- list( file.phe     = file.phe, 
				file.plink.bed = file.plink.bed, 
				file.plink.bim = file.plink.bim, 
				file.plink.fam = file.plink.fam,
				Y.name       = Y.name, 
				covar.names  = covar.names, 
				refit        = refit, 
				add.used     = add.used, 
				dom.used     = dom.used, 
				fgwas.filter = fgwas.filter);
	
	pd <- list();
	r.filter <- list();
	
	if( force.split || !try_load_plink( file.plink.bed,  file.plink.bim, file.plink.fam ) )
	{
		# It is bigdata which need to split it into chromosome unit
		# The following will split the data and force to do fGWAS filter.

		r.filter <- plink_fgwas_bigdata ( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe, plink.command, 
										Y.name, NULL, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "BLS");
		if( r.filter$error ) stop(r.filter$err.info);

		fgwas.filter <- TRUE;

		pd <- list(phe.mat=r.filter$phe.mat, snp.mat=r.filter$snp.mat);
	}	
	else
	{
		pd <- load_plink_binary( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe );
		if( is.null(pd) )
			stop("! Failed to load PLINK dataset.");

		if(fgwas.filter)
		{
			# call FGWAS.R to do FILTER and the bls_snpmat
			r.filter <- plink_fgwas_filter( pd, Y.name, NULL, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "BLS")

			if( r.filter$error ) stop(r.filter$err.info);
		}				
	}
	
	r.bls <- list();

	if( fgwas.filter )
	{
		if( is.null(r.filter$snp.mat)) return( wrap_fgwas_ret( r.filter, options) ); 
	
		subset_op <- function(snpmat, sub.idx)
		{
			return( snpmat[sub.idx,,drop=F] );
		}

		r.bls <- snpmat_parallel(
			NROW( r.filter$snp.mat ),
			subset_op,
			r.filter$snp.mat,
			pd$phe.mat,
			Y.name, 
			NULL,
			covar.names,
			refit,
			add.used,
			dom.used,
			options$nPiecewise.ratio,
			options$nMcmcIter,
			options$fBurnInRound,
			options$fRhoTuning,
			options$fQval.add,
			options$fQval.dom,
			options$debug,
			options$nParallel.cpu,
			"BLS");

	}
	else
	{
		subset_op <- function(snpmat, sub.idx)
		{
			snp.sub <- get_plink_subsnp(snpmat, sub.idx );
			snp.mat <- cbind( snp.sub$info[,c(2,3)], snp.sub$snp )
			return( snp.mat );
		}
		
		r.bls <- snpmat_parallel(
			NCOL( pd$snp.mat$genotypes ),
			subset_op,
			pd$snp.mat,
			pd$phe.mat,
			Y.name, 
			NULL,
			covar.names,
			refit,
			add.used,
			dom.used,
			options$nPiecewise.ratio,
			options$nMcmcIter,
			options$fBurnInRound,
			options$fRhoTuning,
			options$fQval.add,
			options$fQval.dom,
			options$debug,
			options$nParallel.cpu,
			"BLS");
	}
	
	if(!is.null(r.bls) && !is.na(r.bls))
	{
		r <- wrap_BLS_ret( r.bls, r.filter, options );
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
}

bls.plink.tped<-function( file.phe, file.plink.tped, file.plink.tfam, Y.name, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, options=NULL)
{
	cat( "[ BLASSO PLINK.tped ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.plink.tped) || missing(file.plink.tfam) ||
		missing(Y.name) || missing(covar.names) )
		stop("! file.phe, file.plink.tped, file.plink.tfam, Y.name and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.name) && length(Y.name)==1 ) )
		stop("! The parameter of Y.name should be assigned with a outcome column in the phenotypic data.");
	if ( !missing("covar.names") && length(covar.names)>0 && !is.character(covar.names) )
		stop("! The parameter of covar.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(refit) && length(refit)==1 ) )
		stop("! The parameter of refit should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(add.used) && length(add.used)==1 ) )
		stop("! The parameter of add.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(dom.used) && length(dom.used)==1 ) )
		stop("! The parameter of dom.used should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Data File = ",  file.phe, "\n");
	cat("* PLINK TPED File = ",  file.plink.tped, "\n");
	cat("* PLINK TFAM File = ",  file.plink.tfam, "\n");

	show_bls_parameters( Y.name, covar.names, refit, add.used, dom.used, fgwas.filter=FALSE ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
		options0[names(options)] <- options;
		options <- options0;
	}
	
	cat( "Checking the optional items......\n");
	show_options( options);
	
	cat( "Genetic Effect Analysis by BLASSO method......\n");
	r.filter <- list();
	r.bls <- .Call("bls_plink_tped", 
			file.phe,
			file.plink.tped, 
			file.plink.tfam, 
			Y.name, 
			paste(covar.names, collapse=","), 
			refit,
			add.used,
			dom.used,
			options$nMcmcIter,
			options$fBurnInRound,
			options$fRhoTuning,
			options$fQval.add,
			options$fQval.dom,
			ifelse( options$debug, 3, 1) );


	options$params <- list( file.phe     = file.phe, 
				file.plink.tped = file.plink.tped, 
				file.plink.tfam = file.plink.tfam, 
				Y.name       = Y.name, 
				covar.names  = covar.names, 
				refit        = refit, 
				add.used     = add.used, 
				dom.used     = dom.used, 
				fgwas.filter = FALSE);

	if(!is.null(r.bls) && !is.na(r.bls))
	{
		r <- wrap_BLS_ret( r.bls, r.filter, options );
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
}

bls.snpmat<-function(phe.mat, snp.mat, Y.name, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, fgwas.filter=FALSE, options=NULL)
{
	cat( "[ BLASSO SNPMAT ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(phe.mat) || missing(snp.mat) || missing(Y.name) || missing(covar.names) )
		stop("! phe.mat, snp.mat, Y.name and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.name) && length(Y.name)==1 ) )
		stop("! The parameter of Y.name should be assigned with a outcome column in the phenotypic data.");
	if ( !missing("covar.names") && length(covar.names)>0 && !is.character(covar.names) )
		stop("! The parameter of covar.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(refit) && length(refit)==1 ) )
		stop("! The parameter of refit should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(add.used) && length(add.used)==1 ) )
		stop("! The parameter of add.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(dom.used) && length(dom.used)==1 ) )
		stop("! The parameter of dom.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
		stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Matrix = ",  dim(phe.mat), "\n");
	cat("* SNP Matrix = ",  dim(snp.mat), "\n");

	show_bls_parameters( Y.name, covar.names, refit, add.used, dom.used, fgwas.filter ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
		options0[names(options)] <- options;
		options <- options0;
	}
	
	cat( "Checking the optional items......\n");
	show_options( options);

	options$params <- list( Y.name       = Y.name, 
				covar.names  = covar.names, 
				refit        = refit, 
				add.used     = add.used, 
				dom.used     = dom.used, 
				fgwas.filter = fgwas.filter);

	if( class(phe.mat)=="data.frame" )
	{
		cat("Phenotypic data frame is converted to the matrix class.\n");  
		phe.colnames <- colnames(phe.mat); 
		phe.rownames <- rownames(phe.mat); 
		phe.mat <- matrix(as.numeric(as.matrix(phe.mat, rownames.force=NA)), ncol=NCOL(phe.mat))
		colnames(phe.mat) <- phe.colnames;
	 	rownames(phe.mat) <- phe.rownames;
	}
	
	r.bls <- list();
	r.filter <- list();

	if( options$nPiecewise.ratio==0 && !fgwas.filter )
	{
		cat( "Genetic Effect Analysis by BLASSO method......\n");

		r.bls <- .Call("bls_snpmat", 
			as.matrix( phe.mat ),
			as.matrix( snp.mat*1.0 ),
			Y.name, 
			paste(covar.names, collapse=","), 
			refit,
			add.used,
			dom.used,
			options$nMcmcIter,
			options$fBurnInRound,
			options$fRhoTuning,
			options$fQval.add,
			options$fQval.dom,
			ifelse( options$debug, 3, 1) );
	}
	else
	{
		subset_op <- function(snpmat, sub.idx)
		{
			return( snpmat[sub.idx,,drop=F] );
		}

		if(fgwas.filter)
		{
			r.filter <- snpmat_fgwas_filter( phe.mat, snp.mat, Y.name, NULL, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "BLS")

			if( r.filter$error ) stop(r.filter$err.info);
			if( is.null(r.filter$snp.mat)) return( wrap_fgwas_ret( r.filter, options) ); 
		
			r.bls <- snpmat_parallel(
				NROW(r.filter$snp.mat),
				subset_op,
				r.filter$snp.mat,
				phe.mat,
				Y.name, 
				NULL,
				covar.names,
				refit,
				add.used,
				dom.used,
				options$nPiecewise.ratio,
				options$nMcmcIter,
				options$fBurnInRound,
				options$fRhoTuning,
				options$fQval.add,
				options$fQval.dom,
				options$debug,
				options$nParallel.cpu,
				"BLS");
		}
		else
		{
			r.bls <- snpmat_parallel(
				NROW(snp.mat),
				subset_op,
				snp.mat,
				phe.mat,
				Y.name, 
				NULL,
				covar.names,
				refit,
				add.used,
				dom.used,
				options$nPiecewise.ratio,
				options$nMcmcIter,
				options$fBurnInRound,
				options$fRhoTuning,
				options$fQval.add,
				options$fQval.dom,
				options$debug,
				options$nParallel.cpu,
				"BLS");
		}
	}
	
	if(!is.null(r.bls) && !is.na(r.bls))
	{
		r <- wrap_BLS_ret( r.bls, r.filter, options );
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
}

summary.BLS.ret<-function(object, ...)
{
	r.bls <- object;
	r.sum.ret <- list();

	if(!is.null( r.bls$refit_cov ) && NROW( r.bls$refit_cov )>0 )
	{
		re1 <- r.bls$refit_cov;
		r.sum.ret$refit_cov <- data.frame(  Sig=ifelse(re1[,1], "Yes", "---"), 
					Median  = round(re1[,2], digits=3), 
					CI.025  = round(re1[,3], digits=3), 
					CI.975  = round(re1[,4], digits=3) );
		rownames(r.sum.ret$refit_cov) <- rownames( re1 );
	}

	if(!is.null( r.bls$refit ) && NROW( r.bls$refit )>0 )
	{
		re2 <- r.bls$refit;
		r.sum.ret$refit <- data.frame( Chr = re2[,1], 
					Pos        = re2[,2], 
					Add.Sig    = ifelse(re2[,3], "Yes", "---") ,
					Add.Median = round(re2[,4], digits=3), 
					Add.CI.025 = round(re2[,5], digits=3), 
					Add.CI.975 = round(re2[,6], digits=3),
					Dom.Sig    = ifelse(re2[,7], "Yes", "---") ,
					Dom.Median = round(re2[,8], digits=3), 
					Dom.CI.025 = round(re2[,9], digits=3), 
					Dom.CI.975 = round(re2[,10], digits=3),
					H2         = round(re2[,11], digits=3) );
		rownames(r.sum.ret$refit)<- rownames(re2);
	}
	
	if(!is.null( r.bls$varsel_cov ) && NROW( r.bls$varsel_cov )>0 )
	{
		re3 <- r.bls$varsel_cov;
		r.sum.ret$varsel_cov <- data.frame(  Sig=ifelse(re3[,1], "Yes", "---"), 
					Median  = round(re3[,2], digits=3), 
					CI.025  = round(re3[,3], digits=3), 
					CI.975  = round(re3[,4], digits=3) );
		rownames(r.sum.ret$varsel_cov) <- rownames( re3 );
	}
	
	if(!is.null(r.bls$varsel))
	{
		re4 <- r.bls$varsel;
		var.sig <- which( re4[,3]>0 | re4[,7]>0);
		if(length(var.sig)>0)
		{
			re4 <- re4[var.sig,,drop=F];
			
			r.sum.ret$varsel <- data.frame( Chr= re4[ , 1], 
						Pos        = re4[ , 2], 
						Add.Sig    = ifelse(re4[ ,3 ], "Yes", "---") ,
						Add.Median = round( re4[ ,4 ], digits=3 ), 
						Add.CI.025 = round( re4[ ,5 ], digits=3 ), 
						Add.CI.975 = round( re4[ ,6 ], digits=3 ),
						Dom.Sig    = ifelse(re4[ ,7 ], "Yes", "---") ,
						Dom.Median = round( re4[ ,8 ], digits=3 ), 
						Dom.CI.025 = round( re4[ ,9 ], digits=3 ), 
						Dom.CI.975 = round( re4[ ,10], digits=3 ),
						H2         = round( re4[ ,11], digits=3 ) );
			rownames(r.sum.ret$varsel)<- rownames(re4);
		}
	}

	if(!is.null(r.bls$fgwas))
	{
		re7 <- r.bls$fgwas;
		fgwas.sig <- which( re7[,7] <= r.bls$options$fgwas.cutoff );
		if(length(fgwas.sig)>0)
		{
			fgwas_sigs <- re7[ fgwas.sig, , drop=F];
			fgwas.sig.inc <- order(fgwas_sigs[,7]);
			r.sum.ret$fgwas_sig <- fgwas_sigs[fgwas.sig.inc,];
		}
		
		if(!is.null(r.sum.ret$varsel))
			r.sum.ret$varsel <- cbind(r.sum.ret$varsel, fgwas.pvalue=find_fgwas_pvalue( r.bls$fgwas, rownames(r.sum.ret$varsel) ) ) ;

		if(!is.null(r.sum.ret$refit))
			r.sum.ret$refit <- cbind(r.sum.ret$refit, fgwas.pvalue=find_fgwas_pvalue( r.bls$fgwas, rownames(r.sum.ret$refit) ) ) ;
	}

	class(r.sum.ret) <- "sum.BLS.ret";
	
	r.sum.ret
}

print.sum.BLS.ret<-function(x, ...)
{
	r.sum.ret <- x;
	if(!is.null(r.sum.ret$fgwas_sig))
	{
		cat("--- Significant SNPs Estimate by fGWAS method:", NROW(r.sum.ret$fgwas_sig), "SNPs\n");
		if( NROW(r.sum.ret$fgwas_sig)>25 )
		{
			cat("Top 25 SNPs:\n");
			show(r.sum.ret$fgwas_sig[1:25,,drop=F]);
		}
		else	
			show(r.sum.ret$fgwas_sig);
	}

	if(!is.null(r.sum.ret$varsel_cov))
	{
		cat("--- Covariate Estimate in Varsel Procedure:\n");
		show(r.sum.ret$varsel_cov);
	}
	
	if(!is.null(r.sum.ret$varsel))
	{
		cat("--- Variable Selection Result:", NROW(r.sum.ret$varsel), "SNPs\n" );
		if( NROW(r.sum.ret$varsel)>25 )
		{
			cat("Top 25 SNPs:\n");
			show(r.sum.ret$varsel[1:25,,drop=F]);
		}
		else
			show(r.sum.ret$varsel );
	}

	if(!is.null(r.sum.ret$refit_cov))
	{
		cat("--- Covariate Estimate in Refit Procedure:\n");
		show(r.sum.ret$refit_cov);
	}
	
	if(!is.null(r.sum.ret$refit))
	{
		cat("--- Refit Result:", NROW(r.sum.ret$refit), "SNPs\n" );
		show(r.sum.ret$refit);
	}
}
	
#--------------------------------------------------------------
# plot_adh2
# 
# Input:pvs[,1]  snp_name
#       pvs[,2]  chromoseom no
#       pvs[,3]  position
#       pvs[,4]  Additive
#       pvs[,5]  Dominant
#       pvs[,6]  H2
# fig.file: pdf file name
#
# Used by: BLS
#--------------------------------------------------------------
plot.BLS.ret<-function( x, y=NULL, ..., fig.prefix=NULL )
{
	r.bls <- x;
	
	if( missing(fig.prefix)) fig.prefix <- "bls.plot";
	
	if(!is.null(r.bls$fgwas))
	{
		filter.man <- r.bls$fgwas[, c(1,2,5), drop=F]
		draw_man_fgwas( filter.man, fig.prefix, "fgwas" );
	}
	else
		cat("! No fGWAS filter results.\n");		
		
	if(!is.null(r.bls$varsel))
	{
		varsel.man <- r.bls$varsel[, c(1,2,4,8,11), drop=F];
		draw_man_adh2( varsel.man, fig.prefix, "varsel" );
	}
	else
		cat("! No varible selection results.\n");		

	if(!is.null(r.bls$refit))
	{
		refit.man <- r.bls$refit[, c(1,2,4,8,11), drop=F];
		draw_man_adh2( refit.man, fig.prefix, "refit" );
	}
	else
		cat("! No refit results.\n");		
}

wrap_BLS_ret<-function(r.bls, r.filter, options)
{
	cat( "Wrapping the results ......\n");

	if(!is.null(r.bls) && !is.na(r.bls))
	{
		if (!is.null(r.bls$varsel))
			colnames(r.bls$varsel)   <- c("grp", "pos", "add.sig", "add.mu", "add.min", "add.max", "dom.sig", "dom.mu", "dom.min", "dom.max", "h2");
		if (!is.null(r.bls$refit))
			colnames(r.bls$refit)    <- c("grp", "pos", "add.sig", "add.mu", "add.min", "add.max", "dom.sig", "dom.mu", "dom.min", "dom.max", "h2");
		if (!is.null(r.bls$varsel_cov))
		{
			colnames(r.bls$varsel_cov) <- c("cov.sig", "cov.mu", "cov.min", "cov.max");
			rownames(r.bls$varsel_cov) <- c("intercept", options$params$covar.names);
		}
		
		if (!is.null(r.bls$refit_cov))
		{
			colnames(r.bls$refit_cov) <- c("cov.sig", "cov.mu", "cov.min", "cov.max");
			rownames(r.bls$refit_cov) <- c("intercept",  options$params$covar.names);
		}
		
		if(!is.null(r.filter)) r.bls$fgwas <- r.filter$r;
		
		r.bls$options <- options;

		class(r.bls) <- "BLS.ret";
	}
		
	return(r.bls);
}

get_default_options<-function()
{
	options=list(
				nParallel.cpu = 0,    
				nPiecewise.ratio = 2,    
				nMcmcIter = 2000,    
				fBurnInRound = 0.3,   
				fRhoTuning = 0.095,  
				fQval.add  = 0.05,   
				fQval.dom  = 0.09,  
				fgwas.cutoff = 0.05,
				debug      = F ) 

	return(options);	
}

show_options<-function(options)
{
	cat( "* Parallel Computing: ", ifelse( options$nParallel.cpu>1, "Yes,", "No,"), options$nParallel.cpu,  "CPU(s)\n");
	cat( "* Piecewise Ratio: ",  options$nPiecewise.ratio,  "\n");

	cat( "* Threshold of fGWAS filter: ",  options$fgwas.cutoff, "\n");
	cat( "* Iteration of  Markov chain: ",  options$nMcmcIter, "\n");
	cat( "* fBurnInRound: ",  options$fBurnInRound, "\n");
	cat( "* fRhoTuning: ",  options$fRhoTuning, "\n");
	cat( "* fQval.add: ",  options$fQval.add, "\n");
	cat( "* fQval.dom: ",  options$fQval.dom , "\n");
	
	cat( "* Debug Output: ", ifelse( options$debug, "Yes", "No"),"\n");
}

show_bls_parameters<-function( Y.name, covar.names, refit, add.used, dom.used, fgwas.filter ) 
{
	cat( "* Response Variable =",   Y.name, "\n");
	cat( "* Covariate Columns =",  covar.names, "\n");
	cat( "* fGWAS Filter Used =",  ifelse( fgwas.filter, "Yes", "No"), "\n");
	cat( "* Additive Effects Used =",  ifelse( add.used, "Yes", "No"), "\n");
	cat( "* Dominant Effects Used =",  ifelse( dom.used, "Yes", "No"), "\n");
	cat( "* Refit Procedure =",   ifelse( refit, "Yes", "No"), "\n");
}


get_sig_bls_snp <- function( r.bls )
{
	if(is.null(r.bls$varsel)) return( NULL );

	idx.sig <- which( r.bls$varsel[,3]!=0 | r.bls$varsel[,7]!=0 )
	if (length(idx.sig)==0) return(NULL);
	
	return( idx.sig );
}

read_simple_bls_data <- function( file.phe, file.snp, bImputed=T )
{
	tb.phe <- read.csv(file.phe, header=T);
	rownames(tb.phe) <- tb.phe[,1];
	tb.phe <- tb.phe[,-1, drop=F];
	
	tb.snp <- read.csv(file.snp, header=T);

	cat("Checking data files......\n");
	cat("* Individuals:", NCOL(tb.snp)-2, "\n");
	cat("* SNPs:", NROW(tb.snp), "\n");
	
	if(bImputed) tb.snp <- impute_simple_snp(tb.snp);

	return(list(phe.mat=tb.phe, snp.mat=tb.snp));
}

bls.best.qval<-function( r.bls, snp.names )
{
	if( class(r.bls) != "BLS.ret" )	
	{
		cat("! r.bls is NOT a BLS.ret object.\n");
		return(NULL);
	}

	find_Qtable<-function( Q.table, snp.names)
	{
		pv <- rep(NA, length(snp.names));
	
		idx <- match( snp.names, rownames(Q.table));
		if( length(which(!is.na(idx)) ) > 0)
			pv[ !is.na(idx) ] <-  Q.table[ idx[!is.na(idx)], 1 ];
	
		return(pv);
	}
	
	vs.Qval.add <- rep(NA, length(snp.names));
	vs.Qval.dom <- rep(NA, length(snp.names));
	
	if( is.null(r.bls$varsel_Qbest ) )
	{
		cat("! No Qbest matrix for variable selection.\n");
	}
	else
	{
		if(any(r.bls$varsel_Qbest[,1]!=0))
			vs.Qval.add <- find_Qtable( r.bls$varsel_Qbest[,1,drop=F], snp.names );
		if(any(r.bls$varsel_Qbest[,5]!=0))
			vs.Qval.dom <- find_Qtable( r.bls$varsel_Qbest[,5,drop=F], snp.names );
		
	}
	
	refit.Qval.add <- rep(NA, length(snp.names));
	refit.Qval.dom <- rep(NA, length(snp.names));
	if( is.null(r.bls$refit_Qbest ) )
	{
		cat("! No Qbest matrix for refit procedure.\n");
	}
	else
	{
		if(any(r.bls$refit_Qbest[,1]!=0))
			refit.Qval.add <- find_Qtable( r.bls$refit_Qbest[,1,drop=F], snp.names );
		if(any(r.bls$refit_Qbest[,5]!=0))
			refit.Qval.dom <- find_Qtable( r.bls$refit_Qbest[,5,drop=F], snp.names );
	
	}

	fgwas.pv <- find_fgwas_pvalue( r.bls$fgwas, snp.names)

	r.qval <- cbind( fgwas.pv, vs.Qval.add, vs.Qval.dom, refit.Qval.add, refit.Qval.dom);
	colnames(r.qval) <- c("fgwas.pv", "Qval.vs.add", "Qval.vs.dom", "Qval.refit.add", "Qval.refit.dom");
	rownames(r.qval) <- snp.names;

	return(r.qval);
}

bls.qval.cutoff<-function( r.bls, qval.add, qval.dom, refit.select = FALSE )
{
	if( class(r.bls) != "BLS.ret" )	
	{
		cat("! r.bls is NOT a BLS.ret object.\n");
		return(NULL);
	}

	if( is.null(r.bls$varsel_Qbest ) && !refit.select )
	{
		cat("! No Qbest matrix for variable selection.\n");
	}
	else
	{
		idx.vs.add <- c();
		idx.vs.dom <- c();

		if(any(r.bls$varsel_Qbest[,1]!=0))
			idx.vs.add <- which( r.bls$varsel_Qbest[, 1] <= qval.add); 
		if(any(r.bls$varsel_Qbest[,5]!=0))
			idx.vs.dom <- which( r.bls$varsel_Qbest[, 5] <= qval.dom); 
		
		idx.vs <- sort( unique(c(idx.vs.add, idx.vs.dom)) )
		if(length(idx.vs)>0)
			return( rownames(r.bls$varsel_Qbest)[idx.vs])
		else
			return(NULL);			
	}
	
	if( is.null(r.bls$refit_Qbest )  && refit.select )
	{
		cat("! No Qbest matrix for refit procedure.\n");
	}
	else
	{
		idx.refit.add <- c();
		idx.refit.dom <- c();

		if(any(r.bls$refit_Qbest[,1]!=0))
			idx.refit.add <- which( r.bls$refit_Qbest[, 1] <= qval.add); 
		if(any(r.bls$refit_Qbest[,5]!=0))
			idx.refit.dom <- which( r.bls$refit_Qbest[, 5] <= qval.dom); 
		
		idx.vs <- sort( unique(c(idx.refit.add, idx.refit.dom)) )
		if(length(idx.vs)>0)
			return( rownames(r.bls$refit_Qbest)[idx.vs])
		else
			return(NULL);			
	}

	return(NULL);
}
