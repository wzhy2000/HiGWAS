#config param
#
# debug=1
# max_iter=2000
# rho_tuning=0.09
# burn_in_round=0.2
# fQVal=0.05
#simu_a_effect[1] = 1, 1.04, 0.885, -2.055, 0.545
#simu_a_effect[2] = 2, 1.17, -0.20, 0.74, -4.715
#simu_a_effect[3] = 3, 1.40, -2.25, 1.00,  0.00

#simu_d_effect[1] = 3, 1.49, -2.135, 4.82, 1.425
#simu_d_effect[2] = 4, 1.045, 1.320, 1.905,  1.535
#simu_d_effect[3] = 5, 1.265, -1.225, 2.710, -1.96

gls.simulate<-function( file.phe.out, file.snp.out, simu_grp=1, simu_n=500, simu_p=1000, 
			simu_snp_rho   = 0.1, 
			simu_snp_missing= 0.002, 
			simu_rho       = 0.4, 
			simu_sigma2    = 16,
			simu_mu        = c( 13.395, -3.08, 1.875, -3.195 ),
			simu_cov_range = c( -1, 1 ),
			simu_cov_effect= array(c(0,0,0,0), dim=c(1,4)), 
			simu_add_pos   = c( 1,2,3 ), 
			simu_add_effect= array(c( 1.04, 0.885, -2.055, 0.545, 
									1.17, -0.20, 0.74, -4.715,
									1.40, -2.25, 1.00,  0.00), dim=c(3,4)),
			simu_dom_pos   = c( 3,4,5 ), 
			simu_dom_effect= array(c( 1.49, -2.135, 4.82, 1.425, 
									1.045, 1.320, 1.905,  1.535,
									1.265, -1.225, 2.710, -1.96), dim=c(3,4)),
			simu_z_range   = c( 20, 80 ), 
			simu_z_count   = c(  5, 12 ), 
			plink.format	=FALSE,
			debug          = FALSE )
{
	if( !missing(simu_grp) && length(simu_grp) > 1)
		stop("The parameter of simu_grp is not a single valid value.");

	if( !missing(simu_n) && length(simu_n) > 1)
		stop("The parameter of simu_n is not a single valid value.");

	if( !missing(simu_p) && length(simu_p) > 1)
		stop("The parameter of simu_p is not a single valid value.");
	
	if( !missing(simu_snp_rho) && length(simu_snp_rho) > 1)
		stop("The parameter of simu_snp_rho is not a single valid value.");

	if( !missing(simu_snp_missing) && length(simu_snp_missing) > 1)
		stop("The parameter of simu_snp_missing is not a single valid value.");

	if( !missing(simu_rho) && length(simu_rho) > 1)
		stop("The parameter of simu_rho is not a single valid value.");
	
	if( !is.na(simu_sigma2) && length(simu_sigma2) > 1)
		stop("The parameter of simu_sigma2 is not a single valid value.");

	if( length(which(is.na(simu_mu)))>0 || length(simu_mu) != 4 )
		stop("The parameter of simu_mu is not a vector with 4 numeric valid values.");

	if( length(which(is.na(simu_cov_range)))>0 || length(simu_cov_range) != 2 )
		stop("The parameter of simu_cov_range is not a enclosed range.");

	if( length(which(is.na(simu_z_range)))>0 || length(simu_z_range) != 2 )
		stop("The parameter of simu_z_range is not a enclosed range.");

	if( length(which(is.na(simu_z_count)))>0 || length(simu_z_count) != 2 )
		stop("The parameter of simu_z_count is not a enclosed range.");

	if(!is.matrix( simu_cov_effect ) )
		stop("The parameter of simu_cov_effect is not a matrix.");
	if(NCOL(simu_cov_effect)!=4)
		stop("The parameter of simu_cov_effect is not a matrix with 4 columns.");
	if(length(which(is.na(simu_cov_effect)))>0 )
		stop("The parameter of simu_cov_effect has NA values.");
	
	if(length(simu_add_pos)>0 )
	{
		if(!is.matrix(simu_add_effect))
			stop("The parameter of simu_add_effect is not a matrix.");
		if(NCOL(simu_add_effect)!=4)
			stop("The parameter of simu_add_effect is not a matrix with 4 columns.");
		if(length(which(is.na(simu_add_effect)))>0 )
			stop("The parameter of simu_add_effect has NA values.");
		if ( length(simu_add_pos) != NROW(simu_add_effect ) )
			stop("! The length of simu_add_pos is same as simu_add_effect.");
	}
	
	if(length(simu_dom_pos)>0 )
	{
		if(!is.matrix(simu_dom_effect))
			stop("The parameter of simu_dom_effect is not a matrix.");
		if(NCOL(simu_dom_effect)!=4)
			stop("The parameter of simu_dom_effect is not a matrix with 4 columns.");
		if(length(which(is.na(simu_dom_effect)))>0 )
			stop("The parameter of simu_dom_effect has NA values.");
		if ( length(simu_dom_pos) != NROW(simu_dom_effect ) )
			stop("! The length of simu_dom_pos is same as simu_dom_effect.");
	}

	if ( length(simu_add_pos)>0 && length(which(simu_add_pos<=0 | simu_add_pos>simu_p))>0  )
		stop("! The parameter of simu_add_pos should be in correct SNP range.");

	if ( length(simu_dom_pos)>0 && length(which(simu_dom_pos<=0 | simu_dom_pos>simu_p))>0  )
		stop("! The parameter of simu_dom_pos should be in correct SNP range.");

	simu_sig_add <- NROW(simu_add_effect);
	simu_sig_dom <- NROW(simu_dom_effect);
	
	sigp  <-unique(c(simu_add_pos, simu_dom_pos));
	simu_sigp <- length(sigp);
	err   <- 0;
	
	simu_add_mat <- NULL;
	if ( length(simu_add_pos)>0 ) simu_add_mat <- as.matrix( cbind( simu_add_pos, simu_add_effect) );

	simu_dom_mat <- NULL;
	if ( length(simu_dom_pos)>0 ) simu_dom_mat <- as.matrix( cbind( simu_dom_pos, simu_dom_effect) );
	
	out <- .C("gls_simulate", 
			as.character(file.phe.out),				# char* szPhe_out
			as.character(file.snp.out),				# char* szSnp_out
			as.integer(simu_grp),					# int nSimu_grp
			as.integer(simu_n), 					# int nSimu_n
			as.integer(simu_p), 					# int nSimu_p
			as.double(simu_snp_rho),				# double fSimu_snp_rho
			as.double(simu_snp_missing),			# double fSimu_snp_missing
			as.double(simu_rho),					# double fSimu_rho
			as.double(simu_sigma2),					# double fSimu_sigma2
			as.double(as.vector(simu_mu)),			# double* pfSimu_mu
			as.integer(NROW(simu_cov_effect)),		# int nSimu_covar_len
			as.double(as.vector(simu_cov_range)),	# double* pfSimu_covar_range
			as.double(as.matrix(simu_cov_effect)),	# double* pfSimu_covar_effect
			as.integer(simu_sigp),		 	      	# int nSimu_sig_p
			as.integer(simu_sig_add),		  	  	# int nSimu_add_len
			as.double( simu_add_mat ),             	# double* pfSimu_add_effect
			as.integer(simu_sig_dom),			  	# int nSimu_dom_len
			as.double( simu_dom_mat ),             	# double* pfSimu_dom_effect
			as.double(as.vector(simu_z_range)),	  	# double* simu_z_range
			as.integer(as.vector(simu_z_count)),   	# int* pnSimu_z_count
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

gls.simple<-function(file.phe, file.snp, Y.prefix, Z.prefix, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, fgwas.filter=FALSE, options=NULL )  
{
	cat( "[ GLASSO SIMPLE ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.snp) || missing(Y.prefix) || missing(Z.prefix) || missing(covar.names) )
		stop("! file.phe, file.snp, Y.prefix, Z.prefix and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
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

	show_gls_parameters( Y.prefix, Z.prefix, covar.names, refit, add.used, dom.used, fgwas.filter ) ;

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
	
	options$params <- list( 
				file.phe       = file.phe, 
				file.snp       = file.snp, 
				Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covar.names    = covar.names, 
				refit          = refit, 
				add.used       = add.used, 
				dom.used       = dom.used, 
				fgwas.filter   = fgwas.filter);	
	
	r.gls <- list();
	r.filter <- list();

	if( options$nPiecewise.ratio==0 && !fgwas.filter )
	{
		cat( "Genetic Effect Analysis by GLASSO method......\n");

		r.gls <- .Call("gls_simple", 
				file.phe,
				file.snp, 
				Y.prefix, 
				Z.prefix, 
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
		simple <- read_simple_gls_data( file.phe, file.snp );	

		subset_op <- function(snpmat, sub.idx)
		{
			return( snpmat[sub.idx,,drop=F] );
		}
		
		if(fgwas.filter)
		{
			r.filter <- snpmat_fgwas_filter( simple$phe.mat, simple$snp.mat, Y.prefix, Z.prefix, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");

			if( r.filter$error ) stop(r.filter$err.info);
			if( is.null(r.filter$snp.mat) ) return ( wrap_fgwas_ret( r.filter, options) );
		
			r.gls <- snpmat_parallel(
				NROW( r.filter$snp.mat ),
				subset_op,
				r.filter$snp.mat,
				simple$phe.mat,
				Y.prefix, 
				Z.prefix,
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
				"GLS");
		}
		else
		{
			r.gls <- snpmat_parallel(
				NROW( simple$snp.mat ),
				subset_op,
				simple$snp.mat,
				simple$phe.mat,
				Y.prefix, 
				Z.prefix,
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
				"GLS");
		}
	}
	

	if(!is.null(r.gls) && !is.na(r.gls))
	{
		r <- wrap_GLS_ret(r.gls, r.filter, options);
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
}

gls.plink<-function( file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, fgwas.filter=FALSE, options=NULL, force.split=FALSE, plink.command=NULL )        
{
	cat( "[ GLASSO PLINK ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.plink.bed) || missing(file.plink.bim) || missing(file.plink.fam) || 
		missing(Y.prefix) || missing(Z.prefix) || missing(covar.names) )
		stop("! file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
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
	
	show_gls_parameters( Y.prefix, Z.prefix, covar.names, refit, add.used, dom.used, fgwas.filter ) ;

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

	options$params <- list( file.phe       = file.phe, 
				file.plink.bed = file.plink.bed, 
				file.plink.bim = file.plink.bim, 
				file.plink.fam = file.plink.fam,				
				Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covar.names    = covar.names, 
				refit          = refit, 
				add.used       = add.used, 
				dom.used       = dom.used, 
				fgwas.filter   = fgwas.filter);
	
	pd <- list();
	r.filter <- list();
	
	if( force.split || !try_load_plink( file.plink.bed,  file.plink.bim, file.plink.fam ) )
	{
		# It is bigdata which need to split it into chromosome unit
		# The following will split the data and force to do fGWAS filter.

		if( !is.null(options$fgwas.rdata)  && load(options$fgwas.rdata)=="r.filter" )
		{
		}
		else
		{
		r.filter <- plink_fgwas_bigdata ( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe, plink.command, 
										Y.prefix, Z.prefix, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");

		if( r.filter$error ) stop(r.filter$err.info);
		if( !is.null( options$fgwas.rdata ) ) 
			save(r.filter, file=options$fgwas.rdata);
		}
		
		fgwas.filter <- TRUE;
		pd <- list(phe.mat=r.filter$phe.mat, snp.mat=r.filter$snp.mat);
	}	
	else
	{
		pd <- load_plink_binary( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe );
		if( is.null(pd) )
			stop("Failed to load PLINK dataset!");

		if(fgwas.filter)
		{
			# call FGWAS.R to do FILTER and the gls__snpmat
			r.filter <- plink_fgwas_filter( pd, Y.prefix, Z.prefix, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");

			if( r.filter$error ) stop(r.filter$err.info);
			if( !is.null( options$fgwas.rdata ) ) 
				save(r.filter, file=options$fgwas.rdata);
		}				
	}
	
	r.gls <- list();

	if( fgwas.filter)
	{
		if( is.null(r.filter$snp.mat) ) return ( wrap_fgwas_ret( r.filter, options) );

		subset_op <- function(snpmat, sub.idx)
		{
			return( snpmat[sub.idx,,drop=F] );
		}

		r.gls <- snpmat_parallel(
			NROW( r.filter$snp.mat ),
			subset_op,
			r.filter$snp.mat,
			pd$phe.mat,
			Y.prefix, 
			Z.prefix,
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
			"GLS");

	}
	else
	{
		subset_op <- function(snpmat, sub.idx)
		{
			snp.sub <- get_plink_subsnp(snpmat, sub.idx );
			snp.mat <- cbind( snp.sub$info[,c(2,3)], snp.sub$snp )
			return( snp.mat );
		}
		
		r.gls <- snpmat_parallel(
			NCOL( pd$snp.mat$genotypes ),
			subset_op,
			pd$snp.mat,
			pd$phe.mat,
			Y.prefix, 
			Z.prefix,
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
			"GLS");
	}
	
	if(!is.null(r.gls) && !is.na(r.gls))
	{
		r <- wrap_GLS_ret(r.gls, r.filter, options);
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
}

gls.plink.tped<-function( file.phe, file.plink.tped, file.plink.tfam, Y.prefix, Z.prefix, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, options=NULL )       
{
	cat( "[ GLASSO PLINK(tped) ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.plink.tped) || missing(file.plink.tfam) ||
		missing(Y.prefix) || missing(Z.prefix) || missing(covar.names) )
		stop("! file.phe, file.plink.tped, file.plink.tfam, Y.prefix, Z.prefix and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
	if ( !missing("covar.names") && length(covar.names)>0 && !is.character(covar.names) )
		stop("! The parameter of covar.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(refit) && length(refit)==1 ) )
		stop("! The parameter of refit should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(add.used) && length(add.used)==1 ) )
		stop("! The parameter of add.used should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(dom.used) && length(dom.used)==1 ) )
		stop("! The parameter of dom.used should be a logical value(TRUE or FALSE).");
	#if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
	#	stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Data File = ",  file.phe, "\n");
	cat("* PLINK TPED File = ",  file.plink.tped, "\n");
	cat("* PLINK TFAM File = ",  file.plink.tfam, "\n");

	show_gls_parameters( Y.prefix, Z.prefix, covar.names, refit, add.used, dom.used, fgwas.filter=F ) ;

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
	r.gls <- .Call("gls_plink_tped", 
			file.phe,
			file.plink.tped, 
			file.plink.tfam, 
			Y.prefix, 
			Z.prefix, 
			paste(covar.names, collapse=","), 
			refit,
			add.used,
			dom.used,
			options$nMcmcIter,
			options$fBurnInRound,
			options$fRhoTuning,
			options$fQval.add,
			options$fQval.dom,
			ifelse( options$debug,3 , 1) );

	options$params <- list( file.phe = file.phe, 
				file.plink.tped= file.plink.tped, 
				file.plink.tfam= file.plink.tfam, 
				Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covar.names    = covar.names, 
				refit          = refit, 
				add.used       = add.used, 
				dom.used       = dom.used );


	if(!is.null(r.gls) && !is.na(r.gls))
	{
		r <- wrap_GLS_ret(r.gls, r.filter, options);
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}
			
	return(r);		   
}

gls.snpmat<-function( phe.mat, snp.mat, Y.prefix, Z.prefix, covar.names, refit=TRUE, add.used=TRUE, dom.used=TRUE, fgwas.filter=FALSE, options=NULL)     
{
	cat( "[ GLASSO SNPMAT ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(phe.mat) || missing(snp.mat) || missing(Y.prefix) || missing(Z.prefix) || missing(covar.names) )
		stop("! phe.mat, snp.mat, file.plink.tfam, Y.prefix, Z.prefix and covar.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
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

	show_gls_parameters( Y.prefix, Z.prefix, covar.names, refit, add.used, dom.used, fgwas.filter ) ;

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

	options$params <- list( Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covar.names    = covar.names, 
				refit          = refit, 
				add.used       = add.used, 
				dom.used       = dom.used, 
				fgwas.filter   = fgwas.filter);
	
	if( class(phe.mat)=="data.frame" )
	{
		cat("Phenotypic data frame is converted to the matrix class.\n");  
		phe.colnames <- colnames(phe.mat); 
		phe.rownames <- rownames(phe.mat); 
		phe.mat <- matrix(as.numeric(as.matrix(phe.mat, rownames.force=NA)), ncol=NCOL(phe.mat))
		colnames(phe.mat) <- phe.colnames;
		rownames(phe.mat) <- phe.rownames;
	}

	r.gls <- list();
	r.filter <- list();

	if( options$nPiecewise.ratio==0 && !fgwas.filter )
	{
		cat( "Genetic Effect Analysis by GLASSO method......\n");

		r.gls <- .Call("gls_snpmat", 
			as.matrix( phe.mat ),
			as.matrix( snp.mat*1.0 ),
			Y.prefix, 
			Z.prefix, 
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
			r.filter <- snpmat_fgwas_filter( phe.mat, snp.mat, Y.prefix, Z.prefix, covar.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");

			if( r.filter$error ) stop(r.filter$err.info);
			if( is.null(r.filter$snp.mat) ) return ( wrap_fgwas_ret( r.filter, options) );
		
			r.gls <- snpmat_parallel(
				NROW(r.filter$snp.mat),
				subset_op,
				r.filter$snp.mat,
				phe.mat,
				Y.prefix, 
				Z.prefix,
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
				"GLS");
		}
		else
		{
			r.gls <- snpmat_parallel(
				NROW(snp.mat),
				subset_op,
				snp.mat,
				phe.mat,
				Y.prefix, 
				Z.prefix,
				paste(covar.names, collapse=","), 
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
				"GLS");
		}
	}

	if(!is.null(r.gls) && !is.na(r.gls))
	{
		r <- wrap_GLS_ret(r.gls, r.filter, options);
		return(r);		   
	}
	else
	{
		cat("! No results\n");
		return(NULL);		   
	}

	return(r);
}

#summary_output2<-function(re1, re_add, re_dom)
#{
#	cat("(1) Covariate Estimate:\n");
#		
#	for(i in 1:NROW(re1))
#		cat(sprintf("%s \t %s \t %d%d%d%d \t%.3f\t(%.3f,%.3f,%.3f,%.3f)\n", 
#			rownames(re1)[i], ifelse(re1[i,1], "Yes", "---"), 
#			re1[i,1], re1[i,2], re1[i,3], re1[i,4], 
#			
#			re1[i,5], 
#			
#			re1[i,6], re1[i,7], re1[i,8], re1[i,9] ));  
#
#	cat("(2) Significant SNPs Estimate:\n");
#
#	cat("    SNP Name\tGrp/Pos\tAdd\tLR2(Median1,2,3,4)\t\tDom\tLR2(Median1,2,3,4)\n");
#	for(i in 1:NROW(re_add))
#	{
#		cat(sprintf("%12s\t%d/%d\t", rownames(re_add)[i], re_add[i,1], re_add[i,2] ));
#  	        cat(sprintf("%s(%d%d%d%d)\t%.3f(%.3f,%.3f,%.3f,%.3f)", ifelse(sum(re_add[i,3:6])>0, "Yes/", "---/"),
#		    re_add[i,3], re_add[i,4], re_add[i,5], re_add[i,6], 
#		    re_add[i,7], re_add[i,8], re_add[i,9], re_add[i,10], re_add[i,11])) 
#
#  	        cat(sprintf("%s(%d%d%d%d)\t%.3f(%.3f,%.3f,%.3f,%.3f)", ifelse(sum(re_dom[i,3:6])>0, "Yes/", "---/"),
#		    re_dom[i,3], re_dom[i,4], re_dom[i,5], re_dom[i,6], 
#		    re_dom[i,7], re_dom[i,8], re_dom[i,9], re_dom[i,10], re_dom[i,11])) 
#		
#		cat("\n");  
#	}	
#
#	if(!is.null(r.gls$refit_add))
#	{
#		cat("Refit Result\n");
#		summary_output2(r.gls$refit_cov, r.gls$refit_add, r.gls$refit_dom)
#	}
#	else if(!is.null(r.gls$varsel_add))
#	{
#		cat("Variable Selection Result\n");
#		summary_output2(r.gls$varsel_cov, r.gls$varsel_add, r.gls$varsel_dom); 
#	}
#}


merge_add_dom<-function( re_add, re_dom )
{
	if( is.null(re_add) && is.null(re_dom) ) 
		return(NULL);

	idx.add <- c();
	if(!is.null(re_add))
		idx.add <- which( rowSums( re_add[,c(3:6),drop=F] ) > 0  )
	idx.dom <- c();
	if(!is.null(re_dom))
		idx.dom <- which( rowSums( re_dom[,c(3:6),drop=F] ) > 0  )

	idx.sig <- unique(c( idx.add, idx.dom ));	
	
	if( length(idx.sig)==0 ) return(NULL);
	
	sig.add <- NULL;
	if(!is.null(re_add)) sig.add <- re_add[ idx.sig, c(1:11), drop=F];
	sig.dom <- NULL;
	if(!is.null(re_dom)) sig.dom <- re_dom[ idx.sig, c(1:11), drop=F];
	
	sig.mat <- c();
	if( !is.null(sig.dom) && !is.null(sig.add) )
		sig.mat <- data.frame(sig.add[,c(1,2),drop=F], 
				Add.Sig = rowSums( sig.add[,c(3:6),drop=F]), 
				sig.add[,c(7:11),drop=F], 
				Dom.Sig = rowSums( sig.dom[,c(3:6),drop=F]), 
				sig.dom[,c(7:11),drop=F] )
	else if( !is.null(sig.add) )
		sig.mat <- data.frame(sig.add[,c(1,2),drop=F], 
				Add.Sig = rowSums( sig.add[,c(3:6),drop=F]), 
				sig.add[,c(7:11),drop=F], 
				array(NA, dim=c(length(idx.sig), 6)) )
	else
		sig.mat <- data.frame( sig.dom[,c(1,2),drop=F], 
				array( NA, dim=c(length(idx.sig), 6) ), 
				Dom.Sig = rowSums( sig.dom[, c(3:6), drop=F] ), 
				sig.dom[,c(7:11),drop=F] );
				
	colnames(sig.mat) <- c( "chr", "pos", 
				"add.sig", "add.mode", "add.mu1", "add.mu2", "add.mu3", "add.mu4",
				"dom.sig", "dom.mode", "dom.mu1", "dom.mu2", "dom.mu3", "dom.mu4" );

	if( is.null(sig.dom) )				
		rownames(sig.mat) <- rownames( sig.add)
	else	
		rownames(sig.mat) <- rownames( sig.dom);

	return(sig.mat);
}

summary.GLS.ret<-function(object, ...)
{
	r.gls <- object;
	r.sum.ret <- list();

	if(!is.null( r.gls$refit_cov ) && NROW( r.gls$refit_cov )>0 )
	{
		re1 <- r.gls$refit_cov;
		r.sum.ret$refit_cov <- data.frame(  
					Mode  = round(re1[,5], digits=3), 
					L1    = round(re1[,6], digits=3), 
					L2    = round(re1[,7], digits=3), 
					L3    = round(re1[,8], digits=3), 
					L4    = round(re1[,9], digits=3) );
		rownames(r.sum.ret$refit_cov) <- rownames( re1 );
	}

	r.add <- !is.null( r.gls$refit_add ) && NROW( r.gls$refit_add )>0;
	r.dom <- !is.null( r.gls$refit_dom ) && NROW( r.gls$refit_dom )>0;
	if( r.add || r.dom )
		r.sum.ret$refit <- merge_add_dom( r.gls$refit_add, r.gls$refit_dom);
	
	if(!is.null( r.gls$varsel_cov ) && NROW( r.gls$varsel_cov )>0 )
	{
		re3 <- r.gls$varsel_cov;
		r.sum.ret$varsel_cov <- data.frame(  
					Mode  = round(re3[,5], digits=3), 
					L1    = round(re3[,6], digits=3), 
					L2    = round(re3[,7], digits=3), 
					L3    = round(re3[,8], digits=3), 
					L4    = round(re3[,9], digits=3) );
		rownames(r.sum.ret$varsel_cov) <- rownames( re3 );
	}
	
	r.add <- !is.null( r.gls$varsel_add ) && NROW( r.gls$varsel_add )>0;
	r.dom <- !is.null( r.gls$varsel_dom ) && NROW( r.gls$varsel_dom )>0;
	if( r.add || r.dom )
		r.sum.ret$varsel <- merge_add_dom( r.gls$varsel_add, r.gls$varsel_dom);

	if(!is.null(r.gls$fgwas))
	{
		re7 <- r.gls$fgwas;
		fgwas.sig <- which( re7[,7] <= r.gls$options$fgwas.cutoff );
		if(length(fgwas.sig)>0)
		{
			fgwas_sigs <- re7[ fgwas.sig, , drop=F];
			fgwas.sig.inc <- order(fgwas_sigs[,7]);
			r.sum.ret$fgwas_sig <- fgwas_sigs[fgwas.sig.inc,];
		}
		
		if(!is.null(r.sum.ret$varsel))
			r.sum.ret$varsel <- cbind(r.sum.ret$varsel, fgwas.pvalue=find_fgwas_pvalue( r.gls$fgwas, rownames(r.sum.ret$varsel) ) ) ;

		if(!is.null(r.sum.ret$refit))
			r.sum.ret$refit <- cbind(r.sum.ret$refit, fgwas.pvalue=find_fgwas_pvalue( r.gls$fgwas, rownames(r.sum.ret$refit) ) ) ;
		
	}

	class(r.sum.ret) <- "sum.GLS.ret";
	
	r.sum.ret
}

print.sum.GLS.ret<-function(x, ...)
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
		cat("--- Variable Selection Result:", NROW(r.sum.ret$varsel), "SNPs\n");
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
		cat("--- Refit Result:", NROW(r.sum.ret$refit), "SNPs\n");
		show(r.sum.ret$refit);
	}
}

wrap_fgwas_ret<-function( r.filter, options )
{
	cat( "Wrapping the results ......\n");
	
	r.gls <- list();
	if(!is.null(r.filter)) r.gls$fgwas <- r.filter$r;
	r.gls$options <- options;
	class(r.gls) <- "GLS.ret";

	return(r.gls);
}

wrap_GLS_ret<-function(r.gls, r.filter, options )
{
	cat( "Wrapping the results ......\n");

	add_c19 <- c("grp", "pos", "add.m1", "add.m2", "add.m3", "add.m4",
		"add_L2.mu", "add.mu1", "add.mu2", "add.mu3", "add.mu4",
		"add_L2.min", "add.min1", "add.min2", "add.min3", "add.min4",
		"add_L2.max", "add.max1", "add.max2", "add.max3", "add.max4" );
	dom_c19 <- c("grp", "pos", "dom.m1", "dom.m2", "dom.m3", "dom.m4",
		"dom_L2.mu", "dom.mu1", "dom.mu2", "dom.mu3", "dom.mu4",
		"dom_L2.min", "dom.min1", "dom.min2", "dom.min3", "dom.min4",
		"dom_L2.max", "dom.max1", "dom.max2", "dom.max3", "dom.max4" );
	cov_c19 <- c("cov.m1", "cov.m2", "cov.m3", "cov.m4",
		"co_L2.mu", "cov.mu1", "cov.mu2", "cov.mu3", "cov.mu4",
		"co_L2.min", "cov.min1", "cov.min2", "cov.min3", "cov.min4",
		"co_L2.max", "cov.max1", "cov.max2", "cov.max3", "cov.max4" );

	if(!is.null(r.gls) && !is.na(r.gls))
	{
		if (!is.null(r.gls$varsel_add)) colnames(r.gls$varsel_add) <- add_c19;
		if (!is.null(r.gls$varsel_dom)) colnames(r.gls$varsel_dom) <- dom_c19;
		if (!is.null(r.gls$varsel_cov)) colnames(r.gls$varsel_cov) <- cov_c19;
		if (!is.null(r.gls$refit_add))  colnames(r.gls$refit_add) <- add_c19;
		if (!is.null(r.gls$refit_dom))  colnames(r.gls$refit_dom) <- dom_c19;
		if (!is.null(r.gls$refit_cov))  colnames(r.gls$refit_cov) <- cov_c19;
		
		row.cov <- c("Intercept", options$params$covar.names);
		if( !is.null(r.gls$varsel_cov) && NROW(r.gls$varsel_cov)>0)
			rownames(r.gls$varsel_cov) <- row.cov;
		
		if( !is.null(r.gls$refit_cov) && NROW(r.gls$refit_cov)>0)
			rownames(r.gls$refit_cov) <- row.cov;
	
		if(!is.null(r.gls$varsel_PSRF))
		{
			colnames(r.gls$varsel_PSRF) <- r.gls$varsel_PSRF[1,];
			r.gls$varsel_PSRF <- r.gls$varsel_PSRF[-1,]
		}	

		if(!is.null(r.gls$refit_PSRF))
		{
			colnames(r.gls$refit_PSRF) <- r.gls$refit_PSRF[1,];
			r.gls$refit_PSRF <- r.gls$refit_PSRF[-1,]
			
		}	

		if(!is.null(r.filter)) r.gls$fgwas <- r.filter$r;
		
		r.gls$options <- options;

		class(r.gls) <- "GLS.ret";
	}		
	
	return(r.gls);
}

show_gls_parameters<-function( Y.prefix, Z.prefix, covar.names, refit, add.used, dom.used, fgwas.filter ) 
{
	cat( "* Response Variable =",   Y.prefix, "\n");
	cat( "* Time Variable =",   Z.prefix, "\n");
	cat( "* Covariate Columns =",  covar.names, "\n");
	cat( "* fGWAS Filter Used =",  ifelse( fgwas.filter, "Yes", "No"), "\n");
	cat( "* Additive Effects Used =",  ifelse( add.used, "Yes", "No"), "\n");
	cat( "* Dominant Effects Used =",  ifelse( dom.used, "Yes", "No"), "\n");
	cat( "* Refit Procedure =",   ifelse( refit, "Yes", "No"), "\n");
}

get_sig_gls_snp <- function( r.gls )
{
	if( is.null(r.gls$varsel_add) && is.null(r.gls$varsel_dom) ) return( NULL );

	idx.sig.add <- c();
	if( !is.null(r.gls$varsel_add) )
		idx.sig.add <- which( rowSums(r.gls$varsel_add[,c(3:6),drop=F]) > 0 );

	idx.sig.dom <- c();
	if( !is.null(r.gls$varsel_dom) )
		idx.sig.dom <- which( rowSums(r.gls$varsel_dom[,c(3:6),drop=F]) > 0 );

	idx.sig <- unique(c(idx.sig.dom, idx.sig.add));
	if (length( idx.sig )==0) return(NULL);
	
	return(idx.sig);
}


plot.GLS.ret<-function( x, y=NULL, ... , fig.prefix=NULL )
{
	r.gls <- x;

	if( missing(fig.prefix)) fig.prefix <- "gls.plot";

	if(!is.null(r.gls$fgwas))
	{
		filter.man <- r.gls$fgwas[, c(1,2,7), drop=F]
		draw_man_fgwas( filter.man, fig.prefix, "fgwas" );
	}
	else
		cat("! No fGWAS filter results.\n");		
		
	if( !is.null(r.gls$varsel_add) || !is.null(r.gls$varsel_dom))
	{
		if ( !is.null(r.gls$varsel_add) )  varsel <- r.gls$varsel_add[, c(1,2), drop=F]
		if ( !is.null(r.gls$varsel_dom) )  varsel <- r.gls$varsel_dom[, c(1,2), drop=F]

		if ( !is.null(r.gls$varsel_add) ) varsel<- cbind( varsel, r.gls$varsel_add[,7] );
		if ( !is.null(r.gls$varsel_dom) ) varsel<- cbind( varsel, r.gls$varsel_dom[,7] );

		draw_man_adh2( varsel, fig.prefix, "varsel" );
	}
	else
		cat("! No varible selection results.\n");		

	if( !is.null(r.gls$refit_add) || !is.null(r.gls$refit_dom) )
	{
		refit<- merge_add_dom( r.gls$refit_add, r.gls$refit_dom);

		draw_refit_curve( refit, fig.prefix, "curve" );
	}
	else
		cat("! No refit results.\n");		
}

read_simple_gls_data <- function( file.phe, file.snp, bImputed=TRUE )
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

gls.best.qval<-function( r.gls, snp.names )
{
	if( class(r.gls) != "GLS.ret" )	
	{
		cat("! r.bls is NOT a GLS.ret object.\n");
		return(NULL);
	}

	find_Qtable<-function( Q.table, snp.names)
	{
		pv.mat4 <- array( NA, dim=c(length(snp.names),4 ) );
		pv <- rep( NA, length(snp.names) );
	
		idx <- match( snp.names, rownames(Q.table));
		if( length(which(!is.na(idx)) ) > 0)
		{
			pv.mat4[!is.na(idx), ] <-  Q.table[ idx[!is.na(idx)], c(1:4) ];
			pv[!is.na(idx)] <- apply(pv.mat4[!is.na(idx), ], 1, min, na.rm=TRUE)
		}
		
		return(pv);
	}
	
	Qval.vs.add <- rep( NA, length(snp.names) );
	Qval.vs.dom <- rep( NA, length(snp.names) );
	if( is.null(r.gls$varsel_Qbest ) )
	{
		cat("! No Qbest matrix for variable selection.\n");
	}
	else
	{
		if(any(r.gls$varsel_Qbest[,c(1:4)]!=0))
			Qval.vs.add <- find_Qtable( r.gls$varsel_Qbest[,  c(1:4),drop=F], snp.names );
		if(any(r.gls$varsel_Qbest[,c(20:23)]!=0))
			Qval.vs.dom <- find_Qtable( r.gls$varsel_Qbest[,c(20:23),drop=F], snp.names );
	}
	
	Qval.refit.add <- rep( NA, length(snp.names) );
	Qval.refit.dom <- rep( NA, length(snp.names) );
	if( is.null(r.gls$refit_Qbest ) )
	{
		cat("! No Qbest matrix for refit procedure.\n");
	}
	else
	{
		if(any(r.gls$refit_Qbest[,c(1:4)]!=0))
			Qval.refit.add <- find_Qtable( r.gls$refit_Qbest[,c(1:4),drop=F], snp.names );
		if(any(r.gls$refit_Qbest[,c(20:23)]!=0))
			Qval.refit.dom <- find_Qtable( r.gls$refit_Qbest[,c(20:23),drop=F], snp.names );
	}

	fgwas.pv <- find_fgwas_pvalue( r.gls$fgwas, snp.names)

	r.qval <- cbind( fgwas.pv, Qval.vs.add, Qval.vs.dom, Qval.refit.add, Qval.refit.dom);
	colnames(r.qval) <- c("fgwas.pv", "Qval.vs.add", "Qval.vs.dom", "Qval.refit.add", "Qval.refit.dom");
	rownames(r.qval) <- snp.names;
	
	return(r.qval);
}

gls.qval.cutoff<-function( r.gls, qval.add, qval.dom, refit.select = FALSE )
{
	if( class(r.gls) != "GLS.ret" )	
	{
		cat("! r.gls is NOT a GLS.ret object.\n");
		return(NULL);
	}

	find_Qtable<-function( Q.table, qval)
	{
		pv <- apply(Q.table[, c(1:4) ], 1, min, na.rm=TRUE)
		return(which(pv<=qval));
	}
	
	if( is.null(r.gls$varsel_Qbest ) && !refit.select )
	{
		cat("! No Qbest matrix for variable selection.\n");
	}
	else
	{
		idx.vs.add <- c();
		idx.vs.dom <- c();

		if(any(r.gls$varsel_Qbest[,c(1:4)]!=0))
			idx.vs.add <- find_Qtable( r.gls$varsel_Qbest[,  c(1:4),drop=F], qval.add );
		if(any(r.gls$varsel_Qbest[,c(20:23)]!=0))
			idx.vs.dom <- find_Qtable( r.gls$varsel_Qbest[,c(20:23),drop=F], qval.dom );
		
		idx.vs <- sort( unique(c(idx.vs.add, idx.vs.dom)) )
		if(length(idx.vs)>0)
			return( rownames(r.gls$varsel_Qbest)[idx.vs])
		else
			return(NULL);			
	}
	
	if( is.null(r.gls$refit_Qbest )  && refit.select )
	{
		cat("! No Qbest matrix for refit procedure.\n");
	}
	else
	{
		idx.refit.add <- c();
		idx.refit.dom <- c();

		if(any(r.gls$refit_Qbest[,c(1:4)]!=0))
			idx.refit.add <- find_Qtable( r.gls$refit_Qbest[,c(1:4),drop=F], qval.add );
		if(any(r.gls$refit_Qbest[,c(20:23)]!=0))
			idx.refit.dom <- find_Qtable( r.gls$refit_Qbest[,c(20:23),drop=F], qval.dom );

		idx.refit <- sort( unique(c(idx.refit.add, idx.refit.dom)) )
		if(length(idx.refit)>0)
			return( rownames(r.gls$refit_Qbest)[idx.refit])
		else
			return(NULL);			
	}

	return(NULL);
}
