get_nodup_rownames<-function(r.org, r.new)
{
	row.inter <- intersect( rownames(r.new), rownames(r.org) );
	if( length(row.inter)>0 )
	{
		dup.pos <- match( row.inter, rownames(r.new) );
		return( r.new[ -dup.pos,,drop=F] );
	}
	else
		return( r.new );
}

merge_xls_varsel<-function( r.old, r.cluster )
{
	r.new <- list( varsel = c(),
		varsel_add = c(),
		varsel_dom = c(),
		varsel_Qbest = c(),
		varsel_PSRF = c() );

	r.new[ names(r.old) ] <- r.old;

	for(i in 1:length(r.cluster) )
	{
		if ( !is.null(r.cluster[[i]]$varsel) )
		{
			r.nodup <- get_nodup_rownames( r.new$varsel, r.cluster[[i]]$varsel);
			r.new$varsel <- rbind( r.new$varsel,  r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_add) )
		{
			r.nodup <- get_nodup_rownames( r.new$varsel_add, r.cluster[[i]]$varsel_add);
			r.new$varsel_add <- rbind( r.new$varsel_add, r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_dom) )
		{
			r.nodup <- get_nodup_rownames( r.new$varsel_dom, r.cluster[[i]]$varsel_dom);
			r.new$varsel_dom <- rbind( r.new$varsel_dom, r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_Qbest) )
		{
			r.nodup <- get_nodup_rownames( r.new$varsel_Qbest, r.cluster[[i]]$varsel_Qbest);
			r.new$varsel_Qbest <- rbind( r.new$varsel_Qbest, r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_PSRF) )
		{
			r.nodup <- get_nodup_rownames( r.new$varsel_PSRF, r.cluster[[i]]$varsel_PSRF);
			r.new$varsel_PSRF <- rbind( r.new$varsel_PSRF,r.nodup );
		}
	}

	return(r.new);
}


get_snpmat_byname <- function( n.snp, snp.mat, varsel.snpname, f_subset_op)
{
	varsel.snpmat <- c();

	for(i.sect in 1:ceiling(n.snp/1000))
	{
		sub.set <- (i.sect-1)*1000 + c(1:1000);
		sub.set<- sub.set[which( sub.set<=n.snp)];

		sub.snp <- f_subset_op( snp.mat, sub.set );

		sub.snp.idx <- match( varsel.snpname, rownames(sub.snp) );
		sub.snp.idx <- sub.snp.idx[!is.na(sub.snp.idx)]
		if (length(sub.snp.idx) >0 ) varsel.snpmat <- rbind(varsel.snpmat, sub.snp[sub.snp.idx, ,drop=F]);
	}


	return( varsel.snpmat );
}

snpmat_parallel<-function( n.snp,
			f_subset_op,
			snp.mat,
			phe.mat,
			Y.name,
			Z.name,
			covar.names,
			refit,
			gpu.used,
			add.used,
			dom.used,
			op.piecewise.ratio,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			op.debug,
			op.cpu,
			lasso.method)
{
	cat( "Genetic Effect Analysis by BLASSO/GLASSO method......\n");

	n.inv <- NROW(phe.mat);

	if( op.piecewise.ratio == 0 || n.snp < n.inv * op.piecewise.ratio*1.2 )
	{

		cat("* Final LASSO calling.\n");

		snpmat <- f_subset_op(snp.mat, c(1:n.snp) );

		r.xls <- snpmat_call( snpmat,
				phe.mat,
				Y.name,
				Z.name,
				covar.names,
				refit,
				gpu.used,
				add.used,
				dom.used,
				op.nMcmcIter,
				op.fBurnInRound,
				op.fRhoTuning,
				op.fQval.add,
				op.fQval.dom,
				op.debug,
				lasso.method);
		return(r.xls);
	}

	real.task <- ifelse( op.cpu>0, op.cpu*2, 1);
	snp.sect0 <- seq(1, n.snp, n.inv * op.piecewise.ratio);
	snp.sect1 <- c( snp.sect0[-1]-1, n.snp );
	sect.list <- seq( 1, length(snp.sect0), real.task );

	r.cluster.init <- list( varsel = c(), varsel_add = c(), varsel_dom = c(), varsel_Qbest = c(), varsel_PSRF = c() );

	for(i in sect.list )
	{
		snpmat.list <- list();
		for(k in 1:real.task )
		{
			if( i + k - 1 <= length(snp.sect0) )
			{
				idx.snp <- c(snp.sect0[ i + k - 1]:snp.sect1[ i + k - 1 ]);
				if ( length(idx.snp) < n.inv * op.piecewise.ratio )
				{
					n.add <- n.inv * op.piecewise.ratio - length(idx.snp);
					idx.snp0.add <- sample( 1:snp.sect0[ i + k - 1] ) [1:n.add];
					idx.snp0.add <- idx.snp0.add[ !is.na(idx.snp0.add) ];
					idx.snp <- c(idx.snp, idx.snp0.add);
				}

				snpmat.list[[k]] <- f_subset_op(snp.mat, idx.snp );
			}
		}

		r.cluster <- snpmat_parallel_list( phe.mat,
					snpmat.list,
					Y.name,
					Z.name,
					covar.names,
					refit=F,
					gpu.used,
					add.used,
					dom.used,
					op.nMcmcIter,
					op.fBurnInRound,
					op.fRhoTuning,
					op.fQval.add * 1.5, # For 1st run, more SNPs are selected than the pre-Q value
					op.fQval.dom * 1.5,
					op.debug,
					op.cpu,
					lasso.method);

		r.cluster.init <- merge_xls_varsel( r.cluster.init, r.cluster );
	}

	if(lasso.method=="BLS")
		idx.sig <- get_sig_bls_snp( r.cluster.init )
	else
		idx.sig <- get_sig_gls_snp( r.cluster.init );

	if( is.null(idx.sig) )
	{
		cat("! No SNPs are selected in the first run. \n");
		return(r.cluster.init);
	}

	cat("*", length(idx.sig), "SNPs are selected in the first run.\n");

	# ! One reason leads to select the SNP.mat by the name, not easily index(idx.sig)
	# !  Some SNPs are removed by th C code becuase of MAF too small or missing too many.
	# !  so index is not a safe way to proceed the variable selection.

	varsel.snpname <- c();
	# for BLS functions, varsel can be used, but for GLS functions, varsel_add or varsel_dom can be used
	if( !is.null(r.cluster.init$varsel) )
		varsel.snpname <- rownames(r.cluster.init$varsel)[idx.sig]
	else
	if( !is.null(r.cluster.init$varsel_add) )
		varsel.snpname <- rownames(r.cluster.init$varsel_add)[idx.sig]
	else
	if( !is.null(r.cluster.init$varsel_dom) )
		varsel.snpname <- rownames(r.cluster.init$varsel_dom)[idx.sig]

	varsel.snpmat <- get_snpmat_byname( n.snp, snp.mat, varsel.snpname, f_subset_op);

	## DOUBEL check for debug problem
	if(length(idx.sig) != NROW(varsel.snpmat))
	{
		cat("!", NROW(varsel.snpmat), "SNPs are selected in the first run.\n");
		if( NROW(varsel.snpmat) <1 )
		{
			#No SNPs are selected in the first run, exit with the variable selection results.
			return(r.cluster.init);
		}
	}

	R <- 1;
	while( NROW(varsel.snpmat) > n.inv * op.piecewise.ratio*1.2 )
	{
		n.snp0 <- NROW(varsel.snpmat);
		snp.sect0 <- seq(1, n.snp0, n.inv * op.piecewise.ratio);
		snp.sect1 <- c( snp.sect0[-1]-1, n.snp0 );
		sect.list <- seq(1, length(snp.sect0), real.task );

		r.cluster0 <- list( varsel = c(), varsel_add = c(), varsel_dom = c(), varsel_Qbest = c(), varsel_PSRF = c() );
		for(i in sect.list )
		{
			snpmat.list <- list();
			for (k in 1:length(snp.sect0))
			{
				if( i + k - 1 > length(snp.sect0) )
					next;

				idx.snp <- snp.sect0[ i + k - 1]:snp.sect1[i + k -1];
				if ( length(idx.snp) < n.inv * op.piecewise.ratio )
				{
					n.add <- n.inv * op.piecewise.ratio - length(idx.snp);
					idx.snp.add <- sample(1:snp.sect0[ i + k - 1])[1:n.add];
					idx.snp0.add <- idx.snp0.add[ !is.na(idx.snp0.add) ];
					idx.snp <- c(idx.snp, idx.snp.add );
				}

				snpmat.list[[k]] <- varsel.snpmat[ idx.snp, ,drop=F]
			}

			r.cluster <- snpmat_parallel_list( phe.mat,
						snpmat.list,
						Y.name,
						Z.name,
						covar.names,
						refit=F,
						gpu.used,
						add.used,
						dom.used,
						op.nMcmcIter,
						op.fBurnInRound,
						op.fRhoTuning,
						op.fQval.add,
						op.fQval.dom,
						op.debug,
						op.cpu,
						lasso.method);


			r.cluster0 <- merge_xls_varsel( r.cluster0, r.cluster );
		}


		if(lasso.method=="BLS")
			idx.sig <- get_sig_bls_snp( r.cluster0 )
		else
			idx.sig <- get_sig_gls_snp( r.cluster0 );

		if( is.null(idx.sig) )
		{
			cat("! No SNPs are selected in the ", R, "(th) run.\n");
			return(r.cluster.init);
		}

		varsel.snpmat <- varsel.snpmat[idx.sig, ,drop=F];

		R <- R + 1;
		cat("*", NROW(varsel.snpmat), "SNPs are selected in ", R, "(th) run.\n")
	}

	cat("* Final LASSO calling.\n");
	r.xls <- snpmat_call(
			varsel.snpmat,
			phe.mat,
			Y.name,
			Z.name,
			covar.names,
			refit,
			gpu.used,
			add.used,
			dom.used,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			op.debug,
			lasso.method);
	# r.xls$varsel
	# r.xls$varsel_add
	# r.xls$varsel_dom
	# r.xls$varsel_Qbest
	# r.xls$varsel_PSRF
	r.xls[ names(r.cluster.init) ] <- r.cluster.init;

	return(r.xls);
}

snpmat_parallel_list<-function( phe.mat,
				snpmat.list,
				Y.name,
				Z.name,
				covar.names,
				refit=F,
				gpu.used,
				add.used,
				dom.used,
				op.nMcmcIter,
				op.fBurnInRound,
				op.fRhoTuning,
				op.fQval.add,
				op.fQval.dom,
				op.debug,
				op.ncpu,
				lasso.method)

{
	execInSnow <- FALSE;

	cpu.fun<-function( snpmat )
	{
		library(gwas.lasso);
		if (execInSnow)
		{
			library(snowfall);
			sfCat( paste("snpmat[", NROW(snpmat), NCOL(snpmat), "]", "\n") );
		}

		r.xls.i <- snpmat_call(
				#as.matrix( snpmat.list[[ sect ]]),
				as.matrix( snpmat ),
				phe.mat,
				Y.name,
				Z.name,
				covar.names,
				refit=F,
				gpu.used,
				add.used,
				dom.used,
				op.nMcmcIter,
				op.fBurnInRound,
				op.fRhoTuning,
				op.fQval.add,
				op.fQval.dom,
				op.debug,
				lasso.method);

		return(r.xls.i);
	}

	r.cluster <- list();

	if( op.ncpu>1 && require("snowfall") )
	{
		cat("Starting parallel computing, snowfall/snow......\n");
		sfInit(parallel = TRUE, cpus = op.ncpu, type = "SOCK" )

		execInSnow <- TRUE;

		sfExport("phe.mat", "Y.name", "Z.name", "covar.names",
				##"snpmat.list" ,
				"refit",
				"gpu.used",
				"add.used",
				"dom.used",
				"op.nMcmcIter",
				"op.fBurnInRound",
				"op.fRhoTuning",
				"op.fQval.add",
				"op.fQval.dom",
				"op.debug",
				"execInSnow",
				"lasso.method");

		#r.cluster <- sfClusterApplyLB( 1:length(snpmat.list), cpu.fun);
		r.cluster <- sfLapply( snpmat.list, cpu.fun);

		sfStop();

		cat("Stopping parallel computing......\n");
	}
	else
	{
		cat("Starting piecewise analysis......\n");
		for(i in 1:length(snpmat.list))
			r.cluster[[i]] <- cpu.fun( snpmat.list[[i]] );
	}

	return( r.cluster );
}

snpmat_call<-function(  snp.mat,
			phe.mat,
			Y.name,
			Z.name,
			covar.names,
			refit,
			gpu.used,
			add.used,
			dom.used,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			op.debug,
			lasso.method)
{
	if(lasso.method=="BLS")
	{
		r <- .Call("bls_snpmat",
			as.matrix( phe.mat ),
			as.matrix( snp.mat*1.0 ),
			Y.name,
			paste(covar.names, collapse=","),
			refit,
			gpu.used,
			add.used,
			dom.used,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			ifelse( op.debug, 3, 1) );
	}
	else
	{
		r <- .Call("gls_snpmat",
			as.matrix( phe.mat ),
			as.matrix( snp.mat*1.0  ),
			Y.name,
			Z.name,
			paste(covar.names, collapse=","),
			refit,
			gpu.used,
			add.used,
			dom.used,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			ifelse(op.debug, 3, 1));
	}

	return(r);
}

