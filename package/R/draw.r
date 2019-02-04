draw_man_adh2<-function( adh2, fig.prefix=NULL, fig.name=NULL )
{
	n.subfig <- NCOL(adh2) -2;
	pdf.file <- "";
	if (!is.null( fig.prefix ))
	{
		pdf.file <- paste( fig.prefix, fig.name, "pdf", sep=".");
		
		pdf( pdf.file, width=6, height= 2 * n.subfig );
	}

	draw_snplist<-function( values, yLabel, cex=1.0, col="black")
	{
		xlim <- c(0, length(values));
		ylim <- range(values);
		nB <- length(values);
		snps <- c(0:nB);
		plot(1:10,1:10, xlim=xlim, ylim=ylim, type="n", xlab="SNP", ylab=yLabel, cex.axis=cex );
		
		rect(snps[-(nB+1)], 0, snps[-1L], values,  col = col, border = col);
	}

	draw_call<-function()
	{
		if( n.subfig >= 3)
			par(mar=c(4.5, 4, 0.5, 2) + 0.1)
		else
		{
			par(mar=c(3, 3.5, 0.5, 2) + 0.1);
			par(mgp = c(1.6, 0.5, 0 ) );
		}
		
		plot.new();
		par(mfrow=c( NCOL(adh2)-2 , 1 ) );

		par(mfg=c(1, 1));
		draw_snplist( unlist(adh2[,3, drop=F ]), "Additive effects", cex=ifelse( n.subfig >= 3, 1.0, 0.6), col="orange" );

		if( n.subfig >= 2)
		{
			par(mfg=c(2, 1));
			draw_snplist( unlist(adh2[,4, drop=F ]), "Dominant effects", cex=ifelse( n.subfig >= 3, 1.0, 0.6), col="purple" );
		}
	
		if( n.subfig >= 3)
		{
			par(mfg=c(3, 1));
			draw_snplist( unlist(adh2[,5, drop=F ]), "Heritability", cex=ifelse( n.subfig >= 3, 1.0, 0.6) );
		}	
	}

	# sort adh2 by grpup number and position.
	adh2 <- adh2[order(adh2[,1], adh2[,2]),,drop=F]
	
	r <- try( draw_call() );
	if(class(r) == "try-error")
		cat("Failed to draw figure.\n")
	
	if (!is.null( fig.prefix ))
	{
		dev.off();
		cat("* The figure is exported into", pdf.file, "\n");
	}		
}

draw_fgwas_manhattan<-function( res, map.title="", sig.thres, dot.cex, y.max=NA)
{
	pvalues <- -log10(res[,3]);
	if( length( which(is.infinite(pvalues)) )>0 )
	{
		pv.temp <- pvalues[ -which(is.infinite(pvalues)) ];
		pvalues[ which(is.infinite(pvalues)) ] <- max(pv.temp)*1.5;
	}
	
	nrow    <- dim(res)[1];
	log10.max <- round(max(pvalues, na.rm=T))+1;

	if (is.na(y.max))
		y.max <- log10.max;

	if (length(which(pvalues>y.max))>0)
		pvalues[which(pvalues>y.max)] <- y.max;

	par( xaxs="r",  yaxs="r");  # Extend axis limits by 4%
	par( tck = -0.04 );
	par( mgp = c(1.6, 0.5, 0 ) );
	par( mar=c( 2.5, 3.5, 2, 1.5)+0.1);

	plot( 1,1, type="n", xlab="SNP", ylab=expression(-log[10](italic(p))),
		  cex.axis=0.7, xlim=c(1, nrow), ylim=c(0,  y.max), xaxt="n", main=map.title );


	p.lab <-  - c( log10( sig.thres) );

	abline( h= c(p.lab), col="gray", lwd=1, lty="dashed");

	#text( x=0, p.lab[1] + 0.1, "p=0.01", cex=0.8, srt=90, adj=c(0.5, -1)); 

	cols <- c( "darkgreen","black",  "orange",  "red", "blue", "purple");
	
	p.cex <- rep(0.5*dot.cex, length(pvalues));
	p.cex[ which(pvalues>0.4*y.max) ] <- 0.5*dot.cex + ( pvalues[which(pvalues>y.max*0.4)]-0.4*y.max)/y.max*dot.cex;

	points( pvalues, pch=20, col = cols[ ( res[,1] %% 6 + 1 ) ], cex=p.cex);

	x.off <- 0;
	x <- c();
	x.ps <- c();
	for(i in 1:18)
	{
		x.ps<-c( x.ps, length(which(res[,1]==i) ) );
		x <- c(x, x.ps[i]/2 + x.off);
		x.off <- x.off + x.ps[i];
	}

	x.ps<-c( x.ps, length(which(res[,1]>18 ) )  )
	x <- c(x, x.ps[19]*2/3 + x.off);
	
	axis(side=1, at=x, col="black", labels=c(paste("", 1:18, sep=""),"..."), col.axis="black", col.lab="black", cex.axis=0.4, padj=-0.2 )
}	

draw_man_fgwas<-function( r.fgwas, fig.prefix=NULL, fig.name=NULL )
{
	pdf.file <- NULL;

	if (!is.null(fig.prefix ))
	{
		pdf.file <- paste(fig.prefix, fig.name, "pdf", sep=".");
		pdf( pdf.file, width=6, height=3);
	}

	# sort r.fgwas by grpup number and position.
	r.fgwas <- r.fgwas[ order(r.fgwas[,1], r.fgwas[,2]),,drop=F]

	r <- try( draw_fgwas_manhattan( r.fgwas, map.title="fGWAS Analysis", 0.05/NROW(r.fgwas), 0.7 ) );
	if(class(r) == "try-error")
		cat("Failed to draw figure.\n")
	
	if (!is.null( pdf.file ))
	{
		dev.off();
		cat("* The figure is exported into", pdf.file, "\n");
	}		
}

#--------------------------------------------------------------
# draw_refit_curve
# 
#  colnames(refit) <=c( [1] "chr", [2]"pos", 
#			[3] "add.sig", [4] "add.mode", [5]"add.mu0", "add.mu1", "add.mu2", "add.mu3",
#			[9] "dom.sig", [10]"dom.mode", [11]"dom.mu0", "dom.mu1", "dom.mu2", "dom.mu3" );
#
#--------------------------------------------------------------

draw_refit_curve<-function( refit, Z.range=NULL, fig.prefix=NULL, fig.name=NULL, n.lgr = 4, ylim=NULL)
{
    if(is.null(Z.range))  Z.range <- c(-1,1);

	pdf.file <- NULL;
	if (!is.null( fig.prefix ))
	{
		pdf.file <- paste( fig.prefix, fig.name, "pdf", sep=".");
		pdf( pdf.file, width=6, height=6);
	}

	n.row <- 4;
	n.col <- 4;
	n.page.seq <- seq(1, NROW(refit), 16);	
	
	for(p in n.page.seq )
	{
		par( mfrow=c( n.row,n.col ) );

		for (i in 1:n.row )
		for (j in 1:n.col )
		{
			if ( (i-1)*4+j + p-1 > NROW(refit) )
				break;

			n.par <- (i-1)*4+j + p-1 ;
			Add.par <- Dom.par <- NULL;

			## In order to draw negative Additive (For genetype=0) instead of positive Additive(For genetype=2)
			if(refit[ n.par,3]>0) Add.par<- refit[ n.par, c(5:8) ] * -1;
			if(refit[ n.par,9]>0) Dom.par<- refit[ n.par, c(11:14) ];
			
			par(mfg=c(i, j));
			draw_single_curve( rownames(refit)[n.par], Z.range, add=Add.par, dom=Dom.par, bLegend=(i==1 & j==1), ylim=ylim );
		}
		
		if (p != n.page.seq[ NROW(n.page.seq)] ) plot.new();
	}
	
	if (!is.null( pdf.file ))
	{
		dev.off();
		cat("* The figure is exported into", pdf.file, "\n");
	}		
}

draw_single_curve<-function( snp_name, Z.range, add=NULL, dom=NULL, bLegend=FALSE, ylim=NULL )
{
    
	old.p1 <- par( mar=c(2,2,1,1)+0.1);
	on.exit( par(old.p1), add = T);

	tp <- seq(-1, 1, 0.05);
	ui <- cbind( rep(1, length(tp)), tp, (3*tp^2-1)/2, (5*tp^3-3*tp)/2 ) ;
	# ui <- cbind( rep(1, length(tp)), tp, (3*tp^2-1)/2, (5*tp^3-3*tp)/2, (35*tp^4-30*tp^2+3)/8 ) ;
	
	y <- c();
	if (!is.null(add)) y <- cbind(y, ui%*%t(add))
	if (!is.null(dom)) y <- cbind(y, ui%*%t(dom))

    if(is.null(ylim))
    {
        ylim.max <- max(y, na.rm=T);
	    ylim.min <- min(y, na.rm=T);
	    ylim.len <- ylim.max - ylim.min
	    ylim.max <- ylim.max + 0.15*abs(ylim.len);
	    ylim.min <- ylim.min - 0.15*abs(ylim.len);
    }
    else
    {
         ylim.min <- ylim[1];
         ylim.max <- ylim[2];
    }
    
	y.num <- length(tp);

	plot( c(0,0), c(0,0), type="n", xaxt="n", yaxt="s", yaxs="i", main=snp_name, cex.main=0.8, cex.axis=0.8, 
		  xlab="Time", ylab="Y", xlim=c(-1,1.2), ylim=c( ylim.min, ylim.max ) );

    axis(1, at=c(-1, -0.5, 0, 0.5, 1), labels=c(Z.range[1],NA, (Z.range[1] + Z.range[2])/2, NA, Z.range[2]), cex.axis=0.8)

	cur.lab <- c();
	cur.col <- c();

	if (!is.null(add))
	{
		y <- ui%*%t(add);
		lines(tp, y, col="orange");
		cur.lab <- c( cur.lab, "Add");
		cur.col <- c( cur.col, "orange");
	}

	if (!is.null(dom))
	{
		y <- ui%*%t(dom);
		lines(tp, y, col="purple");
		cur.lab <- c( cur.lab, "Dom");
		cur.col <- c( cur.col, "purple");
	}

    if (bLegend)
	legend( "topleft", 
			legend = c("Add", "Dom"),
			text.width = strwidth("ABC"),
			text.col = c("orange", "purple"),
			col = c("orange", "purple"),
			lty=1,
			xjust = 1, 
			yjust = 1,
			cex=0.8)
}

draw_refit_CI_curve<-function( refit, Z.range=NULL, fig.prefix=NULL, fig.name=NULL, n.lgr = 4, ReverseCurve = FALSE, q.probs=0.1, ylim=NULL)
{
    if(is.null(Z.range))  Z.range <- c(-1,1);

	pdf.file <- NULL;
	if (!is.null( fig.prefix ))
	{
		pdf.file <- paste( fig.prefix, fig.name, "pdf", sep=".");
		pdf( pdf.file, width=6, height=6);
	}

	n.row <- 4;
	n.col <- 4;
	n.page.seq <- seq(1, NROW(refit), 16);	
	
	for(p in n.page.seq )
	{
		par( mfrow=c( n.row,n.col ) );

		for (i in 1:n.row )
		for (j in 1:n.col )
		{
			if ( (i-1)*4+j + p-1 > NROW(refit) )
				break;

			n.par <- (i-1)*4+j + p-1 ;
		
			par(mfg=c(i, j));
			draw_single_CI_curve( rownames(refit)[n.par], 
			                    Z.range,
			                    fig.name,
			                    hat= refit[ n.par, c(8:11), drop=F ]*ifelse(ReverseCurve, -1, 1) ,
			                    CI0= refit[ n.par, c(13:16), drop=F ]*ifelse(ReverseCurve, -1, 1),
			                    CI1= refit[ n.par, c(18:21), drop=F ]*ifelse(ReverseCurve, -1, 1),
			                    bLegend=(i==1 & j==1),
			                    q.probs=q.probs,
			                    ylim=ylim);
		}
		
		if (p != n.page.seq[ NROW(n.page.seq)] ) plot.new();
	}
	
	if (!is.null( pdf.file ))
	{
		dev.off();
		cat("* The figure is exported into", pdf.file, "\n");
	}		
}
  
get_curve_range <- function( ui, hat, CI0=NULL, CI1=NULL, q.probs=0.1 )
{
    rowMins <- function(tb) { unlist(apply(tb, 1, function(x){min(x, na.rm=TRUE)}) ) }
    rowMaxs <- function(tb) { unlist(apply(tb, 1, function(x){max(x, na.rm=TRUE)}) ) }
    rowQuantile <- function(tb, thres) { unlist(apply(tb, 1, function(x){quantile(x, thres, na.rm=TRUE)}) ) }

	y <- c();
	if (!is.null(hat)) y <- cbind(y, ui %*% t(hat))
		
	if (!is.null(CI0)) 
	{
        y <- cbind(y, ui %*% t(CI0));
	    for (i in 1:NROW(CI0))
	    {
	        hatx <- hat;
	        hatx[1,i] <- CI0[1, i];
	        y <- cbind(y, ui %*% t(hatx))     ;
	    }     
	}
		
	if (!is.null(CI1)) 
	{
        y <- cbind(y, ui %*% t(CI1));
	    for (i in 1:NROW(CI1))
	    {
	        hatx <- hat;
	        hatx[1,i] <- CI1[1, i];
	        y <- cbind(y, ui %*% t(hatx))     ;
	    }     
	}

	if (!is.null(CI0) && !is.null(CI1)) 
    {
        for(i in 1:1000)
        {
        	hatx <- hat;
        	for(k in 1:4)
        	   hatx[1,k] <- runif(1, min(CI0[1, k], CI1[1, k]), max( CI0[1, k], CI1[1, k]) );
        	y <- cbind(y, ui %*% t(hatx))     ;
        }	
    }     
		
	#return( data.frame(y=y[,1], y.min=rowMins(y), y.max=rowMaxs(y))); 
	return( data.frame(y=y[,1], y.min=rowQuantile(y, q.probs), y.max=rowQuantile(y, 1-q.probs))); 
}


draw_single_CI_curve<-function( snp_name, Z.range, effect.lab, hat, CI0=NULL, CI1=NULL, bLegend=FALSE, q.probs=0.1, ylim=NULL )
{
	old.p1 <- par( mar=c(2,2,1,1)+0.1);
	on.exit( par(old.p1), add = T);

	tp <- seq(-1, 1, 0.05);
	ui <- cbind( rep(1, length(tp)), tp, (3*tp^2-1)/2, (5*tp^3-3*tp)/2 ) ;
	# ui <- cbind( rep(1, length(tp)), tp, (3*tp^2-1)/2, (5*tp^3-3*tp)/2, (35*tp^4-30*tp^2+3)/8 ) ;

	y <- get_curve_range( ui, hat, CI0, CI1, q.probs );
	y[,2] <- smooth.spline(x=tp, y=y[,2], df=10)$y;
	y[,3] <- smooth.spline(x=tp, y=y[,3], df=10)$y;

    if(is.null(ylim))
    {
	   ylim.max <- max(y, na.rm=T);
	   ylim.min <- min(y, na.rm=T);
	   ylim.len <- ylim.max - ylim.min
	   ylim.max <- ylim.max + 0.15*abs(ylim.len);
	   ylim.min <- ylim.min - 0.15*abs(ylim.len);
    }
    else
    {
       ylim.min <- ylim[1];
       ylim.max <- ylim[2];
    }
	y.num <- length(tp);

	plot( c(0,0), c(0,0), type="n", xaxt="n", yaxt="s", yaxs="i", main=snp_name, cex.main=0.8, cex.axis=0.8,  
		  xlab="Time", ylab="Y", xlim=c(-1, 1.2) , ylim=c( ylim.min, ylim.max ) );
    
    axis(1, at=c(-1, -0.5, 0, 0.5, 1), labels=c(Z.range[1],NA, (Z.range[1] + Z.range[2])/2, NA, Z.range[2]), cex.axis=0.8)
    
	cur.col <- ifelse(effect.lab=="Add", "orange", "purple");
	 
	polygon(c(tp, rev(tp)), c(y[,2], rev(y[,3])), col=alpha(cur.col, 0.1), border=alpha(cur.col, 0.1)); 
	lines(tp, y[,1], col=cur.col);

    if (bLegend)
		legend( "topleft", 
			legend = effect.lab,
			text.width = strwidth("ABC"),
			text.col = cur.col,
			col = cur.col,
			lty=1,
			xjust = 1, 
			yjust = 1,
			cex=0.8)
}
