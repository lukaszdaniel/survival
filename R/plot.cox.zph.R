plot.cox.zph <- function(x, resid=TRUE, se=TRUE, df=4, nsmo=40, 
                         var, xlab="Time", ylab="", lty=1:2, col=1, lwd=1,
                         pch= 1, cex=1,
                         hr = FALSE, plot=TRUE, ...) {
    xx <- x$x
    yy <- x$y
    df <- max(df)     # in case df is a vector
    nvar <- ncol(yy)
    pred.x <- seq(from=min(xx), to=max(xx), length=nsmo)
    temp <- c(pred.x, xx)
    lmat <- ns(temp, df=df, intercept=TRUE)
    pmat <- lmat[1:nsmo,]       # for prediction
    xmat <- lmat[-(1:nsmo),]
    if (!is.logical(hr)) stop(gettextf("'%s' argument must be TRUE or FALSE", "hr"))

    if (missing(ylab)) {
        if (hr)  ylab <- paste("HR(t) for", dimnames(yy)[[2]])
            else ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
    }       
    if (missing(var)) var <- 1:nvar
    else {
	if (is.character(var)) var <- match(var, dimnames(yy)[[2]])
	if  (any(is.na(var)) || max(var)>nvar || min(var) <1)
	    stop("Invalid variable requested")
    }
    if (!plot) resid <- FALSE   # only return the curve
    
    # xaxt might be in ...
    if (!is.na(match("xaxt", ...names()))) doxaxis <- FALSE
    else doxaxis <- TRUE

    #
    # Figure out a 'good' set of x-axis labels.  Find 8 equally spaced
    #    values on the 'transformed' axis.  Then adjust until they correspond
    #    to rounded 'true time' values.  Avoid the edges of the x axis, or
    #    approx() may give a missing value
    if (x$transform == 'log') {
	xx <- exp(xx)
	pred.x <- exp(pred.x)
	}
    else if (x$transform != 'identity') {
	xtime <- x$time
        indx <- !duplicated(xx)  #avoid a warning message in R
	apr1  <- approx(xx[indx], xtime[indx], 
                        seq(min(xx), max(xx), length=17)[2*(1:8)])
	temp <- signif(apr1$y,2)
	apr2  <- approx(xtime[indx], xx[indx], temp)
	xaxisval <- apr2$y
	xaxislab <- rep("",8)
	for (i in 1:8) xaxislab[i] <- format(temp[i])
	}
    col <- rep(col, length=2)
    lwd <- rep(lwd, length=2)
    lty <- rep(lty, length=2)

    # Now, finally do the work
    for (i in var) {
        #   Since release 3.1-6, yy can have missing values.  If a covariate is
        # constant within a stratum then it's Shoenfeld residual is identially
        # zero for all observations in that stratum.  These "structural zeros"
        # are marked with an NA.  They contain no information and should not
        # by plotted.  Thus we need to do the spline fit one stratum at a time.
	y <- yy[,i]
        keep <- !is.na(y)
        if (!all(keep)) y <- y[keep]
        
        qmat <- qr(xmat[keep,])
        if (qmat$rank < df) {
            warning(gettextf("spline fit is singular, variable %d skipped", i))
            next
        } 

	yhat <- pmat %*% qr.coef(qmat, y)
	if (resid) yr <-range(yhat, y)
	else       yr <-range(yhat)

	if (se) {
            bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
            xtx <- bk %*% t(bk)
            seval <- ((pmat%*% xtx) *pmat) %*% rep(1, df)
    
	    temp <- 2* sqrt(x$var[i,i]*seval)
	    yup <- yhat + temp
	    ylow<- yhat - temp
	    yr <- range(yr, yup, ylow)
	    }


        if (!plot)  # return the curve
            return(list(x=pred.x, y=cbind(yhat, yup, ylow)))

        if (!hr) {
            if (x$transform=='identity')
                plot(range(xx), yr, type='n', xlab=xlab, ylab=ylab[i], ...)
            else if (x$transform=='log')
                plot(range(xx[keep]), yr, type='n', xlab=xlab, ylab=ylab[i], 
                     log='x', ...)
            else {
                if (doxaxis) {
                    plot(range(xx[keep]), yr, type='n', xlab=xlab, 
                         ylab=ylab[i], xaxt = 'n', ...)
                    axis(1, xaxisval, xaxislab)
                }
                else plot(range(xx[keep]), yr, type='n', xlab=xlab,
                          ylab=ylab[i], ...)
            }
            if (resid) points(xx[keep], y, pch=pch, cex=cex)

            lines(pred.x, yhat, lty=lty[1], col=col[1], lwd=lwd[1])
            if (se) {
                lines(pred.x, yup,  col=col[2], lty=lty[2], lwd=lwd[2])
                lines(pred.x, ylow, col=col[2], lty=lty[2], lwd=lwd[2])
	    }
        } else {    
            if (x$transform=='identity')
                plot(range(xx), exp(yr), type='n', xlab=xlab, ylab=ylab[i],
                     log='y', ...)
            else if (x$transform=='log')
                plot(range(xx[keep]), exp(yr), type='n', xlab=xlab, 
                     ylab=ylab[i], log='xy', ...)
            else {
                if (doxaxis) {
                    plot(range(xx[keep]), exp(yr), type='n', xlab=xlab, 
                         ylab=ylab[i], log='y', xaxt= 'n',...)
                    axis(1, xaxisval, xaxislab)
                }
                else  plot(range(xx[keep]), exp(yr), type='n', xlab=xlab, 
                         ylab=ylab[i], log='y', ...)
	    }
            if (resid) points(xx[keep], exp(y), pch=pch, cex=cex)

            lines(pred.x, exp(yhat), lty=lty[1], col=col[1], lwd=lwd[1])
            if (se) {
                lines(pred.x, exp(yup),  col=col[2], lty=lty[2], lwd=lwd[2])
                lines(pred.x, exp(ylow), col=col[2], lty=lty[2], lwd=lwd[2])
	    }
        }
    }
}
