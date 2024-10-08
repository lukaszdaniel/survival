coxph.detail <-  function(object, riskmat=FALSE, rorder=c("data", "time")) {
    method <- object$method
    if (method!='breslow' && method!='efron')
	stop(gettextf("detailed output is not available for the %s method", method))
    rorder <- match.arg(rorder)
    n <- length(object$residuals)
    temp <- coxph.getdata(object, offset=TRUE)
    weights <- temp$weights        #always present if there are weights
    x <- temp$x
    y <- temp$y
    strat <- temp$strata
    Terms <- object$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of object")
    
    nvar <- ncol(x)
    if (ncol(y)==2) {
	mintime <- min(y[,1])
	if (mintime < 0) y <- cbind( 2*mintime -1, y)
	else 	y <- cbind(-1,y)
	}
    if (is.null(strat)) {
	ord <- order(y[,2], -y[,3])
	newstrat <- rep(0,n)
	}
    else {
	ord <- order(strat, y[,2], -y[,3])
	newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
	}
    newstrat[n] <- 1

    # sort the data
    xnew <- x[ord,]
    ynew <- y[ord,]
    storage.mode(y) <- 'double'
    score <- exp(object$linear.predictors)[ord]
    if (is.null(weights)) weights <- rep(1.0, n)
    else                  weights <- weights[ord]

    ndeath <- sum(y[,3])
    if (riskmat) {
	rmat <- integer(ndeath*n)
	}
    else rmat <- as.integer(1)
    
    
    ff <- .C(Ccoxdetail, as.integer(n),
			  as.integer(nvar),
			  ndeath= as.integer(ndeath),
                          center = object$means,
			  y = ynew,
			  as.double(xnew),
			  index = as.integer(newstrat),
			  event2 =as.double(score),
			  weights = as.double(weights),
			  means= c(method=='efron', double(ndeath*nvar-1)),
			  u = double(ndeath*nvar),
			  i = double(ndeath*nvar*nvar),
	                  rmat = rmat,
	                  nrisk2 = double(ndeath),
			  double(nvar*(3 + 2*nvar)))
    keep <- 1:ff$ndeath
    vname<- dimnames(x)[[2]]
    time <- ynew[ff$index[keep],2]
    names(time) <- NULL
    means<- (matrix(ff$means,ndeath, nvar))[keep,]
    score<-  matrix(ff$u, ndeath, nvar)[keep,]
    var <- array(ff$i, c(nvar, nvar, ndeath))[,,keep]
    if (riskmat) {
	rmat <- matrix(0, n, ff$ndeath)
	rmat[,] <- ff$rmat[1:(n*ff$ndeath)]  #in time order
	dimnames(rmat) <- list(NULL, time)
	}

    if (nvar>1) {
	dimnames(means) <- list(time, vname)
	dimnames(score) <- list(time, vname)
	dimnames(var) <- list(vname, vname, time)
	}
    else {
	names(means) <- time
	names(score) <- time
	names(var) <- time
	}

    dimnames(ff$y) <- NULL
    temp <- list(time = time, means=means, nevent=ff$y[keep,1],
	 nrisk = ff$y[keep,2], hazard= ff$y[keep,3], score= score,  imat=var,
	 varhaz=ff$weights[keep], wtrisk = ff$nrisk2[keep])
    if (rorder == "data") {
        temp$y <- y; temp$x <-x}
    else {temp$y <- ynew; temp$x <- xnew}
    if (length(strat)) temp$strata <- table((strat[ord])[ff$index[keep]])
    if (riskmat) {
        if (rorder=="data") {
            temp$riskmat <- matrix(0, nrow(rmat), ncol(rmat),
                                   dimnames= dimnames(rmat))
            temp$riskmat[ord,] <- rmat
        }
        else {
            temp$riskmat <- rmat
            temp$sortorder <- ord
        }
    }
    if (!all(weights==1)) {
	temp$weights <- weights
	temp$nevent.wt <- ff$event2[keep]
	temp$nrisk.wt  <- ff$nrisk2[keep]
	}
    temp
    }
