# The underlying method for survfit calls
survfit <- function(formula, ...) {
    UseMethod("survfit")
}

# The primary use is survfit(formula,...) which is below, or
#  survfit(coxph-fit-object)
#
# The formula, data, weights, subset, and na.action args are similar to other
#  model functions.  See the help page for the remaining args, and see the
#  methods document for computational details.
survfit.formula <- function(formula, data, weights, subset, 
                            na.action, stype=1, ctype=1, 
                            id, cluster, robust, istate, 
                            timefix=TRUE, etype, model=FALSE, error, 
                            entry=FALSE,  time0=FALSE, ...) {

    Call <- match.call()
    if (missing(formula)) stop(gettextf("'%s' argument is required", "formula"))
    newform <- removeDoubleColonSurv(formula)
    if (!is.null(newform)) {
        formula <- newform$formula
        if (newform$newcall) Call$formula <- formula
    }
    Call[[1]] <- as.name('survfit')  #make nicer printout for the user

    # create a copy of the call that has only the arguments we want,
    #  and use it to call model.frame()
    indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                    'istate', 'id', 'cluster', "etype"), names(Call), nomatch=0)
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(stats::model.frame)
    # make sure we evaluate Surv, strata etc in the survival namespace
    newform <- removeDoubleColonSurv(formula)
    if (!is.null(newform)) {
        temp$formula <- newform$formula
        if (newform$newcall) Call$formula <- formula
    }
    mf <- eval.parent(temp)

    Terms <- terms(formula, c("strata", "cluster"))
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
            stop("Interaction terms are not valid for this function")

    n <- nrow(mf)
    Y <- model.response(mf)
    # The Surv2 approach is experimental, and may not survive.  If a data set is
    #  in that form, covert it to a "standard" layout
    if (inherits(Y, "Surv2")) {
        # this is Surv2 style data
        # if there are any obs removed due to missing, remake the model frame
        if (length(attr(mf, "na.action"))) {
            temp$na.action <- na.pass
            mf <- eval.parent(temp)
        }
        if (!is.null(attr(Terms, "specials")$cluster))
            stop("cluster() cannot appear in the model statement")
        new <- surv2data(mf)
        mf <- new$mf
        istate <- new$istate
        id <- new$id
        Y <- new$y
        if (anyNA(mf[-1])) { #ignore the response variable still found there
            if (missing(na.action)) temp <- get(getOption("na.action"))(mf[-1])
            else temp <- na.action(mf[-1])
            omit <- attr(temp, "na.action")
            mf <- mf[-omit,]
            Y <- Y[-omit]
            id <- id[-omit]
            istate <- istate[-omit]
        }                      
        n <- nrow(mf)
    }       
    else {
        if (!is.Surv(Y)) stop("response must be a survival object")
        id <- model.extract(mf, "id")
        istate <- model.extract(mf, "istate")
    }
    if (n==0) stop("data set has no non-missing observations")

    casewt <- model.extract(mf, "weights")
    if (is.null(casewt)) casewt <- rep(1.0, n)
    else {
        if (!is.numeric(casewt)) stop(gettextf("'%s' must be numeric", "weights"))
        if (any(!is.finite(casewt))) stop("weights must be finite") 
        if (any(casewt <0)) stop("weights must be non-negative")
        casewt <- as.numeric(casewt)  # transform integer to numeric
    }

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    cluster <- model.extract(mf, "cluster")
    temp <- untangle.specials(Terms, "cluster")
    if (length(temp$vars)>0) {
        if (length(cluster) >0) stop("cluster appears as both an argument and a model term")
        if (length(temp$vars) > 1) stop("can not have two cluster terms")
        cluster <- mf[[temp$vars]]
        Terms <- Terms[-temp$terms]
    }

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(mf[ll])

    # Backwards support for the now-depreciated etype argument
    etype <- model.extract(mf, "etype")
    if (!is.null(etype)) {
        if (attr(Y, "type") == "mcounting" ||
            attr(Y, "type") == "mright")
            stop("cannot use both the etype argument and mstate survival type")
        if (length(istate)) 
            stop("cannot use both the etype and istate arguments")
        status <- Y[,ncol(Y)]
        etype <- as.factor(etype)
        temp <- table(etype, status==0)

        if (all(rowSums(temp==0) ==1)) {
            # The user had a unique level of etype for the censors
            newlev <- levels(etype)[order(-temp[,2])] #censors first
        }
        else newlev <- c(" ", levels(etype)[temp[,1] >0])
        status <- factor(ifelse(status==0,0, as.numeric(etype)),
                             labels=newlev)

        if (attr(Y, 'type') == "right")
            Y <- Surv(Y[,1], status, type="mstate")
        else if (attr(Y, "type") == "counting")
            Y <- Surv(Y[,1], Y[,2], status, type="mstate")
        else stop("etype argument incompatable with survival type")
    }
                         
    # Deal with the near-ties problem
    if (!is.logical(timefix) || length(timefix) > 1)
        stop(gettextf("invalid '%s' value", "timefix"))
    if (timefix) newY <- aeqSurv(Y) else newY <- Y
    
    if (missing(robust)) robust <- NULL
    if (!is.logical(time0)) stop(gettextf("'%s' argument must be TRUE or FALSE", "time0"))

    # Call the appropriate helper function, which does the real work
    # Turnbull for interval censored data
    # KM = Kaplan-Meier for survival with a single endpoint
    # AJ = Aalen-Johansen for multi-state models
    if (attr(Y, 'type') == 'left' || attr(Y, 'type') == 'interval')
        temp <-  survfitTurnbull(X, newY, casewt, cluster= cluster,
                                 robust= robust, time0= time0,...)
    else if (attr(Y, 'type') == "right" || attr(Y, 'type')== "counting")
        temp <- survfitKM(X, newY, casewt, stype=stype, ctype=ctype, id=id, 
                          cluster=cluster, robust=robust, entry=entry, 
                          time0=time0, ...)
    else if (attr(Y, 'type') == "mright" || attr(Y, "type")== "mcounting")
        temp <- survfitAJ(X, newY, weights=casewt, stype=stype, ctype=ctype, 
                          id=id, cluster=cluster, robust=robust, 
                          istate=istate, entry=entry, time0=time0, ...)
    else {
        # This should never happen
        stop("unrecognized survival type")
    }

    # If a stratum had no one beyond start.time, the length 0 gives downstream
    #  failure, e.g., there is no sensible printout for summary(fit, time= 100)
    #  for such a curve with no one who makes it to 100.
    temp$strata <- temp$strata[temp$strata >0]  
    if (is.null(temp$states)) class(temp) <- 'survfit'
    else class(temp) <- c("survfitms", "survfit")

    if (!is.null(attr(mf, 'na.action')))
            temp$na.action <- attr(mf, 'na.action')
    if (model) temp$model <- mf
    temp$call <- Call
    temp
    }

# Return the dimension of a survfit object
# Number of groups = right hand side of formula in survfit.formula, or strata in
#  survival from a Cox model
# Number of rows in the newdata data set, for a Cox model
# Number of states, for a multi-state outcome
dim.survfit <- function(x) {
    d1name <- "strata"
    d2name <- "data"
    d3name <- "states"    
    if (is.null(x$strata))  {d1 <- d1name <- NULL} else d1 <- length(x$strata)
    # d3 is present for a survfitms object, null otherwise
    if (is.null(x$states))  {
        d3 <- d3name <- NULL
        if (is.matrix(x$surv)) d2 <- ncol(x$surv)
        else {d2 <- d2name <- NULL}
    } else {
        d3 <- length(x$states) 
        dp <- dim(x$pstate)
        if (length(dp) ==3) d2 <- dp[2]
        else {d2 <- d2name <- NULL}
    }
    
    dd <- c(d1, d2, d3)
    names(dd) <- c(d1name, d2name, d3name)
    dd
}

# The subscript function for a survfit object
# there is a separate subscript function for survfitms objects
# See survfit:subscript in the methods document
"[.survfit" <- function(x, ... , drop=TRUE) {
    nmatch <- function(indx, target) { 
        # This function lets R worry about character, negative, or 
        #  logical subscripts.
        #  It always returns a set of positive integer indices
        temp <- 1:length(target)
        names(temp) <- target
        temp[indx]
    }
    
    if (!inherits(x, "survfit")) stop("[.survfit called on non-survfit object")
    ndots <- ...length()      # the simplest, but not avail in R 3.4
    # ndots <- length(list(...))# fails if any are missing, e.g. fit[,2]
    # ndots <- if (missing(drop)) nargs()-1 else nargs()-2  # a workaround

    dd <- dim(x)
    # for dd=NULL, an object with only one curve, x[1] is always legal
    if (is.null(dd)) dd <- c(strata=1L) # survfit object with only one curve
    dtype <- match(names(dd), c("strata", "data", "states"))

    if (ndots >0 && !missing(..1)) i <- ..1 else i <- NULL
    if (ndots> 1 && !missing(..2)) j <- ..2 else j <- NULL
    
    if (ndots > length(dd)) 
        stop("incorrect number of dimensions")
    if (length(dtype) > 2) stop("invalid survfit object")  # should never happen
    if (is.null(i) && is.null(j)) {
        # called with no subscripts given -- return x untouched
        return(x)
    }
    
    # Code below is easier if "i" is always the strata
    if (dtype[1] !=1) {
        dtype <- c(1, dtype)
        j <- i; i <- NULL
        dd <- c(1, dd)
        ndots <- ndots +1
    }       

   # We need to make a new one
    newx <- vector("list", length(x))
    names(newx) <- names(x)
    for (k in c("logse", "version", "conf.int", "conf.type", "type", "call"))
        if (!is.null(x[[k]])) newx[[k]] <- x[[k]]
    class(newx) <- class(x)
    
    if (ndots== 1 && length(dd)==2) {
        # one subscript given for a two dimensional object
        # If one of the dimensions is 1, it is easier for me to fill in i and j
        if (dd[1]==1) {j <- i; i<- 1}
        else if (dd[2]==1) j <- 1
        else {
            #  the user has a mix of rows/cols
            index <- 1:prod(dd)
            itemp <- matrix(index, nrow=dd[1])
            keep <- itemp[i]   # illegal subscripts will generate an error
            if (length(keep) == length(index) && all(keep==index)) return(x)

            ii <- row(itemp)[keep]
            jj <- col(itemp)[keep]
            # at this point we have a matrix subscript of (ii, jj)
            # expand into a long pair of rows and cols
            temp <- split(seq(along.with=x$time), 
                          rep(1:length(x$strata), x$strata))
            indx1 <- unlist(temp[ii])   # rows of the surv object
            indx2 <- rep(jj, x$strata[ii])
        
            # return with each curve as a separate strata
            newx$n <- x$n[ii]
            for (k in c("time", "n.risk", "n.event", "n.censor", "n.enter"))
                if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[indx1]
            k <- cbind(indx1, indx2)
            for (j in c("surv", "std.err", "upper", "lower", "cumhaz",
                        "std.chaz", "influence.surv", "influence.chaz"))
                if (!is.null(x[[j]])) newx[[j]] <- (x[[j]])[k]
            temp <- x$strata[ii]
            names(temp) <- 1:length(ii)
            newx$strata <- temp
            return(newx)
        }
    }
    
    # irow will be the rows that need to be taken
    #  j the columns (of present)
    if (is.null(x$strata)) {
           if (is.null(i) || all(i==1)) irow <- seq(along.with=x$time)
           else stop("subscript out of bounds")
           newx$n <- x$n
    }
    else { 
        if (is.null(i)) indx <- seq(along.with= x$strata)
        else indx <- nmatch(i, names(x$strata)) #strata to keep
        if (any(is.na(indx))) 
            stop(gettextf("strata %s not matched", paste(i[is.na(indx)], collapse = " ")))
        # Now, indx may not be in order: some can use curve[3:2] to reorder
        #  The list/unlist construct will reorder the data
        temp <- split(seq(along.with =x$time), 
                      rep(1:length(x$strata), x$strata))
        irow <- unlist(temp[indx])
        
        if (length(indx) <=1 && drop) newx$strata <- NULL
        else               newx$strata  <- x$strata[i]

        newx$n <- x$n[indx]
        if (length(indx) ==1 & drop) x$strata <- NULL
        else    newx$strata <- x$strata[indx]
    }

    if (!is.matrix(x[["surv"]])) {  # no j dimension
        for (k in c("time", "n.risk", "n.event", "n.censor", "n.enter",
               "surv", "std.err", "cumhaz", "std.chaz", "upper", "lower",
               "influence.surv", "influence.chaz"))
            if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[irow]
    }
       
    else { # 2 dimensional object
        if (is.null(j)) j <- seq.int(ncol(x$surv))
        # If the curve has been selected by strata and keep has only
        #  one row, we don't want to lose the second subscript too
        if (length(irow)==1)  drop <- FALSE

        for (k in c("time", "n.risk", "n.event", "n.censor", "n.enter"))
                 if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[irow]
        for (k in c("surv", "std.err", "cumhaz", "std.chaz", "upper", "lower",
               "influence.surv", "influence.chaz"))
            if (!is.null(x[[k]])) newx[[k]] <- (x[[k]])[irow, j, drop=drop]
        # for a survfit.coxph object, newdata is a data frame whose rows match j
        if (!is.null(x[["newdata"]])) newx[["newdata"]] <- x[["newdata"]][j,]
    }
    newx
}

# Once upon a time I allowed survfit to be called without the ~1 portion of
# the formula. This was a mistake for multiple reasons, and I later removed it.
# The method below helps give a useful error message in some cases.
survfit.Surv <- function(formula, ...)
    stop("the survfit function requires a formula as its first argument")

# The confidence interval compututation is needed in more than one place, so
# make it a function.  To do: make this user callable?

survfit_confint <- function(p, se, logse=TRUE, conf.type, conf.int,
                            selow, ulimit=TRUE) {
    zval <- qnorm(1- (1-conf.int)/2, 0,1)
    if (missing(selow)) scale <- 1.0
    else scale <- ifelse(selow==0, 1.0, selow/se)  # avoid 0/0 at the origin
    if (!logse) se <- ifelse(se==0, 0, se/p)   # se of log(survival) = log(p)

    if (conf.type=='plain') {
        se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
        if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
        else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
    }
    else if (conf.type=='log') {
        #avoid some "log(0)" messages
        xx <- ifelse(p==0, NA, p)  
        se2 <- zval* se 
        temp1 <- exp(log(xx) - se2*scale)
        temp2 <- exp(log(xx) + se2)
        if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
        else  list(lower= temp1, upper= temp2)
    }
    else if (conf.type=='log-log') {
        xx <- ifelse(p==0 | p==1, NA, p)
        se2 <- zval * se/log(xx)
        temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
        temp2 <- exp(-exp(log(-log(xx)) + se2))
        list(lower = temp1 , upper = temp2)
    }
    else if (conf.type=='logit') {
        xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
        se2 <- zval * se *(1 + xx/(1-xx))
 
        temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
        temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
        list(lower = temp1, upper=temp2)
    }
    else if (conf.type=="arcsin") {
        xx <- ifelse(p==0, NA, p)
        se2 <- .5 *zval*se * sqrt(xx/(1-xx))
        list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
             upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
    }
    else stop("invalid conf.int type")
}
