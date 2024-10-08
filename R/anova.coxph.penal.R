# The first section of this is identical to anova.coxph
anova.coxph.penal <- function (object, ...,  test = 'Chisq') {
    if (!inherits(object, "coxph"))
        stop("argument must be a cox model")

    # All the ... args need to be coxph or coxme fits.  If any of them
    #  have a name attached, e.g., 'charlie=T' we assume a priori
    #  that they are illegal
    #
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
	           rep(FALSE, length(dotargs))
             else (names(dotargs) != "")
    if (any(named)) 
        warning(gettextf("The following arguments to anova.coxph(..) are invalid and dropped: %s",
            paste(deparse(dotargs[named]), collapse = ", ")))
    dotargs <- dotargs[!named]

    # frailties don't work, due to mutability.  
     has.frailty <- function(x) any(x$pterms==2)

    if (length(dotargs) >0) {
        # Check that they are all cox or coxme models
        is.coxmodel <-unlist(lapply(dotargs, function(x) inherits(x, "coxph")))
        is.coxme <- unlist(lapply(dotargs, function(x) inherits(x, "coxme")))
        if (!all(is.coxmodel | is.coxme))
            stop("All arguments must be Cox models")
        
        if (any(sapply(dotargs, has.frailty)))
            stop("anova command does not handle frailty terms")

        if (any(is.coxme)) {
            # We need the anova.coxmelist function from coxme
            # If coxme is not loaded the line below returns NULL
            temp <- getS3method("anova", "coxmelist", optional=TRUE)
            if (is.null(temp)) 
                stop("a coxme model was found and library coxme is not loaded")
            else return(temp(c(list(object), dotargs), test = test))
        }
        else return(anova.coxphlist(c(list(object), dotargs), test = test))
    }

    #
    # The argument is a single Cox model 
    # Show the nested list of models generated by this model.
    # By tradition the sequence is main effects (in the order found in
    #  the model statement), then 2 way interactions, then 3, etc.
    #  One does this by using the "assign" attribute of the model matrix.
    #  This does not work for penalized terms, however, so we use a mixed
    #  strategy.  The penalized terms do not participate in interactions
    #  (which are the terms for which assign is really handy).  Use
    #  the model frame for the penalized terms, and assign for all the
    #  others.
    if (length(object$rscore)>0)
        stop("Can't do anova tables with robust variances")
    if (has.frailty(object))
        stop("anova command does not handle frailty terms")

    has.strata <- !is.null(attr(terms(object), "specials")$strata)
    # The following line causes pspline terms to be re-evaluated correctly
    #  The predvars attr for them does not retrieve the correct penalty
    attr(object$terms, "predvars") <- NULL
    mf <- stats::model.frame(object)  # we must have the model frame
    Y <- model.response(mf)
    X <- model.matrix(object, mf)
    assign <- attr(X, 'assign') 
    if (has.strata) {
        stemp <- untangle.specials(terms(object), "strata")
        if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
        else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
        strats <- as.numeric(strata.keep)
        }

    pname <- names(object$pterms)[object$pterms >0]
    pindex <- match(pname, attr(terms(object), "term.labels"))

    alevels <- sort(unique(assign)) #if there are strata the sequence has holes
    nmodel <- length(alevels)

    df <- integer(nmodel+1)  #this will hold the df vector
    loglik <- double(nmodel+1) #and the loglike vector
    df[nmodel+1] <- if (is.null(object$df)) sum(!is.na(object$coefficients))
                    else sum(object$df)
    loglik[nmodel+1] <- object$loglik[2]
    df[1] <- 0
    loglik[1] <- object$loglik[1]
    
    # Now refit intermediate models
    assign2 <- assign[!(assign %in% pindex)]
    pform   <- paste("mf[['", pname, "']]", sep='')
    for (i in seq.int(1, length.out=nmodel-1)){
        j <- assign2[assign2 <= alevels[i]]
        if (length(j)) form <- "Y ~ X[,j]"
        else           form <- "Y ~"
        
        form <- paste(c(form, pform[pindex <= i]), collapse=" +") 
        if (length(object$offset)) 
            form <- paste(form, " + offset(object$offset")
        if (has.strata) form <- paste(form, " + strata(strats)")
                                      
        tfit <- coxph(as.formula(form))
        df[i+1] <- if (!is.null(tfit$df)) sum(tfit$df)
                       else sum(!is.na(tfit$coefficients))
        loglik[i+1] <- tfit$loglik[2]
    }
                              
    table <- data.frame(loglik=loglik, Chisq=c(NA, 2*diff(loglik)), 
                        Df=c(NA, diff(df))) 

    if (nmodel == 0) #failsafe for a NULL model
             table <- table[1, , drop = FALSE]

    if (length(test) >0 && test[1]=='Chisq') {
        table[['Pr(>|Chi|)']] <- pchisq(table$Chisq, table$Df, lower.tail=FALSE)
        }
    row.names(table) <- c('NULL', attr(terms(object), "term.labels"))

    title <- paste("Analysis of Deviance Table\n Cox model: response is ",
		   deparse(object$terms[[2]]),
		   "\nTerms added sequentially (first to last)\n", 
		   sep = "")
    structure(table, heading = title, class = c("anova", "data.frame"))
}
