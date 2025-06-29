#
# Print out information about a rate table: it's dimensions and keywords
#
summary.ratetable <- function(object, ...) {
    rtable<-object
    if (!inherits(rtable, 'ratetable')) gettextf("'%s' argument is not an object of class %s", "object", dQuote("ratetable"))

    att <- attributes(rtable)
    ncat <- length(dim(rtable))
    cat (" Rate table with", ncat, "dimensions:\n")
    if (is.null(att$dimid)) dimid <- names(dimnames(rtable))
    else dimid <- att$dimid
    for (i in 1:ncat) {
        # One of 'factor' (old style table) or "type" (new style) should exist
        if (!is.null(att$factor)) {
            if (att$factor[i]==0) {
                cat("\t", dimid[i], " ranges from ", 
                    format(min(att$cutpoints[[i]])), " to ", 
                    format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                    " categories\n", sep='')
                }
            else if(att$factor[i]==1) {
                cat("\t", dimid[i], " has levels of: ",
                    paste(att$dimnames[[i]], collapse=' '), "\n", sep='')
                }
            else {
                cat("\t", dimid[i], " ranges from " , 
                    format(min(att$cutpoints[[i]])), " to ", 
                    format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                    " categories,\n\t\tlinearly interpolated in ",
                    att$factor[i], " steps per division\n", sep='')
                }
            }
        else {
            if (att$type[i]==1) {
                cat("\t", dimid[i], " has levels of: ",
                    paste(att$dimnames[[i]], collapse=' '), "\n", sep='')
                }
            else if (att$type[i]>2) { #date
                if (is.numeric(att$cutpoints[[i]])) { #old, numeric
                    # This format is > 5 years out of date
                    #  but some user might keep an old rate table around
                    cat("\t", dimid[i], " ranges from " , 
                        format(as.Date(min(att$cutpoints[[i]]),
                                       origin='1960/01/01')), " to ", 
                        format(as.Date(max(att$cutpoints[[i]]),
                                       origin='1960/01/01')),
                        "; with ", att$dim[i],
                        " categories\n", sep='')
                    }
                else # newer, Date
                    cat("\t", dimid[i], " ranges from " , 
                        format(min(att$cutpoints[[i]])), " to ", 
                        format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                        " categories\n", sep='')
                }

            else {
                cat("\t", dimid[i], " ranges from ", 
                    format(min(att$cutpoints[[i]])), " to ", 
                    format(max(att$cutpoints[[i]])), "; with ", att$dim[i],
                    " categories\n", sep='')
                }
            }
        }
            
    invisible(att)
    }

