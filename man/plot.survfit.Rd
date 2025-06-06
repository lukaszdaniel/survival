\name{plot.survfit}
\alias{plot.survfit}
\title{
  Plot method for \code{survfit} objects 
}
\usage{
\method{plot}{survfit}(x, conf.int=, mark.time=FALSE, 
 pch=3, col=1, lty=1, lwd=1, cex=1, log=FALSE, xscale=1, yscale=1,  
 xlim, ylim, xmax, fun, 
 xlab="", ylab="", xaxs="r",  conf.times, conf.cap=.005,
 conf.offset=.012,
 conf.type = c("log", "log-log", "plain", "logit", "arcsin", "none"),
 mark, noplot="(s0)", cumhaz=FALSE,
 firstx, ymin, cumprob=FALSE, \dots)
}
\arguments{
  \item{x}{
    an object of class \code{survfit}, usually returned by the
    \code{survfit} function. 
  }
  \item{conf.int}{if TRUE/FALSE,
    determines whether pointwise confidence intervals will be plotted.
    The default is to 
    do so if there is only 1 curve, i.e., no strata, using 95\% confidence
    intervals
    Alternatively, this can be a numeric value giving the desired
    confidence level.
  }
  \item{mark.time}{
    controls the labeling of the curves.  If set to \code{FALSE}, no
    labeling is done. 
    If \code{TRUE}, then curves are marked at each censoring time.
    If \code{mark} is a
    numeric vector then curves are marked at the specified time points. 
  }
  \item{pch}{
    vector of characters which will be used to label the curves. 
    The \code{points} help file contains examples of the possible marks.
    A single string such as "abcd" is treated as a vector
    \code{c("a", "b", "c", "d")}.
    The vector is reused cyclically if it is shorter than the number of
    curves.  If it is present this implies \code{mark.time = TRUE}.
  }
  \item{col}{
    a vector of integers specifying colors for each curve. 
    The default value is 1. 
  }
  \item{lty}{
    a vector of integers specifying line types for each curve. 
    The default value is 1. 
  }
  \item{lwd}{
    a vector of numeric values for line widths. The default value is 1. 
  }
  \item{cex}{
    a numeric value specifying the size of the marks. 
    This is not treated as a vector; all marks have the same size.
  }
  \item{log}{
    a logical value, if TRUE the y axis wll be on a log scale. 
    Alternately, one of the standard character strings "x", "y", or "xy"
    can be given to specific logarithmic horizontal and/or vertical axes.
  }
  \item{xscale}{
    a numeric value used like \code{yscale} for labels on the x axis.  
    A value of 365.25 will give labels in years instead of the original days.  
  }
  \item{yscale}{
    a numeric value used to multiply the labels on the y axis. 
    A value of 100, for instance, would be used to give a percent scale. 
    Only the labels are 
    changed, not the actual plot coordinates, so that adding a curve with 
    "\code{lines(surv.exp(...))}", say,  
    will perform as it did without the \code{yscale} argument. 
  }
  \item{xlim,ylim}{optional limits for the plotting region.
  }
  \item{xmax}{
    the maximum horizontal plot coordinate.  This can be used to shrink
    the range of a plot.  It shortens the curve before plotting it, so
    that unlike using the \code{xlim} graphical parameter, warning
    messages about out of bounds points are not generated.
  }
  \item{fun}{
    an arbitrary function defining a transformation of the survival
    (or probability in state, or cumulative hazard) curves. 
    For example \code{fun=log} is an alternative way to draw a log-survival curve 
    (but with the axis labeled with log(S) values), 
    and \code{fun=sqrt} would generate a curve on square root scale. 
    Four often used transformations can be specified with a character 
    argument instead: \code{"S"} gives the usual survival curve,
    \code{"log"} is the same as using the \code{log=T} option, 
    \code{"event"} or \code{"F"} plots the empirical CDF \eqn{F(t)= 1-S(t)}
    (f(y) = 1-y), and
    \code{"cloglog"} creates a complimentary log-log survival plot (f(y) = 
    log(-log(y)) along with log scale for the x-axis).
    The terms \code{"identity"} and \code{"surv"} are
    allowed as synonyms for \code{type="S"}.
    The argument \code{"cumhaz"} causes the cumulative hazard function
    to be plotted.
  }
  \item{xlab}{
    label given to the x-axis. 
  }
  \item{ylab}{
    label given to the y-axis. 
  }
  \item{xaxs}{
    either \code{"S"} for a survival curve or a standard x axis style as
    listed in \code{par}; "r" (regular) is the R default.
    Survival curves have historically been displayed with the curve
    touching the y-axis,
    but not touching the bounding box of the plot on the other 3 sides,
    Type \code{"S"} accomplishes this by manipulating the plot range and
    then using the \code{"i"} style internally.
    The "S" style is becoming increasingly less common, however.
  }
  \item{conf.times}{optional vector of times at which to place a
    confidence bar on the curve(s).  If present, these will be used
    instead of confidence bands.}
  \item{conf.cap}{width of the horizontal cap on top of the confidence
    bars; only used if conf.times is used.  A value of 1 is the width of
    the plot region.}
  \item{conf.offset}{the offset for confidence bars, when there are
    multiple curves on the plot.  A value of 1 is the width of the plot
    region. If this is a single number then each curve's bars are offset
    by this amount from the prior curve's bars, if it is a vector the values are
    used directly.}
  \item{conf.type}{
       One of \code{"plain"}, \code{"log"} (the default),
	\code{"log-log"}, \code{"logit"}, or \code{"arcsin"} or \code{"none"}.
	Only enough of the string to uniquely identify it is necessary.
	The last option causes confidence intervals not to be
	generated.  The plain option gives intervals
	\code{curve +- k *se(curve)}, where k is determined from
	\code{conf.int}.  The log option calculates intervals based on the
	cumulative hazard or log(survival). The log-log option bases the
	intervals on the log hazard or log(-log(survival)), and the
	logit option on log(survival/(1-survival)) The default is
	"plain" for a plot of the cumulative hazard. 
      }
  \item{mark}{a historical alias for \code{pch}}
  \item{noplot}{for multi-state models, curves with this label will not
    be plotted.  (Also see the \code{istate0} argument in
    \code{survcheck}.)}
  \item{cumhaz}{plot the cumulative hazard rather than the probability
    in state or survival.  Optionally, this can be a numeric vector
    specifying which columns of the \code{cumhaz} component to plot.}
  \item{ymin}{ this will normally be given as part of the \code{ylim}
    argument}
  \item{firstx}{this will normally be given as part of the \code{xlim}
    argument.}
  \item{cumprob}{for a multi-state curve, plot the probabilities in
    state 1, (state1 + state2), (state1 + state2 + state3), \ldots.
    If \code{cumprob} is an integer vector the totals will be in the
    order indicated.}
  \item{\dots}{other arguments that will be passed forward to the
    underlying plot method, such as xlab or ylab.}
}
\value{
a list with components \code{x} and \code{y}, containing the coordinates of the last point 
on each of the curves (but not the confidence limits).   
This may be useful for labeling.
If \code{cumprob=TRUE} then \code{y} will be a matrix with one row per
curve and \code{x} will be all the time points.  This may be useful for
adding shading.
}
\description{
A plot of survival curves is produced, one curve for each strata. 
The \code{log=T} option does extra work to avoid log(0), and to try to create a 
pleasing result.  If there are zeros, they are plotted by default at 
0.8 times the smallest non-zero value on the curve(s).

Curves are plotted in the same order as they are listed by \code{print}
(which gives a 1 line summary of each).
This will be the order in which \code{col}, \code{lty}, etc are used.
}
\details{
  If the object contains a cumulative hazard curve, then
  \code{fun='cumhaz'} will plot that curve, otherwise it will plot
  -log(S) as an approximation.  Theoretically, S =
  \eqn{\exp(-\Lambda)}{exp(-H)} where S is the survival and
  \eqn{\Lambda}{H} is the cumulative hazard.  The same relationship
  holds for estimates of S and \eqn{\Lambda}{H} only in special cases,
  but the approximation is often close.

  When the \code{survfit} function creates a multi-state survival curve
  the resulting object also has class `survfitms'.
  Competing risk curves are a common case.
  In this situation the \code{fun} argument is ignored.

  When the \code{conf.times} argument is used, the confidence bars are
  offset by \code{conf.offset} units to avoid overlap.
  The bar on each curve are the confidence interval for the time point
  at which the bar is drawn, i.e., different time points for each curve.
  If curves are steep at that point, the visual impact can sometimes
  substantially differ for positive and negative values of
  \code{conf.offset}.
}

\note{In prior versions the behavior of \code{xscale} and
  \code{yscale} differed: the first changed the scale both for the plot
  and for all subsequent actions such as adding a legend, whereas \code{yscale}
  affected only the axis label.  This was normalized in version 2-36.4,
  and both parameters now only affect the labeling.

  In versions prior to approximately 2.36 a \code{survfit} object did
  not contain the cumulative hazard as a separate result, and the use of
  fun="cumhaz" would plot the approximation -log(surv) to the cumulative
  hazard.  When cumulative hazards were added to the object, the
  \code{cumhaz=TRUE} argument to the plotting function was added.
  In version 2.3-8 the use of fun="cumhaz" became a synonym for
  \code{cumhaz=TRUE}.
  }
\seealso{ \code{\link{points.survfit}},
\code{\link{lines.survfit}},   
\code{\link{par}},  
\code{\link{survfit}}  
}
\examples{
leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml) 
plot(leukemia.surv, lty = 2:3) 
legend(100, .9, c("Maintenance", "No Maintenance"), lty = 2:3) 
title("Kaplan-Meier Curves\nfor AML Maintenance Study") 
lsurv2 <- survfit(Surv(time, status) ~ x, aml, type='fleming') 
plot(lsurv2, lty=2:3, fun="cumhaz", 
	xlab="Months", ylab="Cumulative Hazard") 
}
\keyword{survival}
\keyword{hplot}
