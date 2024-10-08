#
# A tiny multi-state example
#
library(survival)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
mtest <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5, 5),
                    t1= c(0, 4, 9,  0,  2,  0, 2, 8,  1, 3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    st= c(1, 2,  1, 2,  3,  1, 3, 0,  2,  0))

mtest$state <- factor(mtest$st, 0:3, c("censor", "a", "b", "c"))

if (FALSE) {
    # this graph is very useful when debugging
    temp <- survcheck(Surv(t1, t2, state) ~1, mtest, id=id)
    plot(c(0,11), c(1,5.1), type='n', xlab="Time", ylab= "Subject")
    with(mtest, segments(t1+.1, id, t2, id, col=as.numeric(temp$istate)))
    event <- subset(mtest, state!='censor')
    text(event$t2, event$id+.2, as.character(event$state))
}

mtest <- mtest[c(1,3,2,4,5,7,6,10, 9, 8),]  #not in time order

mfit <- survfit(Surv(t1, t2, state) ~ 1, mtest, id=id, time0=FALSE)

# True results
#
#time       state                    probabilities
#         entry  a   b  c         entry  a    b     c
#
#0        124                      1     0    0     0
#1+       1245
#2+       1235   4                3/4   1/4   0     0    4 -> a, add 3
#3+       123    4   5            9/16  1/4  3/16   0    5 -> b
#4+        23    14  5            6/16  7/16 3/16   0    1 -> a
#5+        3     14  5            3/16  7/16 6/16   0    2 -> b, exits
#8+        3     1   5  4         3/16  7/32 6/16  7/32  4 -> c
#9+                  15            0     0  19/32 13/32  1->b, 3->c & exit
# 10+            1   5                19/64 19/64 13/32  1->a

aeq(mfit$n.risk,  matrix(c(4,4,3,2,1,1,0,0,
                                 0,1,1,2,2,1,0,0,
                                 0,0,1,1,1,1,2,1,
                                 0,0,0,0,0,1,0,0), ncol=4))
aeq(mfit$pstate,  matrix(c(24, 18, 12,  6, 6, 0, 0,  0,
                                  8,  8, 14, 14, 7, 0,  9.5, 9.5, 
                                  0,  6,  6, 12, 12,19,9.5, 9.5, 
                                  0,  0,  0,  0, 7, 13, 13, 13)/32, ncol=4))
aeq(mfit$n.transition, matrix(c(1,0,1,0,0,0,0,0,
                           0,0,0,0,0,0,1,0,
                           0,1,0,1,0,0,0,0,
                           0,0,0,0,0,1,0,0,
                           0,0,0,0,0,1,0,0,
                           0,0,0,0,1,0,0,0), ncol=6))
all.equal(mfit$time, c(2, 3, 4, 5, 8, 9, 10, 11))

# Somewhat more complex.
#  Scramble the input data
#  Not everyone starts at the same time or in the same state
#  Case weights
#
tdata <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5,  5),
                    t1= c(0, 4, 9,  1,  2,  0, 2, 8,  1,  3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    st= c(1, 2,  1, 2,  3,  1, 3, 0,  3,  0),
                    i0= c(4, 1,  2, 1,  4,  4, 1, 3,  2,  3),
                    wt= 1:10)

tdata$st <- factor(tdata$st, c(0:3),
                    labels=c("censor", "a", "b", "c"))
tdata$i0 <- factor(tdata$i0, c(4, 1:3),
                    labels=c("entry", "a", "b", "c"))
if (FALSE) {
    #useful picture     
    temp <- survcheck(Surv(t1, t2, st) ~1, tdata, id=id, istate=i0)
    plot(c(0,11), c(1,5.5), type='n', xlab="Time", ylab= "Subject")
    with(tdata, segments(t1+.1, id, t2, id, col=as.numeric(temp$istate)))
    with(subset(tdata, st!= "censor"),
          text(t2, id+.15, as.character(st)))
    with(tdata, text((t1+t2)/2, id+.25, wt))
    with(subset(tdata, !duplicated(id)),
           text(t1, id+.15, as.character(i0)))
    #abline(v=c(2:5, 8:11), lty=3, col='gray')
}

tfun <- function(data=tdata) {
    reorder <- c(10, 9, 1, 2, 5, 4, 3, 7, 8, 6)
    new <- data[reorder,]
    new
}

# These weight vectors are in the order of tdata
# w[9] is the weight for subject 5 at time 1.5, for instance
# p0 is defined as all those at risk just before the first event, which in
#  this data set is entry:a at time 2 for id=4; id 1,2,4,5 at risk

# When the functions below were written, the entry state was listed last.
# Currently the entry state is first, so "[swap]" was added to the aj routines
#  rather than rearranging the formulas
swap <- c(4,1,2,3)
p0 <- function(w) c( w[1]+ w[6], w[4], w[9], 0)/ (w[1]+ w[4] + w[6] + w[9])

#  aj2 = Aalen-Johansen H matrix at time 2, etc.
aj2 <- function(w) {
    #subject 4 moves from entry to 'a'  
    rbind(c(1, 0, 0, 0),  
          c(0, 1, 0, 0),
          c(0, 0, 1, 0),
          c(w[6], 0, 0, w[1])/(w[1] + w[6]))[swap, swap] 
}
aj3 <- function(w) rbind(c(1, 0, 0, 0),   
                         c(0, 0, 1, 0),  # 5 moves from b to c
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))[swap,swap]
aj4 <- function(w) {
    # subject 1 moves from entry to a
    rbind(c(1, 0, 0, 0),
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(w[1], 0, 0, w[5])/(w[1] + w[5])) [swap, swap]
}
aj5 <- function(w) {
    # subject 2 from a to b
    rbind(c(w[2]+w[7], w[4], 0, 0)/(w[2]+ w[4] + w[7]), 
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))[swap, swap]
}
aj8 <- function(w) rbind(c(w[2], 0, w[7], 0)/(w[2]+ w[7]), # 4  a to c
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))[swap, swap]
aj9 <- function(w) rbind(c(0, 1, 0, 0), # 1  a to b
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(0, 0, 1 ,0)) [swap, swap] # 3 entry to c
aj10 <- function(w)rbind(c(1, 0, 0, 0),
                         c(1, 0, 0, 0),  #1 b to a
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))[swap, swap]

#time       state               
#         a   b  c  entry
#
#1        2   5     14       initial distribution
#2        24  5     1        4 -> a, add 3
#3        24     5  13       5 from b to c
#4       124     5   3       1 -> a
#5        14     5   3       2 -> b, exits
#8        1      45  3       4 -> c
#9            1  45          1->b, 3->c & exit
#10       1      45          1->a

# P is a product of matrices
dopstate <- function(w) {
    p1 <- p0(w)
    p2 <- p1 %*% aj2(w)
    p3 <- p2 %*% aj3(w)
    p4 <- p3 %*% aj4(w)
    p5 <- p4 %*% aj5(w)
    p8 <- p5 %*% aj8(w)
    p9 <- p8 %*% aj9(w)
    p10<- p9 %*% aj10(w)
    rbind(p2, p3, p4, p5, p8, p9, p10, p10)
}

# Check the pstate estimate
w1 <- rep(1,10)
mtest2 <- tfun(tdata)  # scrambled order
mfit2 <- survfit(Surv(t1, t2, st) ~ 1, tdata, id=id, istate=i0, 
                 time0=FALSE) # ordered
aeq(mfit2$pstate, dopstate(w1))
aeq(mfit2$p0, p0(w1))

mfit2b <- survfit(Surv(t1, t2, st) ~ 1, mtest2, id=id, istate=i0, time0=FALSE)
aeq(mfit2b$pstate, dopstate(w1))
aeq(mfit2b$p0, p0(w1))

mfit2b$call <- mfit2$call <- NULL
all.equal(mfit2b, mfit2) 
aeq(mfit2$transitions, c(2,0,1,0, 0,2,0,0, 1,1,1,0, 0,0,0,2))

# Now the harder one, where subjects change weights
mfit3  <- survfit(Surv(t1, t2, st) ~ 1, tdata, id=id, istate=i0,
                  weights=wt, influence=TRUE, time0=FALSE)
aeq(mfit3$p0, p0(1:10))
aeq(mfit3$pstate, dopstate(1:10))
    

# The derivative of a matrix product AB is (dA)B + A(dB) where dA is the
#  elementwise derivative of A and etc for B.
# dp0 creates the derivatives of p0 with respect to each subject, a 5 by 4
#  matrix
# All the functions below are hand coded for a weight vector that is in
#  exactly the same order as the rows of mtest.
# Since p0 = (w[1]+ w[6], w[4], w[9], 0)/ (w[1]+ w[4] + w[6] + w[9])
#      and subject id is 1,1,1, 2, 3, 4,4,4, 5,5
#   we get the derivative below
# 

dp0 <- function(w) {   # influence just before the first event
  p <- p0(w)
  w0 <- w[c(1,4,6,9)]  # the 4 obs at the start, subjects 1, 2, 4, 5
  rbind(c(1,0, 0, 0) - p,   # subject 1 affects p[entry]
        c(0,1, 0, 0) - p,   # subject 2 affects p[a]
        0,                   # subject 3 affects none
        c(1, 0, 0, 0) - p,   # subject 4 affect p[entry]
        c(0, 0, 1, 0) - p)/   # subject 5 affects p[b]
      sum(w0)
}
  

dp2 <- function(w) {
    h2 <- aj2(w)   # H matrix at time 2
    part1 <- dp0(w) %*% h2

    # 1 and 4 in entry, obs 4 moves from entry to a
    mult  <- p0(w)[1]/(w[1] + w[6])  #p(t-) / weights in state
    part2 <- rbind((c(1,0,0,0)- h2[1,]) * mult,
                   0,
                   0,
                   (c(0,1,0,0) - h2[1,]) * mult,
                   0)
    part1 + part2
}

dp3 <- function(w) {
    dp2(w) %*% aj3(w)
}

dp4 <- function(w) {
    h4 <- aj4(w)   # H matrix at time 4
    part1 <- dp3(w) %*% h4

    # subjects 1 and 3 in state entry (obs 1 and 5) 1 moves to a
    mult <- dopstate(w)[2,1]/ (w[1] + w[5])   # p_1(time 4-0) / wt
    part2 <- rbind((c(0,1,0,0)- h4[1,]) * mult,
                   0,
                   (c(1,0,0,0)- h4[1,]) * mult,
                   0,
                   0)
    part1 + part2
}
dp5 <- function(w) {
    h5 <- aj5(w)   # H matrix at time 5
    part1 <- dp4(w) %*% h5

    # subjects 124 in state a (obs 2,4,7), 2 goes to b
    mult <- dopstate(w)[3,2]/ (w[2] + w[4] + w[7]) 
    part2 <- rbind((c(0,1,0,0)- h5[2,]) * mult,
                   (c(0,0,1,0)- h5[2,]) * mult,
                   0,
                   (c(0,1,0,0)- h5[2,]) * mult,
                   0)
    part1 + part2
}
dp8 <- function(w) {
    h8 <- aj8(w)   # H matrix at time 8
    part1 <- dp5(w) %*% h8

    # subjects 14 in state a (obs 2 &7), 4 goes to c
    mult <- dopstate(w)[4, 2]/ (w[2] + w[7]) 
    part2 <- rbind((c(0,1,0,0)- h8[2,]) * mult,
                   0,
                   0,
                   (c(0,0,0,1)- h8[2,]) * mult,
                   0)
    part1 + part2
}
dp9 <- function(w) dp8(w) %*% aj9(w)
dp10<- function(w) dp9(w) %*% aj10(w)

#
# Feb 4 2024: discovered that the variance computation above is incorrect.
# Let U = influence for phat, with one row per observation in the data
# The weighted per subject influence is Z(t)= BDU(t) where 
#  B= rbind(c(1,1,1,0,0,0,0,0,0,0), 
#           c(0,0,0,1,0,0,0,0,0,0),
#           c(0,0,0,0,1,0,0,0,0,0),
#           c(0,0,0,0,0,1,1,1,0,0),
#           c(0,0,0,0,0,0,0,0,1,1))
#  and D= diag(1:10)
# which can be summarized as "weight each row, then add over subjects". 
# The variance at time t is the column sums of Z^2(t) (elementwise squares)
#
# The code above for dp0, dp2, etc returns BU, which matches the computation of
#  the influence in survfitci.c. If the weight for a given subject is constant 
#  over time, then BDU= WBU where W is the diagonal matrix of per-subject 
#  weights: survfitci.c implicitly made this assumption, and was correct
#  in this case. It returned U as the influence, which matches dp0 etc.
#
# survfitci.c has been replaced by survfitaj.c, which uses the careful
# derivations in the methods vignette, and returns BDU.
# The checks below have been changed to a case with constant weights per
#  subject.  R code to test for general weights is in mstate2.R
#
w1 <- tdata$id 
mfit4 <- survfit(Surv(t1, t2, st) ~1, tdata, id=id, weights=id, istate=i0,
                 influence=TRUE, time0= FALSE)
aeq(mfit4$influence[,1,], 1:5*dp2(w1))  #time 2
aeq(mfit4$influence[,2,], 1:5*dp3(w1))
aeq(mfit4$influence[,3,], 1:5*dp4(w1))
aeq(mfit4$influence[,4,], 1:5*dp5(w1))
aeq(mfit4$influence[,5,], 1:5*dp8(w1)) # time 8
aeq(mfit4$influence[,6,], 1:5* dp9(w1))
aeq(mfit4$influence[,7,], 1:5* dp10(w1))
aeq(mfit4$influence[,8,], 1:5* dp10(w1)) # no changes at time 11

ssq <- function(x) sqrt(sum(x^2))
temp2 <- apply(mfit4$influence.pstate, 2:3, ssq)
aeq(temp2, mfit4$std.err)  

if (FALSE) { # old test, survfitci returned the time 0 influence as well
    w1 <- 1:10
    aeq(mfit3$influence[,1,], dp0(w1))
    aeq(mfit3$influence[,2,], dp2(w1))
    aeq(mfit3$influence[,3,], dp3(w1))
    aeq(mfit3$influence[,4,], dp4(w1))
    aeq(mfit3$influence[,5,], dp5(w1))
    aeq(mfit3$influence[,6,], dp8(w1))
    aeq(mfit3$influence[,7,], dp9(w1))
    aeq(mfit3$influence[,8,], dp10(w1))
    aeq(mfit3$influence[,9,], dp10(w1)) # no changes at time 11
} # end of if (FALSE)

# The cumulative hazard at each time point is remapped from a matrix
#  into a vector (in survfit)
# First check out the names
nstate <- length(mfit4$states)
temp <- matrix(0, nstate, nstate)
indx1 <- match(rownames(mfit4$transitions), mfit4$states)
indx2 <- match(colnames(mfit4$transitions), mfit4$states, nomatch=0)
temp[indx1, indx2] <- mfit4$transitions[, indx2>0]
# temp is an nstate by nstate version of the transitions matrix
from <- row(temp)[temp>0]
to   <- col(temp)[temp>0]
all.equal(colnames(mfit4$cumhaz), paste(from, to, sep=':'))

# check the cumulative hazard
temp <- mfit4$n.risk[,from]
hazard <- mfit4$n.transition/ifelse(temp==0, 1, temp)
aeq(apply(hazard, 2, cumsum), mfit4$cumhaz)

