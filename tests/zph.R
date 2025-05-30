library(survival)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

test1 <- data.frame(time=  c(9, 3,1,1,6,6,8),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0))

# Verify that cox.zph computes a score test
# First for the Breslow estimate
r <- (3 + sqrt(33))/2   # actual MLE for log(beta)
U <- c(1/(r+1), 3/(r+3), -r/(r+3), 0)   # score statistic
imat <- c(r/(r+1)^2, 3*r/(r+3)^2, 3*r/(r+3)^2, 0)  # information matrix
g = c(1, 6, 6, 9)  # death times

u2 <- c(sum(U), sum(g*U))  # first derivative
i2 <- matrix(c(sum(imat), sum(g*imat), sum(g*imat), sum(g^2*imat)),
               2,2)  # second derivative
sctest <- solve(i2, u2) %*% u2

# Verify that centering makes no difference for the test (though i2 changes)
g2 <- g - mean(g)
u2b <- c(sum(U), sum(g2*U))
i2b <- matrix(c(sum(imat), sum(g2*imat), sum(g2*imat), sum(g2^2*imat)),
               2,2) 
sctest2 <- solve(i2b, u2b) %*% u2b
all.equal(sctest, sctest2)

# Now check the program
fit1 <- coxph(Surv(time, status) ~ x, test1, ties='breslow')
aeq(fit1$coef, log(r))
zp1 <- cox.zph(fit1, transform='identity', global=FALSE)
aeq(zp1$table[,1], sctest)
aeq(zp1$y, resid(fit1, type="scaledsch"))

dummy <- rep(0, nrow(test1))
fit1b <- coxph(Surv(dummy, time, status) ~ x, test1, ties='breslow')
aeq(fit1b$coef, log(r))
zp1b <- cox.zph(fit1b, transform='identity', global=FALSE)
aeq(zp1b$table[,1], sctest)
# the pair of tied times gets reversed in the zph result
#  but since the 'y' values are only used to plot it doesn't matter
aeq(zp1b$y[c(1,3,2,4)], resid(fit1b, type="scaledsch"))  

# log time check
g3 <- log(g) - mean(log(g))
u3 <- c(sum(U), sum(g3*U))  # first derivative
i3 <- matrix(c(sum(imat), sum(g3*imat), sum(g3*imat), sum(g3^2*imat)),
               2,2)  # second derivative
sctest3 <- solve(i3, u3) %*% u3
zp3 <- cox.zph(fit1, transform='log', global=FALSE)
aeq(zp3$table[,1], sctest3)

# Efron approximation
phi <- acos((45/23)*sqrt(3/23))
r <- 2*sqrt(23/3)* cos(phi/3)   # actual MLE for log(beta)
U <- c(1/(r+1), 3/(r+3), -r/(r+5), 0)   # score statistic
imat <- c(r/(r+1)^2, 3*r/(r+3)^2, 5*r/(r+5)^2, 0)  # information matrix

u4 <- c(sum(U), sum(g3*U))  # first derivative
i4 <- matrix(c(sum(imat), sum(g3*imat), sum(g3*imat), sum(g3^2*imat)),
               2,2)  # second derivative
sctest4 <- solve(i4, u4) %*% u4

fit4 <- coxph(Surv(time, status) ~ x, test1, ties='efron')
aeq(fit4$coef, log(r))
zp4 <- cox.zph(fit4, transform='log', global=FALSE)
aeq(zp4$table[,1], sctest4)
aeq(zp4$y, resid(fit4, type="scaledsch"))

fit5 <-  coxph(Surv(dummy, time, status) ~ x, test1, ties="efron")
aeq(fit5$coef, log(r))
zp5 <- cox.zph(fit5, transform="log", global=FALSE)
aeq(zp5$table[,1], sctest4)

# Artificial stratification
test2 <- rbind(test1, test1)
test2$group <- rep(letters[1:2], each=nrow(test1))
# U, imat, and sctest will all double
dummy <- c(dummy, dummy)
fit6 <- coxph(Surv(dummy, time, status) ~ x + strata(group), test2)
aeq(fit6$coef, log(r))
zp6 <- cox.zph(fit6, transform="log", globa=FALSE)
aeq(zp6$table[,1], 2*sctest4)

# A multi-state check, 2 covariates
#  Verify that the multi-state result = the single state Cox models
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))
table(event)

ct1 <- coxph(Surv(etime, event) ~ sex + age, mgus2, id=id, ties='efron')
ct2 <- coxph(Surv(etime, event=='pcm') ~ sex + age, mgus2)
ct3 <- coxph(Surv(etime, event=='death') ~ sex + age, mgus2)

zp1 <- cox.zph(ct1, transform='identity')
zp2 <- cox.zph(ct2, transform='identity')
zp3 <- cox.zph(ct3, transform='identity')
aeq(zp1$table[1:2,], zp2$table[1:2,])
aeq(zp1$table[3:4,], zp3$table[1:2,])

# Now add a starting time of zero
dummy <- rep(0, nrow(mgus2))
ct4 <- coxph(Surv(dummy, etime, event) ~ sex + age, mgus2, id=id, ties='efron')
ct5 <- coxph(Surv(dummy, etime, event=='pcm') ~ sex + age, mgus2)
ct6 <- coxph(Surv(dummy, etime, event=='death') ~ sex + age, mgus2)
zp4 <- cox.zph(ct4, transform='identity')
zp5 <- cox.zph(ct5, transform='identity')
zp6 <- cox.zph(ct6, transform='identity')
aeq(zp4$table[1:2,], zp5$table[1:2,])
aeq(zp4$table[3:4,], zp6$table[1:2,])


# Direct check of a multivariate model with start, stop data
p1 <- pbcseq[!duplicated(pbcseq$id), 1:6]
pdata <- tmerge(p1[, c("id", "trt", "age", "sex")], p1, id=id,
                death = event(futime, status==2))
pdata <- tmerge(pdata, pbcseq, id=id, bili=tdc(day, bili),
                edema = tdc(day, edema), albumin=tdc(day, albumin),
                protime = tdc(day, protime))
pfit <- coxph(Surv(tstart, tstop, death) ~ log(bili) + albumin + edema +
                  age + log(protime), data = pdata, ties='efron')
zp7  <- cox.zph(pfit, transform='log', global=FALSE)

direct <- function(fit) {
    nvar <- length(fit$coef)
    dt <- coxph.detail(fit)
    gtime <- log(dt$time) - mean(log(dt$time))
    # key idea: at any event time I have a first deriviative vector
    #     c(dt$score[i,],  gtime[i]* dt$score[i,]) 
    # and second derivative matrix
    #     dt$imat[,,i]              gtime[i]   * dt$imat[,,i]
    #     gtime[i]*dt$imat[,,i]     gtime[i]^2 * dt$imat[,,i]
    # for the expanded model, where imat[,,i] is symmetric, 
    # and colSums(dt$score) =0  (since the model converged)
    # 
    # Create score tests for adding one time-dependent variable 
    #  gtime * x[,j] at a time: first derivative of this test is
    #     c(dt$score[i,],  gtime[i]* dt$score[i,j]) 
    #  and etc.
    unew  <-   colSums(gtime * dt$score)
    temp1 <- apply(dt$imat, 1:2, sum)
    temp2 <- apply(dt$imat, 1:2, function(x) sum(x*gtime))
    temp3 <- apply(dt$imat, 1:2, function(x) sum(x * gtime^2))

    score <- double(nvar)
    smat  <- matrix(0., nvar+1, nvar+1)  # second deriv matrix for the test
    smat[1:nvar, 1:nvar] <- temp1
    for (i in 1:nvar) {
        smat[nvar+1,] <- c(temp2[i,], temp3[i,i])
        smat[,nvar+1] <- c(temp2[,i], temp3[i,i])
        utemp <- c(rep(0,nvar), unew[i])
        score[i] <- solve(smat, utemp) %*% utemp
    }  
    list(sctest = score, u= c(colSums(dt$score), unew), 
         imat=cbind(rbind(temp1, temp2), rbind(temp2, temp3)))
}   

aeq(zp7$table[,1], direct(pfit)$sctest)

# Last, make sure that NA coefficients are ignored
d1 <- survSplit(Surv(time, status) ~ ., veteran, cut=150, episode="epoch")
fit <- coxph(Surv(tstart, time, status) ~ celltype:strata(epoch) + age, d1)
zz <- cox.zph(fit)

fit2 <- coxph(Surv(tstart, time, status) ~ celltype:strata(epoch) + age, d1,
              x=TRUE)
zz2 <- cox.zph(fit2)

x2 <- fit2$x[, !is.na(fit$coefficients)][,-1]
fit3 <- coxph(Surv(tstart, time, status) ~ age + x2, d1)
all.equal(fit3$loglik, fit2$loglik)
zz3 <- cox.zph(fit3)

all.equal(unclass(zz)[1:7], unclass(zz2)[1:7]) #ignore the call component
all.equal(as.vector(zz$table), as.vector(zz3$table)) # variable names change
