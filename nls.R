###############
# Housekeeping
###############
rm(list=ls()) # clear variables

###############
#Simulation
###############
set.seed(20160227)
N <- 1000

P1.x  <-runif(N)
noise <- rnorm(N,mean=1,sd=0.7)
P2.x  <- P1.x + 0.01*noise
P2.x[P2.x<0] <- 0
P2.x[P2.x>1] <- 1
print(head(cbind(P1.x,P2.x)))
print(cor(P1.x,P2.x))
hist(P1.x)
hist(P2.x)

#P2.x <- pbeta(P1.x,1,2)
#print(head(cbind(P1.x,P2.x)))
#print(cor(P1.x,P2.x))
#hist(P1.x)
#hist(P2.x)

###############
# NLS
###############
nls.control(maxiter = 500, tol = 1e-04, minFactor = 1/1024,
            printEval = FALSE, warnOnly = FALSE)

# nls without  starting values for the parameters even if it throw a warning
#nls.1 <- nls(P1.x ~ 1-((1-((P2.x)^(1/b)))^(1/a)), control = list(maxiter = 1000))

a_start = 1
b_start = 1
# nls with starting values for the parameters
nls.2 <-nls(P1.x ~ 1-((1-((P2.x)^(1/b)))^(1/a)),start=list(a=a_start,b=b_start),control = list(maxiter = 10000))

coef <- coef(nls.2)

# Simple goodness of fit
# print(summary(nls.1))
# print(cor(P1.x,predict(nls.1)))
print(summary(nls.2))
print(cor(P1.x,predict(nls.2)))

xseq = seq(100000)/100000
true = xseq
#true = pbeta(xseq,3,1)
aest = coef[1]
best = coef[2]
est  = 1-((1-((xseq)^(1/best)))^(1/aest))

###############
# Diagnositics
###############

# Plot Truth vs Estimate
plot(xseq,true,type='l',ylab="CDF",xlab="x-type",main="True vs. Estimated CDF")
lines(xseq,est,col='red')

# Kolmogorov-Smirnov Test for equality in distributions
print(ks.test(true,est))

