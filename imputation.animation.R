
source("~/Desktop/Research/scalemixture/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture/scalemix_priors.R")
source("~/Desktop/Research/scalemixture/generic_samplers.R")
source("~/Desktop/Research/scalemixture/scalemix_utils.R")
setwd("~/Desktop/Research/")

library(fields)   # For rdist
# Simulation settings

n.s <- 200        # Number of sites
n.t <- 4          # Number of time points
tau <-9      # Nugget SD
delta <- 0.7      # For R
lambda <- 0.5     # Powered exponential range
gamma <-  1       # Powered exponential smoothness
rho <- 0.1

# Threshold for fitting
thresh <- 11

# Generate fake data

S     <- cbind(seq(0, 1, length=n.s), rep(1, n.s))
# Cor   <- corr.fn(rdist(S), lambda = lambda, gamma = gamma)
Cor   <- corr.fn(rdist(S), rho)
C.Cor <- chol(Cor)

set.seed(3333)
u<-rep(NA,n.t)
R<-rep(NA,n.t)
for(t in 1:n.t){
  u[t] <-runif(1,0,1)
  R[t]<-pow(1/(1-u[t]),delta/(1-delta))
}


X <- matrix(NA, n.s, n.t)
X.s <- matrix(NA, n.s, n.t)
for(t in 1:n.t) {
  Z.t<-crossprod(C.Cor, rnorm(n.s))
  Z.to.W.s<-1/(1-pnorm(Z.t))
  X.s[ ,t] <- R[t]*Z.to.W.s
  X[ ,t] <- X.s[ ,t] + sqrt(tau)*rnorm(n.s)
  plot(X[ ,t], type="l", lwd=3)
  lines(X.s[ ,t], col=2, lwd=2)
  abline(h=thresh, lty=2, lwd=3, col="gray80")
  # readline()
}

true.sigma.s <- R

prob.below <- 0.8
theta.gpd <- c(thresh, 1, 0)


thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau, delta = delta) 
sum(X < thresh.X) / length(X)
cen <- X < thresh.X  

plot(X[ ,1], type="l", lwd=3)
abline(h=thresh.X, lty=2, lwd=3, col="gray80")


## Marginal transformation

library(evd)

Y <- X
# Y[cen] <- scalemix.me.2.gpd(X[cen], theta.mix, theta.gaussian, theta.gpd)
# Y[!cen] <- scalemix.me.2.gpd(X[!cen], theta.mix, theta.gaussian, theta.gpd)
Y[cen] <- NA
Y[!cen] <- scalemix.me.2.gpd(x = X[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below)

# Check the transformation
plot(X[!cen] - gpd.2.scalemix.me(Y[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below))


j <- 1
par(mfrow=c(2, 1), mar=c(1,3,1,1))
plot(X[ ,j], type="l", lwd=3)
abline(h=thresh.X, lty=2, col="gray80", lwd=3)

plot(Y[ ,j], type="l", lwd=3)
abline(h=thresh, lty=2, col="gray80", lwd=3)
par(mfrow=c(1, 1), mai=c(1.02,0.82,0.82,0.42))






################################################################################
n.updates <- 100
thin <- 25

library(scales)

#**********   Eigendecomposition of the correlation matrix   ****************#
eig.Sigma <- eigen(Cor, symmetric=TRUE)
V <- eig.Sigma$vectors
d <- eig.Sigma$values
X.imp <- X
X.star <- X.s

accepted <- rep(0, n.updates)

j <- 1
curr.draw <- X.star[ ,j]
X.s.trace <- matrix(NA, n.s, n.updates)
X.trace <- matrix(NA, n.s, n.updates)

par(mfrow=c(1,1), mar=c(1, 3, 1, 1))
plot(S[ ,1], X.s[ ,j], type="n", lwd=3, xaxt="n", yaxt="n")
abline(h=thresh.X, lty=2, col="gray80", lwd=3)
clip(min(S[ ,1]), max(S[ ,1]), thresh.X, 10)
lines(S[ ,1], X.s[ ,j],  lwd=3)
points(S[ ,1], X[ ,j], col=2, pch=20)
par(mfrow=c(1, 1), mai=c(1.02,0.82,0.82,0.42))



for (i in 1:n.updates) {
  for (k in 1:thin) {
    curr.draw <- X.star[ ,j]
    # Should NOT have asked for X but we are using it here for proposal.
    X.star <- X.s.update.mixture.me.par(R, Y, X.imp, X.star, cen, 
                                    prob.below, theta.gpd, delta,
                                    tau, V, d, v.q=2, n.chain = 100, thresh.X = thresh.X)$X.s
    X.imp[cen] <- update.censored.obs.mixture.me(X.s = X.star, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
  }
  if (max(curr.draw-X.star[ ,j]) > 0.0001) {
    accepted[i] <- 1
    curr.draw <- X.star[ ,j]
  }
  cat("Done with", i, "updates,\n")
  if (i > 50) cat(" Acceptance rate is", mean(accepted[(i-49):i]), "\n")
  X.s.trace[ ,i] <- curr.draw
  X.trace[ ,i] <- X.imp[ ,j]
  
  if (i > 1) {
    points(S[ ,1], X.imp[ ,j], pch=c(20, 4)[cen[ ,j]+1], col=2*(cen[ ,j]+1), cex=0.5)
    dev.off()
  }
  
  png(file=sprintf("./plots/ani_%4.4i-01.png", i), height=5, width=6, units="in", res=150)
  plot(S[ ,1], X[ ,j], type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(S[ ,1], X.s.trace, col=alpha("gray80", 0.05), type="l", add=TRUE)
  points(S[ ,1], X.imp[ ,j], pch=c(20, 4)[cen[ ,j]+1], col=2*(cen[ ,j]+1), cex=0.5)
  abline(h=thresh.X, lty=2, col="gray80", lwd=3)
  dev.off()
  
  png(file=sprintf("./plots/ani_%4.4i-02.png", i), height=5, width=6, units="in", res=150)
  plot(S[ ,1], X[ ,j], type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(S[ ,1], X.s.trace, col=alpha("gray80", 0.05), type="l", add=TRUE)
  points(S[ ,1], X.imp[ ,j], pch=c(20, 4)[cen[ ,j]+1], col=2*(cen[ ,j]+1), cex=0.5)
  lines(S[ ,1], curr.draw, col="gray50", lwd=3)
  abline(h=thresh.X, lty=2, col="gray80", lwd=3)
  dev.off()
  
  
  png(file=sprintf("./plots/ani_%4.4i-03.png", i), height=5, width=6, units="in", res=150)
  plot(S[ ,1], X[ ,j], type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(S[ ,1], X.s.trace, col=alpha("gray80", 0.05), type="l", add=TRUE)
  lines(S[ ,1], curr.draw, col="gray50", lwd=3)
  abline(h=thresh.X, lty=2, col="gray80", lwd=3)
  dev.off()  
  
  
  png(file=sprintf("./plots/ani_%4.4i-04.png", i), height=5, width=6, units="in", res=150)
  plot(S[ ,1], X[ ,j], type="n", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="")
  matplot(S[ ,1], X.s.trace, col=alpha("gray80", 0.05), type="l", add=TRUE)
  lines(S[ ,1], curr.draw, col="gray50", lwd=3)
  abline(h=thresh.X, lty=2, col="gray80", lwd=3)
  
}

dev.off()

par(mfrow=c(1, 1), mai=c(1.02,0.82,0.82,0.42))

save(X.trace, X.s.trace, S, file="ani_files/traces.RData")








