library(doParallel)
library(foreach)




source("~/Desktop/Research/scalemixture_tproc/scalemix_utils.R")
source("~/Desktop/Research/scalemixture_tproc/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture_tproc/scalemix_priors.R")
source("~/Desktop/Research/scalemixture_tproc/generic_samplers.R")
source("~/Desktop/Research/scalemixture_tproc/scalemix_sampler_02.R")

load("./ffwi.RData")

n.cores <- 20
cl <- makePSOCKcluster(n.cores)
# library(snow);library(Rmpi)
# cl <- parallel::makeCluster(n.cores, type="MPI")

# clusterEvalQ(cl, Rcpp::sourceCpp('~/Desktop/Research/scalemixture_reg/integration.cpp'))
clusterEvalQ(cl, source("~/Desktop/Research/scalemixture_tproc/scalemix_likelihoods.R"))
clusterEvalQ(cl, source("~/Desktop/Research/scalemixture_tproc/scalemix_utils.R"))
clusterEvalQ(cl, source("~/Desktop/Research/scalemixture_tproc/generic_samplers.R"))
clusterEvalQ(cl, source("~/Desktop/Research/scalemixture_tproc/scalemix_priors.R"))

registerDoParallel(cl)

library(fields)   # For rdist
set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")


# -------------------------------------------------
# ------------ 1. Simulation settings -------------
# -------------------------------------------------

n.s <- nrow(ffwi.locs)            # Number of sites
n.t <- ncol(ffwi.mat.max)         # Number of time points

rho <- 1; nu <- 3/2
Cor   <- corr.fn(rdist(ffwi.locs), rho, nu)
C.Cor <- chol(Cor)



# ------------------------------------------
# -------------- 2. Y to X.s -----------------
# ------------------------------------------

Y <- ffwi.mat.max
lon <- matrix(rep(ffwi.locs$lon,n.t),ncol=n.t)

# (a) thresh and theta.gpd
require(gPdtest)
prob.below <- 0.98
thresh <- quantile(Y, prob.below)
cen <- Y < thresh; mean(cen)
gpd <- gpd.fit(Y[Y>=thresh]-thresh, method='amle')
theta.gpd <- c(thresh, a, b) #scale, shape
#all(a+b*lon[!cen] > 0)


# (b) Pick initial values for delta and tau ==> MARGINAL BEHAVIOR of X
#     Note that pmixture.me fails at low quantile
tau = 5; gamma = 6; sigma =2
thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau, gamma = gamma, sigma = sigma)
num_samples <- 200000
Y_sample <- rep(NA,num_samples)
for(j in 1:num_samples){
  R_tmp <- sigma*sqrt(gamma*invgamma::rinvgamma(1, shape=gamma/2, scale=2))
  Y_sample[j] <- R_tmp*rnorm(1)+rnorm(1, sd=sqrt(tau))
}
X.cdf.emp <- ecdf(Y_sample)
# plot(X.cdf.emp, xlim=c(-2,quantile(Y_sample,0.80)), main="CDF of X")
# points(100, pmixture_me(100,tau,delta), pch=20, col='red')
# points(30, pmixture_me(30,tau,delta), pch=20, col='red')



# (c) Marginal transformation from Y to X
X <- matrix(NA, n.s, n.t)
U <- matrix(NA, n.s, n.t)

## ----- Censored part
for(i in 1:n.s){
  tmp_cdf <-ecdf(Y[i,])
  unifs <- tmp_cdf(Y[i, cen[i,]]); U[i, cen[i,]] <- unifs
}
x.vals <- seq(-50, thresh.X, length=1000)  # design points for pmixture.me
cdf.vals <- X.cdf.emp(x.vals); x.vals[which(cdf.vals>0)[1]]
X[cen] <- qmixture.me.interp(U[cen], tau_sqd = tau, gamma = gamma, sigma = sigma, x.vals = x.vals, cdf.vals = cdf.vals)
# qmixture.me.interp(min(unifs), tau_sqd = tau, delta = delta, x.vals = x.vals, cdf.vals = cdf.vals) close to 1
# X[, which(colSums(Y==min(Y))==1)]

## ----- Uncensored part
require(evd)
Scale <- a+b*lon[!cen]
unifs <- (1-prob.below) * pgpd(Y[!cen], loc=theta.gpd[1], scale=Scale, shape=shape) + prob.below; U[!cen]<-unifs
X[!cen] <- qmixture.me.interp(unifs, tau_sqd = tau, gamma = gamma, sigma = sigma, n.x=500)

# which(colSums(X<1)>0)
# col <- 3
# plot(Y[,col],type='l',ylab="HOT")
# abline(h=thresh, lty=2,col="red")
# plot(ffwi.locs, pch=15, col=gray((Y[,col]-min(Y[,col]))/(max(Y[,col])-min(Y[,col]))), main="HOT", xlab='x1', ylab='x2')
# plot(ffwi.locs[cen[,col],],pch=15, main="HOT", xlab='x1', ylab='x2')
# 
# plot(X[,col], type='p', ylab='X')
# abline(h=thresh.X, lty=2,col="red")



# (d) Get good initial values for R
require(mvtnorm)
Tr_inverse<-function(w){return(qnorm(1-(1/w)))}
lik<-function(params,X.s){
  R<-params[1]
  n<-length(X.s)
  return(dmvnorm(X.s/R,mean=rep(0,n),sigma=Cor,log=TRUE)-n*log(R))
}

R<-rep(NA,ncol(X))
for(i in 1:n.t){
  X.s.tmp<-X[,i]
  
  temp<-tryCatch(optim(par=max(min(X.s.tmp)/2,1),lik,X.s=X.s.tmp,method="L-BFGS-B",lower=0,
                       upper=Inf,control = list(pgtol=10e-20,trace=TRUE))$par, error= function(e){
                         Ri <- seq(-min(X.s.tmp), min(X.s.tmp), length=100)
                         Lik<-rep(NA, 100)
                         for(i in 1:100){
                           Lik[i]<-lik(Ri[i],X.s.tmp)
                         }
                         cat("Optim failed\n")
                         return(Ri[which.min(Lik)])
                       })
  
  R[i]<-temp
}

# (e) Get good initial values for X.s
# require('doParallel'); require('foreach'); registerDoParallel(cores=20)
load("./ffwi_org.RData")
Y.org <- ffwi.Mat  # Original Y with NA's denoting missing observation
missing <- is.na(Y.org)

v.q<-matrix(2.4^2,nrow=n.s, ncol=n.t)
eig.Sigma <- eigen(Cor, symmetric=TRUE)
V <- eig.Sigma$vectors
d <- eig.Sigma$values
Res<-foreach(t = 1:n.t, .combine = "cbind", .inorder=FALSE, .noexport = "update_X_s_onetime")%dopar%{
  X.s.tmp<-X[,t]
  
  res<- update_X_s_onetime(Y[,t], X[,t], X.s.tmp, cen[,t], prob.below, theta.gpd, shape, lon[,t], gamma, sigma, 
                           tau, thresh.X, v.q[,t], R[t], V=V, d=d, missing = missing)
  c(res$X.s, res$accept)
}
Accepted<-Res[(n.s+1):(2*n.s),]
X.s.first<-Res[1:n.s,]

for(i in 1:10){
  thin = 10
  Accepted<-matrix(0, nrow=n.s, ncol=n.t)
  for(iter in 1:thin){
    Res<-foreach(t = 1:n.t, .combine = "cbind", .inorder=FALSE, .noexport = "update_X_s_onetime")%dopar%{
      
      res<- update_X_s_onetime(Y[,t], X[,t], X.s.first[,t], cen[,t], prob.below, theta.gpd, shape, lon[,t], gamma, sigma, 
                               tau, thresh.X, v.q[,t], R[t], V, d, missing = missing)
      c(res$X.s, res$accept)
    }
    Accepted<-Accepted + Res[(n.s+1):(2*n.s),]
    X.s.first<-Res[1:n.s,]
  }
  r.hat.X.s <- Accepted/thin
  c.0 <- 10; c.1 <- 0.8 ; k <- 3 ; gamma1 <- c.0 / (i + k)^(c.1)
  v.q <- exp(log(v.q) + gamma1*(r.hat.X.s - 0.41))
}


# col <- 11  #21, 11
# par(mfrow=c(2,2))
# plot(W[,col],type='l',ylab="W")
# plot(Y.s[,col],type='l',ylab="HOT")
# points(Y[,col], pch=20)
# abline(h=thresh, lty=2,col="red")
# plot(X[,col], type='p', ylab='X')
# abline(h=thresh.X, lty=2,col="red")
# lines(X.s.first[,col], col='red',lwd=2)
# par(mfrow=c(1,1))

X.s <- X.s.first
Y_tmp <- Y.org; Y_tmp[cen] <- NA
Y <- Y_tmp




# ----------------------------------------------------------
## --------------- 4. Running Metropolis -------------------
# ----------------------------------------------------------

initial.values <- list(gamma = gamma, sigma = sigma, rho=rho, nu=nu, tau=tau, theta.gpd=theta.gpd, shape=shape, prob.below=prob.below, X.s=X.s, R=R)
# true.params <- initial.values
n.updates <- 50000
thin <- 10
echo.interval <- 50

save(Y, X, Y.org, thresh.X, cen, missing, initial.values, file="Initial.RData")

# Calculate
Res <- adaptive.metr(z = R, starting.theta = theta.gpd[2:3],
                     likelihood.fn = theta.gpd.update.mixture.me.likelihood, prior.fn = unif.a.b,
                     hyper.params = max(lon[!cen]), n.updates = 30000, prop.Sigma = NULL, adapt.cov = TRUE, return.prop.Sigma.trace = FALSE,
                     r.opt = .234, c.0 = 10, c.1 = .8, K = 10,
                     shape=shape, lon=lon, Y =Y, X.s = X.s, cen = cen, prob.below = prob.below,gamma=gamma, sigma=sigma, tau_sqd = tau, loc = thresh, 
                     thresh.X=thresh.X, missing = missing)
prop.Sigma.theta <- cov(Res$trace[10000:30000,])
sd.ratio <- sqrt(prop.Sigma.theta[1,1]/prop.Sigma.theta[2,2])
prop.Sigma <- list(gpd.corr=cor(Res$trace[10000:30000,])[1,2], theta.gpd=prop.Sigma.theta/prop.Sigma.theta[1,1])

Res<-adaptive.metr(z = R, starting.theta = c(rho, nu),
            likelihood.fn = rho.update.mixture.me.likelihood, prior.fn = rho.prior,
            hyper.params = 1, n.updates = 30000, prop.Sigma = NULL, adapt.cov=TRUE, return.prop.Sigma.trace = FALSE,
            r.opt = .234, c.0 = 10, c.1 = .8, K = 10,
            X.s = X.s, R = R, S = ffwi.locs)
prop.Sigma.rho <- cov(Res$trace[10000:30000,])
sd.ratio.rho <- sqrt(prop.Sigma.rho[1,1]/prop.Sigma.rho[2,2])
prop.Sigma <- c(prop.Sigma, list(rho.corr=cor(Res$trace[10000:30000,])[1,2], rho=prop.Sigma.rho/prop.Sigma.rho[1,1]))

Res<-adaptive.metr(z = R, starting.theta = c(gamma, sigma), 
                   likelihood.fn = delta.update.mixture.me.likelihood, prior.fn = interval.unif,
                   hyper.params = 1, n.updates = 30000, prop.Sigma = NULL, adapt.cov=TRUE, return.prop.Sigma.trace = FALSE,
                   r.opt = .234, c.0 = 10, c.1 = .8, K = 10,
                   Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, theta.gpd = theta.gpd, shape = shape, lon = lon, tau_sqd = tau,
                   missing = missing)
prop.Sigma.delta <- cov(Res$trace[10000:30000,])
sd.ratio.delta <- sqrt(prop.Sigma.delta[1,1]/prop.Sigma.delta[2,2])
prop.Sigma <- c(prop.Sigma, list(delta.corr=cor(Res$trace[10000:30000,])[1,2], delta=prop.Sigma.delta/prop.Sigma.delta[1,1]))


sigma.m<-list(theta.gpd=(2.4/2)^2)
cat("sd.ratio=",sd.ratio,"\n")
cat("corr=",prop.Sigma$gpd.corr,"\n")
cat("sigma11=",prop.Sigma.theta[1,1],"\n")

if(prop.Sigma.theta[2,2]>0  & sd.ratio<100) {
  scalemix.sampler.02.DA(Y=Y, S=ffwi.locs, cen=cen, lon=lon, thresh=thresh,
                         initial.values=initial.values,
                         n.updates=n.updates, thin=thin,
                         experiment.name="Huser-wadsworth-DA",
                         echo.interval=echo.interval,
                         sigma.m=sigma.m, prop.Sigma=prop.Sigma, missing = missing,
                         sd.ratio=sd.ratio, sd.ratio.rho = sd.ratio.rho,sd.ratio.delta = sd.ratio.delta, lower.prob.lim=0.5,
                         holdout=TRUE, S_full = ffwi.locs_full, n.holdout = 5)
}

stopCluster(cl)
# scalemix.sampler.03(Y=Y, S=S, cen=cen, thresh=thresh,
#                     initial.values=initial.values,
#                     n.updates=n.updates, thin=thin,
#                     experiment.name="Huser-wadsworth-scaleshape",
#                     echo.interval=echo.interval,
#                     sigma.m=NULL, prop.Sigma=NULL, 
#                     true.params=true.params, lower.prob.lim=0.5)

