################################################################################
## The log likelihood of the data, where the data comes from a scale mixture
## of Gaussians, transformed to GPD (matrix/vector input)
##
## NOT ACTUALLY depending on X, just calculated Fx^-1(Fy(Y)) ahead of time
##
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X ................................. a matrix of Y, transformed to HOT 
##                                     scale mixture of Gaussian
## X.s ............................... X, but without the measurement error
## cen 
## prob.below
## theta.gpd=(u,a,b,shape)
## delta
## tau_sqd
## lon
##
marg.transform.data.mixture.me.likelihood <- function(Y, X, X.s, cen, prob.below,
                                                      theta.gpd, shape, lon, gamma, sigma,
                                                      tau_sqd, thresh.X=NULL, missing=NULL) {
  if (!is.matrix(Y)) Y <- matrix(Y, ncol=1)
  if (!is.matrix(X)) X <- matrix(X, ncol=1)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  if(is.null(missing)) missing <- matrix(FALSE, nrow=nrow(X),ncol=ncol(X))
  
  ll <- matrix(0, nrow(Y), ncol(Y))
  
  loc <- theta.gpd[1]
  a <- theta.gpd[2]
  b <- theta.gpd[3]
  Scale <- a+b*lon
  
  if (is.null(thresh.X)) thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau_sqd, gamma = gamma, sigma = sigma)
  
  if (sum(cen) > 0)  {
    ll[cen] <- pnorm(thresh.X, mean=X.s[cen], sd=sqrt(tau_sqd), log.p=TRUE)
  }
  if (sum(!cen) > 0) {
    ll[!cen & !missing] <- dnorm(X[!cen & !missing], mean=X.s[!cen& !missing], sd=sqrt(tau_sqd), log=TRUE) +
      dgpd(Y[!cen & !missing], loc=loc, scale=Scale[!cen & !missing], shape=shape, log=TRUE)  -
      dmixture.me(X[!cen & !missing], tau_sqd = tau_sqd, gamma=gamma, sigma=sigma, log=TRUE)
    ll[!cen & missing] <- pnorm(thresh.X, mean=X.s[!cen & missing], sd=sqrt(tau_sqd), lower.tail = FALSE, log.p=TRUE)
  }
  
  which <- is.na(ll)
  if(any(which)) ll[which] <- -Inf  # Normal density larger order than marginal density of scalemix
  return(sum(ll)) 
}

#                                                                              #
################################################################################





################################################################################
## The log likelihood of X.s, when it is conditioned on scaling factor R and Σ(λ,γ)
##
## X.s ............................... X, but without the measurement error
## R 
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

X.s.likelihood.conditional<-function(X.s, R, V, d){
  
    X.s.to.Z <- X.s/R
    loglik <- -0.5*as.vector(eig2inv.quadform.vector(V, 1/d, X.s.to.Z))-0.5*sum(log(d))-length(X.s)*log(R)
    return(loglik)
  
}


X.s.likelihood.conditional.on.X<-function(X.s, X, R, V, d, tau_sqd){
  if(any(X.s<R)) return(-Inf) else{
    X.s.to.Z <- qnorm(1-R/X.s)
    loglik <- -0.5*as.vector(eig2inv.quadform.vector(V, 1/d, X.s.to.Z))+0.5*sum(X.s.to.Z^2)-2*sum(log(X.s))-0.5*sum(log(d))+length(X.s)*log(R)
    loglik <- loglik - 0.5*sum((X.s-X)^2)/tau_sqd
    return(loglik)
  }
}

var.at.a.time.update.X.s <- function(X, R, V, d, tau_sqd, v.q=0.5, n.chain=100){
  X <- as.vector(X)
  n.s <- length(X)
  X.s <- X
  X.s[which(X.s<R)] <- R + abs(rnorm(length(which(X.s<R)),sd=0.5)) # X.s has to be larger than R
  accept <- rep(0, n.s)
  
  for(i in 1:n.chain){
    for(iter in 1:n.s){
      X.s.update<-X.s
      X.s.update[iter]<-X.s[iter]+rnorm(1,0,sd=v.q)
      log.num <- X.s.likelihood.conditional.on.X(X.s.update, X=X, R=R, V=V, d=d, tau_sqd=tau_sqd)
      log.denom <- X.s.likelihood.conditional.on.X(X.s, X=X, R=R, V=V, d=d, tau_sqd=tau_sqd)
      
      r <- exp(log.num - log.denom)
      
      if(runif(1) < r){    
        X.s <- X.s.update
        accept[iter]<-accept[iter]+1
      }
    }
  }
  
  return(list(X.s=X.s, accept=accept))
}

#                                                                              #
################################################################################





################################################################################
## For the generic Metropolis sampler
## Samples from the scaled Gaussian process (update the smooth process).
## The mixing distribution comes from from the Huser-wadsworth scale mixing distribution.
## The PROPOSAL will be the conjugate update from the model that treats the
## X process (i.e. X.s, but with measurement error) as the response, ignoring
## the marginal transformation.  Then the Metropolis part either rejects or
## accepts the entire draw.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters (theta in dhuser.thibaud):
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## theta.mix
## tau_sqd
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##


X.s.update.mixture.me.update.par.once.without.X.par <- function(R, Y, X, X.s, cen,
                                                            prob.below, theta.gpd, shape, lon, gamma, sigma,
                                                            tau_sqd, V, d, v.q=NULL, 
                                                            thresh.X=NULL, missing = NULL){
  # library(doParallel)
  # library(foreach)
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  # cores <- detectCores()
  # registerDoParallel(cores=cores)
  
  accepted <- matrix(0,nrow=n.s, ncol=n.t) 
  if(is.null(v.q)) v.q<-matrix(2.4^2,nrow=n.s, ncol=n.t)
  if(is.null(missing)) missing <- matrix(FALSE, nrow=nrow(X),ncol=ncol(X))
  if(is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, gamma = gamma, sigma = sigma)
  
  Res<-foreach(t = 1:n.t, .combine = "cbind", .inorder=FALSE, .noexport = "update_X_s_onetime")%dopar%{
    res<- update_X_s_onetime(Y[,t], X[,t], X.s[,t], cen[,t], prob.below, theta.gpd, shape, lon[,t], gamma, sigma, 
                             tau_sqd, thresh.X, v.q[,t], R[t], V, d, missing[,t])
    c(res$X.s, res$accept, t)
  }
  
  X.s <- Res[1:n.s, ]
  accepted <-Res[(n.s+1):(2*n.s),]
  O<-order(Res[2*n.s+1,])
  X.s <- X.s[,O]
  accepted <- accepted[,O]
  
  return(list(X.s=X.s, accepted=accepted))
}


#                                                                              #
################################################################################





################################################################################
## Updates the censored Xs, just by drawing from the log likelihood of the
## data, censored at the threshold (on the Y scale), where the data comes
## from a Huser-wadsworth scale mixture
## of transformed Gaussians, transformed to GPD
##
## X.s ............................... X, but without the measurement error
## cen 
## prob.below
## theta.mix
## theta.gaussian
##
update.censored.obs.mixture.me <- function(X.s, cen, tau_sqd, thresh.X) {
  
  ll <- matrix(0, nrow(X.s), ncol(X.s))
  sd=sqrt(tau_sqd)
  
  # Draws a collection of univariate truncated normals
  B <- pnorm(thresh.X, X.s[cen], sd=sd)
  U <- B * runif(sum(cen))
  X.update <- qnorm(U, mean=X.s[cen], sd=sd)
  
  return(X.update) 
}

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the parameters of the mixing distribution, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. delta (from Huser and Wadsworth 2017)
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## tau_sqd
##

delta.update.mixture.me.likelihood <- function(data, params, Y, X.s, cen, 
                                                   prob.below, theta.gpd,
                                                   shape, lon, tau_sqd, missing = NULL) {
  R <- data
  gamma <- params[1]
  sigma <- params[2]
  if(gamma < 0 || sigma < 0) return(-Inf)
  if(is.null(missing)) missing <- matrix(FALSE, nrow=nrow(X.s),ncol=ncol(X.s))
  
  X <- NA * Y
  X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd=tau_sqd, gamma = gamma, sigma = sigma, 
                               theta.gpd=theta.gpd, shape=shape, lon=lon[!cen& !missing], prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
               theta.gpd, shape, lon, gamma, sigma, tau_sqd, missing = missing) + dhuser.wadsworth(R, gamma, sigma, log=TRUE)

  return(ll)
}
# delta.update.mixture.me.likelihood(R, delta, Y, X.s, cen, prob.below, theta.gpd, tau)

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the measurement error variance (on the X scale), for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
## Just a wrapper for marg.transform.data.mixture.me.likelihood
##
##   *********** If we do end up updating the prob.below parameter, this
##   *********** is a good place to do it.
##
## data............................... a n.t vector of scaling factors
## params............................. tau_sqd
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## delta
##
tau.update.mixture.me.likelihood <- function(data, params, Y, X.s, cen, 
                                                 prob.below, gamma, sigma, 
                                                 theta.gpd, shape, lon, missing = NULL) {
  
  if(is.null(missing)) missing <- matrix(FALSE, nrow=nrow(X.s),ncol=ncol(X.s))
  R <- data
  tau_sqd <- params

  X <- NA * Y
  X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd=tau_sqd, gamma=gamma, sigma=sigma, 
                               theta.gpd=theta.gpd, shape=shape, lon=lon[!cen & !missing], prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
                                                  theta.gpd, shape, lon, gamma, sigma, tau_sqd, missing = missing)
  
  return(ll)
}
#tau.update.mixture.me.likelihood(R, tau, Y, X.s, cen, prob.below, delta, theta.gpd)

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the parameters of the GPD response distribution, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
## Just a wrapper for marg.transform.data.mixture.me.likelihood
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters:
##                                     theta[1] = GPD scale
##                                     theta[2] = GPD shape
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## delta
## tau_sqd
##
theta.gpd.update.mixture.me.likelihood <- function(data, params, shape, lon, Y, X.s, cen, 
                                                   prob.below, gamma, sigma,
                                                   tau_sqd, loc, thresh.X=NULL, missing = NULL) {
  
  R <- data
  theta.gpd <- c(loc, params)
  if(is.null(missing)) missing <- matrix(FALSE, nrow=nrow(X.s),ncol=ncol(X.s))
  
  a <- params[1]
  b <- params[2]
  Scale <- a+b*lon
  if (shape >= 0) max.support <- matrix(Inf, nrow=nrow(lon), ncol=ncol(lon))  else max.support <- loc - Scale/shape
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if (any(Y[!cen & !missing] > max.support[!cen & !missing]) | any(Scale[!cen]<0)) return(-Inf)
  
  X <- NA * Y
  X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd=tau_sqd, gamma=gamma, sigma=sigma, 
                               theta.gpd=theta.gpd, shape=shape, lon=lon[!cen & !missing], prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
                                  theta.gpd, shape, lon, gamma, sigma, tau_sqd, thresh.X = thresh.X, missing = missing)
  
  return(ll)
}


shape.gpd.update.mixture.me.likelihood <- function(data, params, theta.gpd, lon, Y, X.s, cen, 
                                                   prob.below, gamma, sigma,
                                                   tau_sqd, loc, thresh.X=NULL, missing = NULL) {
  
  if(is.null(missing)) missing <- matrix(FALSE, nrow=nrow(X.s),ncol=ncol(X.s))
  R <- data
  loc <- theta.gpd[1]
  a <- theta.gpd[2]
  b <- theta.gpd[3]
  
  shape <- params[1]
  Scale <- a+b*lon
  if (shape >= 0) max.support <- matrix(Inf, nrow=nrow(lon), ncol=ncol(lon))  else max.support <- loc - Scale/shape
  
  # If the parameters imply support that is not consistent with the data,
  # then reject the parameters.
  if (any(Y[!cen & !missing] > max.support[!cen & !missing])) return(-Inf)
  
  X <- NA * Y
  X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd=tau_sqd, gamma=gamma, sigma=sigma, 
                               theta.gpd=theta.gpd, shape=shape, lon=lon[!cen & !missing], prob.below = prob.below)
  
  ll <- marg.transform.data.mixture.me.likelihood(Y, X, X.s, cen, prob.below,
                                        theta.gpd, shape, lon, gamma, sigma, tau_sqd, thresh.X = thresh.X, missing = missing)
  
  return(ll)
}


#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the parameters of the underlying Gaussian process, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters:
##                                     lambda = params[1]
##                                     gamma  = params[2]
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## R
## S
## V, d
##
lam.gam.update.mixture.me.likelihood <- function(data, params, X.s, R, S, 
                                                 V=NULL, d=NULL) {
  
  library(doParallel)
  library(foreach)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  n.t <- ncol(X.s)
  registerDoParallel(cores=n.t)
  
  R <- data
  lambda <- params[1]
  gamma <- params[2]
  if(lambda<0 || gamma<0 || gamma>2)  return(-Inf)
  
  if(is.null(V)){
    Cor   <- corr.fn(rdist(S), lambda = lambda, gamma = gamma)
    eig.Sigma <- eigen(Cor, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values
  }
  
  ll<-foreach(i = 1:n.t, .combine = "c") %dopar% {
    X.s.likelihood.conditional(X.s[,i], R[i], V, d)
    # dmvn.eig(qnorm(1-R[i]/X.s[, i]), V = V, d.inv = 1/d)
  }
  
  return(sum(ll))
}

rho.update.mixture.me.likelihood <- function(data, params, X.s, R, S, 
                                                 V=NULL, d=NULL) {
  
  # library(doParallel)
  # library(foreach)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  n.t <- ncol(X.s)
  # registerDoParallel(cores=n.t)
  
  R <- data
  rho <- params[1]
  nu <- params[2]

  if(is.null(V)){
    Cor   <- corr.fn(rdist(S), rho = rho, nu = nu)
    eig.Sigma <- eigen(Cor, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values
  }
  
  ll<-foreach(i = 1:n.t, .combine = "c") %dopar% {
    X.s.likelihood.conditional(X.s[,i], R[i], V, d)
    # dmvn.eig(qnorm(1-R[i]/X.s[, i]), V = V, d.inv = 1/d)
  }
  
  return(sum(ll))
}
# lam.gam.update.mixture.me.likelihood(R, c(0.5, 1), X.s, R, S)
# lam.gam.update.mixture.me.likelihood(R, c(0.5, 1), X.s, R, S, V, d)

#                                                                              #
################################################################################



################################################################################
## For the generic Metropolis sampler
## Samples from the scaling factors, for the scale 
## mixture of Gaussians, where the   
## mixing distribution comes from Huser-wadsworth scale mixture.
##
## data............................... a n.t vector of scaling factors
## params............................. R[t]
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error (vector)
## V, d
## delta
##
Rt.update.mixture.me.likelihood <- function(data, params, X.s, gamma, sigma, 
                                                 V=NULL, d=NULL) {
  R <- data
  Rt <- params
  if(Rt < 1)  return(-Inf)
  
  ll <- X.s.likelihood.conditional(X.s, Rt, V, d)
  return(as.vector(ll))
}

# Rt.update.mixture.me.likelihood(R, R[1], X.s[,1], delta, V, d)

#                                                                              #
################################################################################




