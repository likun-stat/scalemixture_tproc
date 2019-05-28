Rcpp::sourceCpp('~/Desktop/Research/scalemixture_tproc/integration.cpp')
#Rcpp::sourceCpp('~/Desktop/Research/scalemixture/likelihood.cpp')


pow<-function(x,y){
  return(x^y)
}


################################################################################
################################################################################
## Compute the Matern nu-3/2 correlation function from a matrix of pairwise distances
## and a vector of parameters
##
## d ................................. a square matrix of distances
## lambda ............................ powered exponential
## gamma ............................. powered exponential 
##
## theta ............................. a 2-vector of correlation parameters:
##                                     theta[1] = log(rho), where rho is the 
##                                     range parameter
##                                     theta[2] is ignored

# corr.fn <- function(d, lambda, gamma) {
#   
#   n<- nrow(d)
#   Sigma<-matrix(NA,n,n)
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       Sigma[i,j] <- exp(- pow(d[i,j]/lambda, gamma))
#       Sigma[j,i] <- Sigma[i,j]
#     }
#   }
#   
#   for (k in 1:n) {
#     Sigma[k,k] <- 1
#   }
#   
#   return(Sigma)
# }

# corr.fn <- function(d, theta) {
#   
#   rho <- exp(theta[1])
#   
#   Sigma <- (1 + d/rho) * exp(-d/rho) 
#   return(Sigma)
# }

# corr.fn <- function(d, rho) {
#   
#   Sigma <- (1 + d/rho) * exp(-d/rho) 
#   return(Sigma)
# }

corr.fn <- function(d, rho, nu){
  Sigma <- fields::Matern(d, range=rho, smoothness = nu)
  return(Sigma)
}

##
################################################################################
################################################################################





################################################################################
################################################################################
## Approximates the marginal quantile function by taking values of the
## marginal CDF of X and doing linear interpolation.  If no values of the CDF
## are supplied, it computes n.x of them, for x in (lower, upper).
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
## mu ................................ The mean of the normal that gets mixed
##                                     This should probably always be zero
##                                     if the model is used as a copula.
## cdf.vals .......................... If the CDF has already been calculated,
##                                     then they can be passed in and re-used.
## x.vals ............................ Design points for the numerical
##                                     approximation.  Probably best to keep
##                                     as NULL and let the function figure out
##                                     good design points.
## n.x ............................... The number of design points to evaluate
##                                     the CDF.
## lower ............................. The best guess as the smallest design
##                                     point that is necessary.  If this is
##                                     too large, it will get fixed in the 
##                                     function.
## upper ............................. The best guess as the largest design
##                                     point that is necessary.  If this is
##                                     too small, it will get fixed in the 
##                                     function.
##
qmixture.me.interp <- function(p, tau_sqd, gamma, sigma, mu=0, cdf.vals = NULL, x.vals = NULL,
                               n.x=10000, lower=5, upper=20, method="fmm") {
  
  if (is.null(x.vals)) {
    x.range <- find_xrange_pmixture_me(min(p), max(p), c(lower, upper), 
                                       tau_sqd, gamma, sigma, relerr = 1e-10)
    # x.vals <- seq(x.range[1], x.range[2], length=n.x)
    if(x.range[2]==Inf) {x.range[2] <- 10^20; large.delta.large.x <- TRUE}
    x.vals <- exp(seq(log(x.range[1]), log(x.range[2]), length=n.x))
    cdf.vals <- pmixture_me(x.vals, tau_sqd, gamma, sigma)
  } else {
    if (is.null(cdf.vals)) {
      cdf.vals <- pmixture_me(x.vals, tau_sqd, gamma, sigma)
    }
  }
  
  q.vals <- spline(x=cdf.vals, y=x.vals, xout=p, method=method)$y 
  return(q.vals)
  
}
##
################################################################################
################################################################################






################################################################################
################################################################################
## The thing that gets integrated dr to result in the marginal CDF of X
##
## 
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
##
mix.distn.integrand <- function(x, xval, tau_sqd, gamma, sigma) {
  prod <- rep(0, length(x))

  half_result = pt((xval-x)/sigma, gamma)
  prod <- dnorm(x, 0.0, sqrt(tau_sqd)) * half_result
  return(prod)
}

mix.dens.integrand <- function(x, xval, tau_sqd, gamma, sigma) {
  prod <- rep(0, length(x))
  
  half_result = dt((xval-x)/sigma, gamma)
  prod <- dnorm(x, 0.0, sqrt(tau_sqd)) * half_result
  return(prod)
}
##
################################################################################
################################################################################






################################################################################
################################################################################
## Integrates mix.distn.integrand dr to result in the marginal CDF of X
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
##
pmixture.uni <- function(xval, tau_sqd, gamma, sigma) {
  integ <- integrate(mix.distn.integrand, -Inf, Inf,  xval=xval, tau_sqd = tau_sqd, gamma=gamma, sigma=sigma, rel.tol = 1e-10)
  return(pnorm(xval-1, 0.0, sqrt(tau_sqd))-integ$value)
}
pmixture <- Vectorize(pmixture.uni, "xval")

dmixture.uni <- function(xval, tau_sqd, gamma, sigma) {
  integ <- integrate(mix.dens.integrand, -Inf, Inf,  xval=xval, tau_sqd = tau_sqd, gamma=gamma, sigma=sigma, rel.tol = 1e-10)
  return(integ$value)
}
dmixture <- Vectorize(dmixture.uni, "xval")
##
################################################################################
################################################################################



################################################################################
################################################################################
## Evaluates the marginal PDF of X.  If the vactor x has length greater than
## max.x, it will do the integral at max.x points, and interpolate them to
## the values of x.
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
##
dmixture.me <- function(x, tau_sqd, gamma, sigma, log=FALSE, max.x=8000) {
  if (length(x) < max.x) {
    dens <- dmixture_me(x, tau_sqd, gamma, sigma)
  } else {
    design.x <- seq(min(x), max(x), length=max.x)
    design.y <- dmixture_me(design.x, tau_sqd, gamma, sigma) 
    dens <- spline(x=design.x, y=design.y, xout=x)$y
  }
  
  # if(any(x<tau_sqd/5000)) {
  #   dens[x<tau_sqd/5000] <- tryCatch(dmixture(x[x<tau_sqd/5000], tau_sqd, gamma, sigma), 
  #                                    error=function(e){cat("delta=",delta,"tau=",tau_sqd,"\n");
  #                                      rep(1e-20,length(x<tau_sqd/5000))})
  # }
  if (!log) return(dens) else return(log(dens))
}

##
################################################################################
################################################################################






################################################################################
################################################################################
## Transforms observations from a Gaussian scale mixture to a GPD
##
## x ................................. A vector of observations from a HOT
##                                     scale mixture
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
## theta.gpd ......................... A vector, (thresh, scale, shape)
## 
##
scalemix.me.2.gpd <- function(x, tau_sqd, gamma, sigma, theta.gpd, shape, lon, prob.below=0) {
  require(evd)
  
  thresh <- theta.gpd[1]
  a <- theta.gpd[2]
  b <- theta.gpd[3]
  Scale <- a+b*lon
  
  unifs <- pmixture_me(x, tau_sqd, gamma, sigma)
  # gpds <- qgpd(unifs, loc=thresh, scale=scale, shape=shape)
  gpds <- qgpd((unifs-prob.below) / (1-prob.below), loc=thresh, scale=Scale, shape=shape)
  
  return(gpds)
}
##
################################################################################
################################################################################


# 
# u = (1-p) * p(Y) + P
# (u-p) / (1-p) = p(y)



################################################################################
################################################################################
## Transforms observations from a GPD to a Gaussian scale mixture
##
## y ................................. A vector of observations from a GPD
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
## theta.gaussian .................... A parameter, mu
## 
##
gpd.2.scalemix.me <- function(y, tau_sqd, gamma, sigma, theta.gpd, shape, lon, prob.below=0) {
  require(evd)
  
  thresh <- theta.gpd[1]
  a <- theta.gpd[2]
  b <- theta.gpd[3]
  Scale <- a+b*lon
  
  # unifs <- pgpd(y, loc=thresh, scale=scale, shape=shape)
  unifs <- (1-prob.below) * pgpd(y, loc=thresh, scale=Scale, shape=shape) + prob.below
  scalemixes <- qmixture.me.interp(unifs, tau_sqd = tau_sqd, gamma=gamma, sigma=sigma, n.x=10000)
  
  return(scalemixes)
}
##
################################################################################
################################################################################



## -------------------------------------------------------------------------------------------
## --------------------------------  For Matric Algebra  -------------------------------------
## -------------------------------------------------------------------------------------------



################################################################################
################################################################################
## Multivariate normal log density of R, where each column of 
## R iid N(0,V'DV), where V'DV is the covariance matrix
## It essentially computes the log density of each column of R, then takes
## the sum.  Faster than looping over the columns, but not as transparent.
##
## Is there a faster way to to this? i.e. by using the quadform function?
##
## R ................................. a (s x t) matrix, where each column is
##                                     multivariate normal, independently from
##                                     the other columns
## V ................................. a (s x s) matrix of eigenvectors of the
##                                     covariance matrix of the columns of R
## d.inv ............................. a (s) vector of the inverted eigenvalues,
##                                     d.inv_i = 1/D_ii
##
dmvn.eig <- function(R, V, d.inv){
  if (is.vector(R)) n.rep <- 1 else n.rep <- ncol(R)
  return(-0.5*n.rep*eig2logdet(1/d.inv) - 0.5 * sum(R * eig2inv.times.vector(V, d.inv, R)))
} 
##
################################################################################
################################################################################


################################################################################
################################################################################
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
## Computes A^{-1}x
##
## V ................................. a matrix of eigenvectors, as above
## d.inv ............................. a vector of the inverted eigenvalues, with
##                                     d.inv_i = 1/D_ii
## x ................................. a vector x, as above
##
eig2inv.times.vector <- function(V, d.inv, x) {
  return(V %*% (d.inv * crossprod(V, x)))
}
##
################################################################################
################################################################################



################################################################################
################################################################################
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
##
## Computes x'A^{-1}x
##
## V ................................. a matrix of eigenvectors, as above
## d.inv ............................. a vector of the inverted eigenvalues, with
##                                     d.inv_i = 1/D_ii
## x ................................. a vector x, as above
##
eig2inv.quadform.vector <- function(V, d.inv, x) {
  cp <- crossprod(V, x)
  return(crossprod(cp, d.inv*cp))
}
##
################################################################################
################################################################################



################################################################################
################################################################################
## Assumes that A = VDV', where D is a diagonal vector of eigenvectors of A, and
## V is a matrix of normalized eigenvectors of A.
##
## log(|A|)
##
## d ................................. a vector of the eigenvalues, with
##                                     d_i = D_ii
##
eig2logdet <- function(d) {
  return(sum(log(d)))
}
##
################################################################################
################################################################################



################################################################################
################################################################################
## Compute the density of the mixing distribution in Huser-Wadsworth (2017).
## delta in the notation in the manuscript.
## The support of r is [0, infinity)
##
## R ................................. value at which to evaluate the density
## delta 
##
dhuser.wadsworth <- function(R, gamma, sigma, log=TRUE) {
  require(invgamma)
  if(!all(R > 0)) return(-Inf)
  n.t <- length(R)
  
  if (log) {
    dens <- sum(dinvgamma(R^2/sigma^2, shape=gamma/2, scale=2, log=TRUE)+log(2*R/sigma^2))
  } else {
    dens <- prod(dinvgamma(R^2/sigma^2, shape=gamma/2, scale=2)*2*R/sigma^2)
  }
  return(dens)
}
##
################################################################################
################################################################################

