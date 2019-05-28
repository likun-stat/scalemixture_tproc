

################################################################################
#  Computes two independent half Cauchy densities with scale hyper.params.     #

half.cauchy.half.cauchy <- function(params, hyper.params) {
  theta.1 <- params[1]
  theta.2 <- params[2]
  scale.theta.1 <- hyper.params[1]
  scale.theta.2 <- hyper.params[2]
  
  if (min(params)<0) return(-Inf)
  return(log(2) + dcauchy(theta.1, scale=scale.theta.1, log=TRUE) + 
         log(2) + dcauchy(theta.2, scale=scale.theta.2, log=TRUE))
}
#                                                                              #
################################################################################





################################################################################
#  Computes a logit normal prior for a scalar parameter.                       #

normal.scalar <- function(params, hyper.params) {
  theta <- params[1]
  mean.theta <- hyper.params[1]
  sd.theta <- hyper.params[2]
  
  return(dnorm(theta,mean= mean.theta, sd=sd.theta, log=TRUE))
}
#                                                                              #
################################################################################



################################################################################
#  Computes a half Cauchy densities with scale hyper.params for the first      #
#  of the params, and a uniform(-0.5, 0.5) for the second.                     #

half.cauchy.scale.unif.shape <- function(params, hyper.params) {
  scale <- params[1]
  shape <- params[2]
  scale.theta.1 <- hyper.params[1]
  
  if (scale < 0) return(-Inf)
  if ((shape < -0.5) | (shape > 0.5)) return(-Inf)
  return(log(2) + dcauchy(scale, scale=scale.theta.1, log=TRUE))
}
#                                                                              #
################################################################################



################################################################################
#  Computes a half Cauchy densities with b for the slope (b<0),                #
#  and a uniform a for the intercept.                                          #

unif.a.b <- function(params, hyper.params) {
  a <- params[1]
  b <- params[2]
  maxX <- hyper.params[1]
  
  if (b > 0 | maxX>= -a/b) return(-Inf)
  return(1)
}

unif.shape <- function(params, hyper.params=1) {
  shape <- params[1]
  scale.theta.1 <- hyper.params[1]
  
  if ((shape < -0.5) | (shape > 0.5)) return(-Inf)
  return(1)
}
#                                                                              #
################################################################################



################################################################################
#  Computes a half Cauchy density with scale hyper.params.                     #

half.cauchy <- function(params, hyper.params) {
  sigma <- params[1]
  scale.theta.1 <- hyper.params[1]
  
  if (sigma < 0) return(-Inf)
  return(log(2) + dcauchy(sigma, scale=scale.theta.1, log=TRUE))
}
#                                                                              #
################################################################################


################################################################################
#  Computes the Huser-wadsworth prior for R. (one dim)                         #

huser.wadsworth.prior <- function(params, hyper.params) {
  R <- params[1]
  gamma <- hyper.params[1]
  sigma <- hyper.params[2]
  return(dhuser.wadsworth(R, gamma, sigma, log=TRUE))
}
#                                                                              #
################################################################################


################################################################################
#  Computes a (0,w) unif with hyper.params w.                                  #

interval.unif <- function(params, hyper.params) {
  gamma <- params[1]
  sigma <- params[2]
  scale <- hyper.params
  
  if (gamma<0 | sigma <0) return(-Inf)
  return(log(2) + dcauchy(gamma, scale=scale, log=TRUE)  + dcauchy(sigma, scale=scale, log=TRUE))
  
}
#                                                                              #
################################################################################


################################################################################
#  Computes a (0,∞)×(0,w) unif with hyper.params w.                            #

lam.gam.prior <- function(params, hyper.params) {
  lambda <- params[1]
  gamma <- params[2]
  w <- hyper.params[1]
  
  if(lambda<0 || gamma<0 || gamma>w)  return(-Inf)
  
  return(-log(w))
}
#                                                                              #
################################################################################


################################################################################
#  Computes inverse gamma distribution with hyper.params w.                    #

tau.sqd.prior <- function(params, hyper.params) {
  tau_sqd <- params[1]
  alpha <- hyper.params[1]
  beta <- hyper.params[2]
  
  if(tau_sqd < 0 || tau_sqd > 1e6)  return(-Inf)
  
  return((-alpha-1)*log(tau_sqd)-beta/tau_sqd)
}


#                                                                              #
################################################################################


################################################################################
#  Computes independent half Cauchy density and Cauchy density                 #

gpd.prior <- function(params, hyper.params) {
  scale <- params[1]
  shape <- params[2]
  w <- hyper.params[1]

  if (scale<w) return(-Inf)
  return(0)
}

#  Computes a half Cauchy densities with scale hyper.params for the first      #
#  of the params, and a uniform(-0.5, 0.5) for the second.                     #

half.cauchy.scale.unif.shape <- function(params, hyper.params) {
  scale <- params[1]
  shape <- params[2]
  scale.theta.1 <- hyper.params[1]
  
  if (scale < 0) return(-Inf)
  if ((shape < -0.5) | (shape > 0.5)) return(-Inf)
  return(log(2) + dcauchy(scale, scale=scale.theta.1, log=TRUE))
}
#                                                                              #
################################################################################


################################################################################
#  Computes prior for rho                                                      #

log.rho.prior <- function(params, hyper.params) {
  rho <- params
  w <- hyper.params
  
  return(0)
}

rho.prior <- function(params, hyper.params) {
  rho <- params[1]
  nu <- params[2]
  scale <- hyper.params

  if (rho<0 | nu <0) return(-Inf)
  return(log(2) + dcauchy(rho, scale=scale, log=TRUE)  + dcauchy(nu, scale=scale, log=TRUE))
}

#                                                                              #
################################################################################


################################################################################
#  Computes a uniform (0, 1) for delta, a half Cauchy densities with scale     #
#  of the params, and a uniform(-0.5, 0.5) for the shape  .                    #

delta.gpd.prior <- function(params, hyper.params) {
  delta <- params[1]
  scale <- params[2]
  shape <- params[3]
  scale.theta.1 <- hyper.params[1]
  
  if((delta<0) | (delta>1)) return(-Inf)
  if (scale < 0) return(-Inf)
  if ((shape < -0.5) | (shape > 0.5)) return(-Inf)
  return(log(2) + dcauchy(scale, scale=scale.theta.1, log=TRUE))
}

#                                                                              #
################################################################################
