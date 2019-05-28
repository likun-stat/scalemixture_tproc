library(fields)
library(scales)


# For data analysis, save less entries for X.s.trace

scalemix.sampler.02.DA <- function(Y, S, cen, lon, thresh,
                                   initial.values,
                                   n.updates, thin=10,
                                   experiment.name="Huser-wadsworth",
                                   echo.interval=50,
                                   sigma.m=NULL, prop.Sigma=NULL, missing = NULL,
                                   true.params=NULL, sd.ratio=NULL, sd.ratio.rho=NULL,sd.ratio.delta=NULL, lower.prob.lim=0.5,
                                   holdout=FALSE, S_full = NULL, n.holdout = NULL) {
  
  #library(doParallel)
  #library(foreach)
  save.bit <- TRUE
  
  # Constants to control how many Metropolis steps get executed at each
  # iteration of the main loop
  n.metr.updates.delta <- 4
  n.metr.updates.rho <- 4
  n.metr.updates.theta.gpd <- 4
  n.metr.updates.shape <- 4
  n.metr.updates.prob.below <- 4
  n.metr.updates.tau <- 2
  n.metr.updates.R <- 2
  
  # Constants to control adaptation of the Metropolis sampler
  c.0 <- 10
  c.1 <- 0.8
  k <- 3  # the iteration offset
  metr.opt.1d <- 0.41
  metr.opt.2d <- 0.35
  
  # Hyper parameters for the prior of the mixing distribution parameters and 
  # the correlation parameters
  hyper.params.delta <- 1
  hyper.params.rho <- 1
  hyper.params.theta.gpd <- max(lon[!cen])
  hyper.params.shape <- 1
  hyper.params.tau <- c(0.1,0.1)
  # hyper.params.prob.below <- c(0, 10)
  
  # A small number
  eps <- 1e-06
  
  # Bookkeeping
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  h <- rdist(S)
  diag(h)  <- 0
  if(holdout) {h_full <- rdist(S_full); diag(h_full) <- 0}
  #registerDoParallel(cores=n.t)
  if(is.null(missing)) missing <- matrix(FALSE, nrow=n.s,ncol=n.t)
  
  # Load initial values
  gamma <- initial.values$gamma
  sigma <- initial.values$sigma
  rho <- initial.values$rho
  nu <- initial.values$nu
  tau <- initial.values$tau
  a <- initial.values$theta.gpd[2]
  b <- initial.values$theta.gpd[3]
  shape <- initial.values$shape
  prob.below <- initial.values$prob.below
  # logit.prob.below <- logit(prob.below, c(lower.prob.lim, 1))
  X.s <- initial.values$X.s
  R <- initial.values$R
  
  X.s.accept <- rep(1, n.t)
  
  theta.gpd <- c(thresh, a, b)
  
  
  #**********   Eigendecomposition of the correlation matrix   ****************#
  Sigma <- corr.fn(h, rho, nu)
  eig.Sigma <- eigen(Sigma, symmetric=TRUE)
  V <- eig.Sigma$vectors
  d <- eig.Sigma$values
  if(holdout) Sigma_full <- corr.fn(h_full, rho, nu)
  
  # Initialize X and X.s
  thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, gamma=gamma, sigma=sigma)
  X <- Y
  X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd = tau, gamma=gamma, sigma=sigma, theta.gpd = theta.gpd, shape=shape, lon=lon[!cen & !missing], prob.below=prob.below)
  X[cen] <- update.censored.obs.mixture.me(X.s, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
  
  # Initialize trace objects
  accepted <- rep(0, n.updates)
  
  X.s.trace <- array(NA, dim=c(n.updates, n.s, 9))
  X.trace <- array(NA, dim=c(n.updates, n.s, 9))
  X.s.accept.trace <- matrix(NA, n.updates, n.t)
  rho.trace <- matrix(NA, n.updates, 2)
  tau.trace <- rep(NA, n.updates)
  delta.trace <- matrix(NA, n.updates,2)
  theta.gpd.trace <- matrix(NA, n.updates, 2)
  shape.trace <- rep(NA, n.updates)
  # prob.below.trace <- rep(NA, n.updates)
  R.trace <- matrix(NA, n.updates, n.t)
  
  sd.ratio.trace <- rep(NA, n.updates)
  
  # Column names for trace objects
  colnames(theta.gpd.trace) <- c("a", "b")
  if(holdout) X.s.holdout.trace <- array(NA, dim=c(n.updates, n.holdout, n.t))
  
  # Fill in trace objects with initial values
  subset.replicates <- ceiling(seq(1,n.t,length.out = 9))
  X.s.trace[1, , ] <- X.s[,subset.replicates]
  X.trace[1, , ] <- X[,subset.replicates]
  X.s.accept.trace[1, ] <- X.s.accept
  rho.trace[1,] <- c(rho,nu)
  tau.trace[1] <- tau
  delta.trace[1,] <- c(gamma, sigma)
  theta.gpd.trace[1, ] <- theta.gpd[c(2,3)]
  shape.trace[1] <- shape
  # prob.below.trace[1] <- prob.below
  R.trace[1, ] <- R
  # if(is.null(prop.Sigma$theta.gpd))  prop.Sigma$theta.gpd<-diag(2)
  if(is.null(prop.Sigma$gpd.corr))  prop.Sigma$gpd.corr<-0
  if(is.null(prop.Sigma$rho.corr))  prop.Sigma$rho.corr<-0
  if(is.null(prop.Sigma$delta.corr))  prop.Sigma$delta.corr<-0
  sd.ratio.trace[1] <- prop.Sigma$theta.gpd[2, 2]
  
  # For tuning Metropolis updates of theta
  if (is.null(sigma.m$delta)) sigma.m$delta <- (2.4/2)^2
  if (is.null(sigma.m$rho)) sigma.m$rho <- (2.4/2)^2
  if (is.null(sigma.m$theta.gpd)) sigma.m$theta.gpd <- (2.4/2)^2
  if (is.null(sigma.m$shape)) sigma.m$shape <- 2.4^2
  # if (is.null(sigma.m$prob.below)) sigma.m$prob.below <- 1
  if (is.null(sigma.m$tau)) sigma.m$tau <- 2.4^2
  if (is.null(sigma.m$R)) sigma.m$R <- rep(2.4^2, n.t)
  if (is.null(sigma.m$X.s)) sigma.m$X.s <- matrix(2.4^2,nrow=n.s,ncol=n.t)
  
  r.hat.delta <- NA
  r.hat.rho <- NA
  r.hat.theta.gpd <- NA
  r.shape <- NA
  r.hat.prob.below <- NA
  r.hat.tau <- NA
  # r.hat.R <- rep(NA, n.t)
  
  if (is.null(sd.ratio))   sd.ratio <- 7
  if (is.null(sd.ratio.rho))   sd.ratio.rho <- 7
  if (is.null(sd.ratio.delta))   sd.ratio.delta <- 7
  
  for (i in 2:n.updates) {
    
    ################################################################
    ## Update Metropolis adaptation parameters
    ################################################################
    gamma1 <- c.0 / (i + k)^(c.1)
    gamma2 <- 1 / (i + k)^(c.1)
    Accepted<-matrix(0, nrow=n.s, ncol=n.t) #For X.s
    
    for (j in 1:thin) {
      
      ################################################################
      ## Update X -- branching: use X.s
      ################################################################
      # X[cen] <- update.censored.obs.mixture.me(X.s = X.s, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
      X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd = tau, gamma=gamma, sigma=sigma, theta.gpd = theta.gpd, 
                                   shape = shape, lon = lon[!cen & !missing], prob.below=prob.below)
      
      ################################################################
      ## Update X.star  *Parallel
      ################################################################
      # X.s.res <- X.s.update.mixture.me.par(R, Y, X, X.s, cen, 
      #                                      prob.below, theta.gpd, delta,
      #                                      tau, V, d, v.q=2, n.chain = 100, thresh.X = thresh.X)
      X.s.res<-X.s.update.mixture.me.update.par.once.without.X.par(R, Y, X, X.s, cen,
                                                                   prob.below, theta.gpd, shape, lon, gamma, sigma,
                                                                   tau, V, d, v.q=sigma.m$X.s, thresh.X=thresh.X, missing = missing)
      X.s <- X.s.res$X.s
      X.s.accept <- apply(X.s.res$accepted==1,2,any)
      Accepted <- tryCatch(Accepted + X.s.res$accepted, error=function(e){cat("The dim of X.s.res$accepted is (",dim(X.s.res$accepted),")","\n")})
      
      
      
      
      ################################################################
      ## Update rho
      ################################################################
      
      
      metr.out.rho <- static.metr(z = R, starting.theta = c(rho, nu),
                                  likelihood.fn = rho.update.mixture.me.likelihood, prior.fn = rho.prior,
                                  hyper.params = hyper.params.rho, n.updates = n.metr.updates.rho, prop.Sigma = prop.Sigma$rho, sigma.m=sigma.m$rho, verbose=FALSE,
                                  X.s = X.s, R = R, S = S)
      r.hat.rho <- metr.out.rho$acc.prob
      tmp <- metr.out.rho$trace[n.metr.updates.rho,]
      rho <- tmp[1]; nu <- tmp[2]
      sigma.m$rho <- exp(log(sigma.m$rho) + gamma1*(r.hat.rho - metr.opt.2d))
          
      
      ## Re-create covariance matrix and eigenvectors/eigenvalues
      Sigma   <- corr.fn(h, rho, nu)
      if(holdout) Sigma_full <- corr.fn(h_full, rho, nu)
      eig.Sigma <- eigen(Sigma, symmetric=TRUE)
      V <- eig.Sigma$vectors
      d <- eig.Sigma$values
      
      
      
      ################################################################
      ## Update R  *Parallel
      ################################################################
      Metr.R<-foreach(t = 1:n.t, .combine = "rbind") %dopar% {
        
        metr.out.R <- static.metr(z = R, starting.theta = R[t],
                                  likelihood.fn = Rt.update.mixture.me.likelihood, prior.fn = huser.wadsworth.prior,
                                  hyper.params = c(gamma, sigma), n.updates = n.metr.updates.R, prop.Sigma = 1, sigma.m=sigma.m$R[t], verbose=FALSE,
                                  X.s = X.s[,t], gamma=gamma, sigma=sigma, V = V, d = d)
        c(metr.out.R$trace[n.metr.updates.R],
          exp(log(sigma.m$R[t]) + gamma1*(metr.out.R$acc.prob - metr.opt.1d)))
      }
      
      R <- Metr.R[, 1]
      sigma.m$R <- Metr.R[, 2]
      
      
      
      ################################################################
      ## Update tau
      ################################################################
      metr.out.tau <-  static.metr(z = R, starting.theta = tau, 
                                   likelihood.fn = tau.update.mixture.me.likelihood, prior.fn = tau.sqd.prior,
                                   hyper.params = hyper.params.tau, n.updates = n.metr.updates.tau, prop.Sigma = 1, sigma.m=sigma.m$tau, verbose=FALSE, 
                                   Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, gamma=gamma, sigma=sigma, theta.gpd = theta.gpd, shape = shape, lon = lon, 
                                   missing=missing)
      tau <- metr.out.tau$trace[n.metr.updates.tau]
      r.hat.tau <- metr.out.tau$acc.prob
      sigma.m$tau <- exp(log(sigma.m$tau) + gamma1*(r.hat.tau - metr.opt.1d))
      
      
      
      ################################################################
      ## Update delta
      ################################################################
      metr.out.delta <- static.metr(z = R, starting.theta = c(gamma, sigma), 
                                    likelihood.fn = delta.update.mixture.me.likelihood, prior.fn = interval.unif,
                                    hyper.params = hyper.params.delta, n.updates = n.metr.updates.delta, prop.Sigma = prop.Sigma$delta, sigma.m=sigma.m$delta, verbose=FALSE, 
                                    Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, theta.gpd = theta.gpd, shape = shape, lon = lon, tau_sqd = tau,
                                    missing = missing)
      tmp <- metr.out.delta$trace[n.metr.updates.delta,]
      gamma <- tmp[1]; sigma <- tmp[2]
      r.hat.delta <- metr.out.delta$acc.prob
      sigma.m$delta <- exp(log(sigma.m$delta) + gamma1*(r.hat.delta - metr.opt.2d))
      
         
      
      
      ################################################################
      ## Update prob.below
      ################################################################
      #    metr.out.prob.below <- static.metr(R, logit.prob.below,
      #                                       logit.prob.below.update.mixture.me.likelihood,
      #                                       normal.scalar,
      #                                       hyper.params=hyper.params.prob.below,
      #                                       n.updates=n.metr.updates.prob.below,
      #                                       prop.Sigma=1,
      #                                       sigma.m=sigma.m$prob.below,
      #                                       Y=Y, X.s=X.s, cen=cen,
      #                                       theta.mix=theta.mix, theta.gaussian=theta.gaussian,
      #                                       theta.gpd=theta.gpd, lower.prob.lim=lower.prob.lim)
      #    logit.prob.below <- metr.out.prob.below$trace[n.metr.updates.prob.below]
      #    prob.below <- ilogit(logit.prob.below, c(lower.prob.lim, 1))
      #    r.hat.prob.below <- metr.out.prob.below$acc.prob
      #    sigma.m$prob.below <- exp(log(sigma.m$prob.below) +
      #                    gamma1*(r.hat.prob.below - metr.opt.1d))
      
      # Re-calculate the threshold on the X-scale    
      thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, gamma = gamma, sigma = sigma)
      
      
      
      ################################################################
      ## Update theta.gpd
      ################################################################
      metr.out.theta.gpd <- static.metr(z = R, starting.theta = theta.gpd[2:3],
                                        likelihood.fn = theta.gpd.update.mixture.me.likelihood, prior.fn = unif.a.b,
                                        hyper.params = hyper.params.theta.gpd, n.updates = n.metr.updates.theta.gpd, prop.Sigma = prop.Sigma$theta.gpd, sigma.m=sigma.m$theta.gpd, verbose=FALSE, 
                                        shape = shape, lon = lon, Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, gamma = gamma, sigma = sigma, tau_sqd = tau, loc = thresh, thresh.X=thresh.X,
                                        missing = missing)
      
      theta.gpd[c(2,3)] <- metr.out.theta.gpd$trace[n.metr.updates.theta.gpd, ]
      r.hat.theta.gpd <- metr.out.theta.gpd$acc.prob
      sigma.m$theta.gpd <- exp(log(sigma.m$theta.gpd) + gamma1*(r.hat.theta.gpd - metr.opt.2d))
      
      ################################################################
      ## Update shape
      ################################################################
      metr.out.shape <- static.metr(z = R, starting.theta = shape, 
                                    likelihood.fn = shape.gpd.update.mixture.me.likelihood, prior.fn = unif.shape,
                                    hyper.params = hyper.params.shape, n.updates = n.metr.updates.shape, prop.Sigma = 1, sigma.m=sigma.m$shape, verbose=FALSE, 
                                    theta.gpd = theta.gpd, lon = lon, Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, gamma = gamma, sigma = sigma, tau_sqd = tau, loc = thresh, 
                                    thresh.X=thresh.X, missing = missing)
      shape <- metr.out.shape$trace[n.metr.updates.shape]
      r.hat.shape <- metr.out.shape$acc.prob
      sigma.m$shape <- exp(log(sigma.m$shape) + gamma1*(r.hat.shape - metr.opt.1d))
      
    }
    
    
    
    # ---------------- Attempt to adapt proposal covariance ------------------
    if (r.hat.theta.gpd > 0) {
      sd.ratio.hat <- sd(metr.out.theta.gpd$trace[ ,1]) / sd(metr.out.theta.gpd$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio <- exp(log(sd.ratio) + gamma1*(log(sd.ratio.hat) - log(sd.ratio)))
    prop.Sigma$theta.gpd <-  matrix(c(1, prop.Sigma$gpd.corr/sd.ratio, prop.Sigma$gpd.corr/sd.ratio, 1/sd.ratio^2), 2, 2)
    
    
    if (r.hat.rho > 0) {
      sd.ratio.hat <- sd(metr.out.rho$trace[ ,1]) / sd(metr.out.rho$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio.rho <- exp(log(sd.ratio.rho) + gamma1*(log(sd.ratio.hat) - log(sd.ratio.rho)))
    prop.Sigma$rho <-  matrix(c(1, prop.Sigma$rho.corr/sd.ratio.rho, prop.Sigma$rho.corr/sd.ratio.rho, 1/sd.ratio.rho^2), 2, 2)
    
    
    if (r.hat.delta > 0) {
      sd.ratio.hat <- sd(metr.out.delta$trace[ ,1]) / sd(metr.out.delta$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio.delta <- exp(log(sd.ratio.delta) + gamma1*(log(sd.ratio.hat) - log(sd.ratio.delta)))
    prop.Sigma$delta <-  matrix(c(1, prop.Sigma$delta.corr/sd.ratio.delta, prop.Sigma$delta.corr/sd.ratio.delta, 1/sd.ratio.delta^2), 2, 2)
    
    
    # ---------------- Adapt sigma.m$X.s ------------------
    r.hat.X.s <- Accepted/thin
    sigma.m$X.s <- exp(log(sigma.m$X.s) + gamma1*(r.hat.X.s - metr.opt.1d))
    
    
    
    # ------------------------- Fill in trace objects ---------------------------
    X.s.trace[i, , ] <- X.s[,subset.replicates]
    X.trace[i, , ] <- X[,subset.replicates]
    X.s.accept.trace[i, ] <- X.s.accept
    rho.trace[i,] <- c(rho,nu)
    tau.trace[i] <- tau
    delta.trace[i,] <- c(gamma, sigma)
    theta.gpd.trace[i, ] <- theta.gpd[c(2,3)]
    shape.trace[i] <- shape
    # prob.below.trace[i] <- prob.below
    R.trace[i, ] <- R
    sd.ratio.trace[i] <- sd.ratio
    
    if(holdout) {
      X.s.holdout<-foreach(t = 1:n.t, .combine='cbind')%dopar%{
        X.s.to.Z <- qnorm(1-R[t]/X.s[,t])
        Z.tmp<-condMVNorm::rcmvnorm(1, mean=rep(0,n.s+n.holdout), sigma = Sigma_full, dependent.ind = (n.s+1):(n.s+n.holdout), given.ind = 1:n.s, X.given = X.s.to.Z)
        R[t]/(1-pnorm(Z.tmp))
      }
      X.s.holdout.trace[i-1, , ] <- X.s.holdout
    }
    
    
    # ----------------------------- Echo progress --------------------------------
    cat("Done with", i, "updates,\n")
    if ((i %% echo.interval) == 0) {
      cat(" Acceptance rate for X^* is", mean(X.s.accept.trace[(i-echo.interval):i, ]), "\n")
      pdf(file=sprintf("%s_progress.pdf", experiment.name))
      par(mfrow=c(3,3))
      plot(rho.trace[,1], type="l", ylab=expression(rho))
      if (!is.null(true.params)) abline(h=true.params$rho, lty=2, col=2, lwd=3)
      plot(rho.trace[,2], type="l", ylab=expression(nu))
      if (!is.null(true.params)) abline(h=true.params$nu, lty=2, col=2, lwd=3)
      plot(tau.trace, type="l", ylab=expression(tau^2))
      if (!is.null(true.params)) abline(h=true.params$tau, lty=2, col=2, lwd=3)
      plot(delta.trace[,1], type="l", ylab=expression(gamma))
      if (!is.null(true.params)) abline(h=true.params$gamma, lty=2, col=2, lwd=3)
      plot(delta.trace[,2], type="l", ylab=expression(sigma))
      if (!is.null(true.params)) abline(h=true.params$sigma, lty=2, col=2, lwd=3)
      # plot(prob.below.trace, type="l", ylab=expression(pi))
      # if (!is.null(true.params)) abline(h=true.params$prob.below, lty=2, col=2, lwd=3)
      plot(theta.gpd.trace[ ,1], type="l", ylab="a")
      if (!is.null(true.params)) abline(h=true.params$theta.gpd[2], lty=2, col=2, lwd=3)
      plot(theta.gpd.trace[ ,2], type="l", ylab="b")
      if (!is.null(true.params)) abline(h=true.params$theta.gpd[3], lty=2, col=2, lwd=3)
      plot(shape.trace, type="l", ylab="shape")
      if (!is.null(true.params)) abline(h=true.params$shape, lty=2, col=2, lwd=3)
      par(mfrow=c(3, 3))
      for (j in 1:min(18, n.t)) {
        plot(R.trace[ , j], type="l", main=paste("Replication", j), 
             ylab=paste0("R[", j, "]"),  ylim=range(R.trace[ , 1:min(18, n.t)], na.rm=TRUE))
        if (!is.null(true.params)) abline(h=true.params$R[j], lty=2, lwd=3, col="gray80")
      }
      # For each year, plot the location with an exceedance that happens to be first in the matrix
      for (j in 1:min(9, n.t)) { 
        plot.loc <- which(!cen[ ,subset.replicates[j]])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          if (!is.null(true.params)) abline(h=true.params$X.s[plot.loc, j], lty=2, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      # For each year, plot the location with a non-exceedance that happens to be first in the matrix
      for (j in 1:min(9, n.t)) { 
        plot.loc <- which(cen[ ,subset.replicates[j]])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylim=c(min(X.s.trace[ , plot.loc, j], na.rm=TRUE), thresh.X),
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          abline(h=thresh.X, lty=1, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      if(holdout) {
        # ylim<-range(X.s.holdout.trace[ , 1, 3],na.rm=TRUE)
        # if(ylim[2]>1e7) ylim[2]<-1e7
        # if(ylim[1]<-1e7) ylim[1]<--1e7
        ylim=c(-100, 10000)
        plot(X.s.holdout.trace[ , 1, 3], type="l", ylab=paste0("X*holdout[n1,t3]"),ylim=ylim)} 
      dev.off()
      
      # lines(curr.draw, col=alpha("gray80", 0.05))
      # # points(X.imp[ ,j], col=2, pch=20, cex=0.5)
      # Sys.sleep(0.1)
      state <- list(Y=Y, cen=cen, S=S, lon=lon, thresh=thresh,
                    experiment.name="Huser-wadsworth",
                    i=i, sigma.m=sigma.m, prop.Sigma=prop.Sigma, sd.ratio.rho=sd.ratio.rho,sd.ratio.trace=sd.ratio.trace,
                    X=X, X.s=X.s, theta.gpd=theta.gpd, shape=shape, prob.below=prob.below,
                    gamma=gamma, sigma=sigma, R=R, tau=tau, rho=rho, nu =nu, missing = missing)
      out.obj <- list(Y=Y, cen=cen, S=S, thresh=thresh,
                      i=i, missing = missing,
                      X.s.trace=X.s.trace,
                      X.trace=X.trace,
                      X.s.accept.trace=X.s.accept.trace,
                      rho.trace=rho.trace,
                      tau.trace=tau.trace,
                      delta.trace=delta.trace,
                      theta.gpd.trace=theta.gpd.trace,
                      shape.trace=shape.trace,
                      # prob.below.trace=prob.below.trace,
                      R.trace=R.trace,
                      sd.ratio.trace=sd.ratio.trace)
      if(holdout) out.obj<-c(out.obj,list(X.s.holdout.trace=X.s.holdout.trace))
      save(state, out.obj, file=sprintf("%s_progress_%1.1s.RData", experiment.name, save.bit))
      save.bit <- !save.bit
    }
  }
  
  return(out.obj)
  
}



scalemix.sampler.02.cont.DA <- function(Y, S, cen, lon, thresh,
                                        initial.values , i_prev=0,
                                        n.updates, thin=10,
                                        experiment.name="Huser-wadsworth",
                                        echo.interval=50,
                                        sigma.m=NULL, prop.Sigma=NULL, missing = NULL,
                                        true.params=NULL, sd.ratio=NULL, sd.ratio.rho=NULL, sd.ratio.delta=NULL, lower.prob.lim=0.5, out.obj,
                                        holdout=FALSE, S_full = NULL, n.holdout = NULL) {
  
  #library(doParallel)
  #library(foreach)
  save.bit <- TRUE
  
  # Constants to control how many Metropolis steps get executed at each
  # iteration of the main loop
  n.metr.updates.delta <- 4
  n.metr.updates.rho <- 4
  n.metr.updates.theta.gpd <- 4
  n.metr.updates.shape <- 2
  n.metr.updates.prob.below <- 4
  n.metr.updates.tau <- 2
  n.metr.updates.R <- 2
  
  # Constants to control adaptation of the Metropolis sampler
  c.0 <- 10
  c.1 <- 0.8
  k <- 3  # the iteration offset
  metr.opt.1d <- 0.41
  metr.opt.2d <- 0.35
  
  # Hyper parameters for the prior of the mixing distribution parameters and 
  # the correlation parameters
  hyper.params.delta <- 1
  hyper.params.rho <- 1
  hyper.params.theta.gpd <- max(lon[!cen])
  hyper.params.shape <- 1
  hyper.params.tau <- c(0.1,0.1)
  # hyper.params.prob.below <- c(0, 10)
  
  # A small number
  eps <- 1e-06
  
  # Bookkeeping
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  h <- rdist(S)
  diag(h)  <- 0
  if(holdout) {h_full <- rdist(S_full); diag(h_full) <- 0}
  #registerDoParallel(cores=n.t)
  subset.replicates <- ceiling(seq(1,n.t,length.out = 9))
  if(is.null(missing)) missing <- matrix(FALSE, nrow=n.s, ncol=n.t)
  
  # Load current values
  gamma <- initial.values$gamma
  sigma <- initial.values$sigma
  rho <- initial.values$rho
  nu <- initial.values$nu
  tau <- initial.values$tau
  a <- initial.values$theta.gpd[2]
  b <- initial.values$theta.gpd[3]
  shape <- initial.values$shape
  prob.below <- initial.values$prob.below
  # logit.prob.below <- logit(prob.below, c(lower.prob.lim, 1))
  X.s <- initial.values$X.s
  R <- initial.values$R
  
  X.s.accept <- rep(1, n.t)
  
  theta.gpd <- c(thresh, a, b)
  
  
  #**********   Eigendecomposition of the correlation matrix   ****************#
  Sigma <- corr.fn(h, rho, nu)
  eig.Sigma <- eigen(Sigma, symmetric=TRUE)
  V <- eig.Sigma$vectors
  d <- eig.Sigma$values
  if(holdout) Sigma_full <- corr.fn(h_full, rho, nu)
  
  # Initialize X and X.s
  thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, gamma = gamma, sigma = sigma)
  X <- Y
  X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd = tau, gamma = gamma, sigma=sigma, theta.gpd = theta.gpd, shape=shape, lon=lon[!cen & !missing], prob.below=prob.below)
  X[cen] <- update.censored.obs.mixture.me(X.s, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
  
  # Initialize trace objects
  accepted <- rep(0, n.updates)
  
  X.s.trace <- out.obj$X.s.trace
  X.trace <- out.obj$X.trace
  X.s.accept.trace <- out.obj$X.s.accept.trace
  rho.trace <- out.obj$rho.trace
  tau.trace <- out.obj$tau.trace
  delta.trace <- out.obj$delta.trace
  theta.gpd.trace <- out.obj$theta.gpd.trace
  shape.trace <- out.obj$shape.trace
  # prob.below.trace <- rep(NA, n.updates)
  R.trace <- out.obj$R.trace
  sd.ratio.trace <- out.obj$sd.ratio.trace
  X.s.holdout.trace <- out.obj$X.s.holdout.trace
  
  if(is.null(prop.Sigma$gpd.corr))  prop.Sigma$gpd.corr<-0
  if(is.null(prop.Sigma$rho.corr))  prop.Sigma$rho.corr<-0
  if(is.null(prop.Sigma$delta.corr))  prop.Sigma$delta.corr<-0
  
  # For tuning Metropolis updates of theta
  if (is.null(sigma.m$delta)) sigma.m$delta <- (2.4/2)^2
  if (is.null(sigma.m$rho)) sigma.m$rho <- (2.4/2)^2
  if (is.null(sigma.m$theta.gpd)) sigma.m$theta.gpd <- (2.4/2)^2
  if (is.null(sigma.m$shape)) sigma.m$shape <- 2.4^2
  # if (is.null(sigma.m$prob.below)) sigma.m$prob.below <- 1
  if (is.null(sigma.m$tau)) sigma.m$tau <- 2.4^2
  if (is.null(sigma.m$R)) sigma.m$R <- rep(2.4^2, n.t)
  if (is.null(sigma.m$X.s)) sigma.m$X.s <- matrix(2.4^2,nrow=n.s,ncol=n.t)
  
  r.hat.delta <- NA
  r.hat.rho <- NA
  r.hat.theta.gpd <- NA
  r.hat.shape <- NA
  r.hat.prob.below <- NA
  r.hat.tau <- NA
  # r.hat.R <- rep(NA, n.t)
  
  if (is.null(sd.ratio))   sd.ratio <- 7
  if (is.null(sd.ratio.rho))   sd.ratio.rho <- 7
  if (is.null(sd.ratio.delta))   sd.ratio.delta <- 7
  
  for (i in 1:n.updates) {
    
    ################################################################
    ## Update Metropolis adaptation parameters
    ################################################################
    gamma1 <- c.0 / (i + i_prev + k)^(c.1)
    gamma2 <- 1 / (i + i_prev + k)^(c.1)
    Accepted<-matrix(0, nrow=n.s, ncol=n.t) #For X.s
    
    for (j in 1:thin) {
      
      ################################################################
      ## Update X -- branching: use X.s
      ################################################################
      # X[cen] <- update.censored.obs.mixture.me(X.s = X.s, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
      X[!cen & !missing] <- gpd.2.scalemix.me(Y[!cen & !missing], tau_sqd = tau, gamma=gamma, sigma=sigma, theta.gpd = theta.gpd, shape = shape, lon = lon[!cen& !missing], prob.below=prob.below)
      
      ################################################################
      ## Update X.star  *Parallel
      ################################################################
      # X.s.res <- X.s.update.mixture.me.par(R, Y, X, X.s, cen, 
      #                                      prob.below, theta.gpd, delta,
      #                                      tau, V, d, v.q=2, n.chain = 100, thresh.X = thresh.X)
      X.s.res<-X.s.update.mixture.me.update.par.once.without.X.par(R, Y, X, X.s, cen,
                                                                   prob.below, theta.gpd, shape, lon, gamma,sigma,
                                                                   tau, V, d, v.q=sigma.m$X.s, thresh.X=thresh.X, missing = missing)
      X.s <- X.s.res$X.s
      X.s.accept <- apply(X.s.res$accepted==1,2,any)
      Accepted <- tryCatch(Accepted + X.s.res$accepted, error=function(e){cat("The dim of X.s.res$accepted is (",dim(X.s.res$accepted),")","\n")})
      
      
      
      ################################################################
      ## Update rho
      ################################################################
      
      metr.out.rho <- static.metr(z = R, starting.theta = c(rho, nu),
                                  likelihood.fn = rho.update.mixture.me.likelihood, prior.fn = rho.prior,
                                  hyper.params = hyper.params.rho, n.updates = n.metr.updates.rho, prop.Sigma = prop.Sigma$rho, sigma.m=sigma.m$rho, verbose=FALSE,
                                  X.s = X.s, R = R, S = S)
      r.hat.rho <- metr.out.rho$acc.prob
      tmp <- metr.out.rho$trace[n.metr.updates.rho,]
      rho <- tmp[1]; nu <- tmp[2]
      sigma.m$rho <- exp(log(sigma.m$rho) + gamma1*(r.hat.rho - metr.opt.2d))
      
      ## Re-create covariance matrix and eigenvectors/eigenvalues
      Sigma   <- corr.fn(h, rho, nu)
      if(holdout) Sigma_full <- corr.fn(h_full, rho, nu)
      eig.Sigma <- eigen(Sigma, symmetric=TRUE)
      V <- eig.Sigma$vectors
      d <- eig.Sigma$values
      
      
      
      ################################################################
      ## Update R  *Parallel
      ################################################################
      Metr.R<-foreach(t = 1:n.t, .combine = "rbind") %dopar% {
        
        metr.out.R <- static.metr(z = R, starting.theta = R[t],
                                  likelihood.fn = Rt.update.mixture.me.likelihood, prior.fn = huser.wadsworth.prior,
                                  hyper.params = c(gamma,sigma), n.updates = n.metr.updates.R, prop.Sigma = 1, sigma.m=sigma.m$R[t], verbose=FALSE,
                                  X.s = X.s[,t],  gamma=gamma, sigma=sigma, V = V, d = d)
        c(metr.out.R$trace[n.metr.updates.R],
          exp(log(sigma.m$R[t]) + gamma1*(metr.out.R$acc.prob - metr.opt.1d)))
      }
      
      R <- Metr.R[, 1]
      sigma.m$R <- Metr.R[, 2]
      
      
      
      ################################################################
      ## Update tau
      ################################################################
      metr.out.tau <-  static.metr(z = R, starting.theta = tau, 
                                   likelihood.fn = tau.update.mixture.me.likelihood, prior.fn = tau.sqd.prior,
                                   hyper.params = hyper.params.tau, n.updates = n.metr.updates.tau, prop.Sigma = 1, sigma.m=sigma.m$tau, verbose=FALSE, 
                                   Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, gamma=gamma, sigma=sigma, theta.gpd = theta.gpd, shape = shape, lon = lon,
                                   missing = missing)
      tau <- metr.out.tau$trace[n.metr.updates.tau]
      r.hat.tau <- metr.out.tau$acc.prob
      sigma.m$tau <- exp(log(sigma.m$tau) + gamma1*(r.hat.tau - metr.opt.1d))
      
      
      
      ################################################################
      ## Update delta
      ################################################################
      metr.out.delta <- static.metr(z = R, starting.theta = c(gamma, sigma), 
                                    likelihood.fn = delta.update.mixture.me.likelihood, prior.fn = interval.unif,
                                    hyper.params = hyper.params.delta, n.updates = n.metr.updates.delta, prop.Sigma = prop.Sigma$delta, sigma.m=sigma.m$delta, verbose=FALSE, 
                                    Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, theta.gpd = theta.gpd, shape = shape, lon = lon, tau_sqd = tau,
                                    missing = missing)
      tmp <- metr.out.delta$trace[n.metr.updates.delta,]
      gamma <- tmp[1]; sigma <- tmp[2]
      r.hat.delta <- metr.out.delta$acc.prob
      sigma.m$delta <- exp(log(sigma.m$delta) + gamma1*(r.hat.delta - metr.opt.2d))      
      
      
      ################################################################
      ## Update prob.below
      ################################################################
      #    metr.out.prob.below <- static.metr(R, logit.prob.below,
      #                                       logit.prob.below.update.mixture.me.likelihood,
      #                                       normal.scalar,
      #                                       hyper.params=hyper.params.prob.below,
      #                                       n.updates=n.metr.updates.prob.below,
      #                                       prop.Sigma=1,
      #                                       sigma.m=sigma.m$prob.below,
      #                                       Y=Y, X.s=X.s, cen=cen,
      #                                       theta.mix=theta.mix, theta.gaussian=theta.gaussian,
      #                                       theta.gpd=theta.gpd, lower.prob.lim=lower.prob.lim)
      #    logit.prob.below <- metr.out.prob.below$trace[n.metr.updates.prob.below]
      #    prob.below <- ilogit(logit.prob.below, c(lower.prob.lim, 1))
      #    r.hat.prob.below <- metr.out.prob.below$acc.prob
      #    sigma.m$prob.below <- exp(log(sigma.m$prob.below) +
      #                    gamma1*(r.hat.prob.below - metr.opt.1d))
      
      # Re-calculate the threshold on the X-scale    
      thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, gamma = gamma, sigma = sigma)
      
      
      
      ################################################################
      ## Update theta.gpd
      ################################################################
      metr.out.theta.gpd <- static.metr(z = R, starting.theta = theta.gpd[2:3],
                                        likelihood.fn = theta.gpd.update.mixture.me.likelihood, prior.fn = unif.a.b,
                                        hyper.params = hyper.params.theta.gpd, n.updates = n.metr.updates.theta.gpd, prop.Sigma = prop.Sigma$theta.gpd, sigma.m=sigma.m$theta.gpd, verbose=FALSE, 
                                        shape = shape, lon = lon, Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, gamma = gamma, sigma = sigma, tau_sqd = tau, loc = thresh, thresh.X=thresh.X, 
                                        missing = missing)
      
      theta.gpd[c(2,3)] <- metr.out.theta.gpd$trace[n.metr.updates.theta.gpd, ]
      r.hat.theta.gpd <- metr.out.theta.gpd$acc.prob
      sigma.m$theta.gpd <- exp(log(sigma.m$theta.gpd) + gamma1*(r.hat.theta.gpd - metr.opt.2d))
      
      ################################################################
      ## Update shape
      ################################################################
      metr.out.shape <- static.metr(z = R, starting.theta = shape, 
                                    likelihood.fn = shape.gpd.update.mixture.me.likelihood, prior.fn = unif.shape,
                                    hyper.params = hyper.params.shape, n.updates = n.metr.updates.shape, prop.Sigma = 1, sigma.m=sigma.m$shape, verbose=FALSE, 
                                    theta.gpd = theta.gpd, lon = lon, Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, gamma = gamma, sigma = sigma,  tau_sqd = tau, loc = thresh, 
                                    thresh.X=thresh.X, missing = missing)
      shape <- metr.out.shape$trace[n.metr.updates.shape]
      r.hat.shape <- metr.out.shape$acc.prob
      sigma.m$shape <- exp(log(sigma.m$shape) + gamma1*(r.hat.shape - metr.opt.1d))
      
    }
    
    
    
    # ---------------- Attempt to adapt proposal covariance ------------------
    if (r.hat.theta.gpd > 0) {
      sd.ratio.hat <- sd(metr.out.theta.gpd$trace[ ,1]) / sd(metr.out.theta.gpd$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio <- exp(log(sd.ratio) + gamma1*(log(sd.ratio.hat) - log(sd.ratio)))
    prop.Sigma$theta.gpd <-  matrix(c(1, prop.Sigma$gpd.corr/sd.ratio, prop.Sigma$gpd.corr/sd.ratio, 1/sd.ratio^2), 2, 2)
    
    if (r.hat.rho > 0) {
      sd.ratio.hat <- sd(metr.out.rho$trace[ ,1]) / sd(metr.out.rho$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio.rho <- exp(log(sd.ratio.rho) + gamma1*(log(sd.ratio.hat) - log(sd.ratio.rho)))
    prop.Sigma$rho <-  matrix(c(1, prop.Sigma$rho.corr/sd.ratio.rho, prop.Sigma$rho.corr/sd.ratio.rho, 1/sd.ratio.rho^2), 2, 2)
    
    if (r.hat.delta > 0) {
      sd.ratio.hat <- sd(metr.out.delta$trace[ ,1]) / sd(metr.out.delta$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio.delta <- exp(log(sd.ratio.delta) + gamma1*(log(sd.ratio.hat) - log(sd.ratio.delta)))
    prop.Sigma$delta <-  matrix(c(1, prop.Sigma$delta.corr/sd.ratio.delta, prop.Sigma$delta.corr/sd.ratio.delta, 1/sd.ratio.delta^2), 2, 2)
    
    
    # ---------------- Adapt sigma.m$X.s ------------------
    r.hat.X.s <- Accepted/thin
    sigma.m$X.s <- exp(log(sigma.m$X.s) + gamma1*(r.hat.X.s - metr.opt.1d))
    
    
    
    # ------------------------- Fill in trace objects ---------------------------
    X.s.trace[i + i_prev, , ] <- X.s[,subset.replicates]
    X.trace[i + i_prev, , ] <- X[,subset.replicates]
    X.s.accept.trace[i + i_prev, ] <- X.s.accept
    rho.trace[i + i_prev, ] <- c(rho,nu)
    tau.trace[i + i_prev] <- tau
    delta.trace[i + i_prev,] <- c(gamma, sigma)
    theta.gpd.trace[i + i_prev, ] <- theta.gpd[c(2,3)]
    shape.trace[i + i_prev] <- shape
    # prob.below.trace[i + i_prev] <- prob.below
    R.trace[i + i_prev, ] <- R
    sd.ratio.trace[i + i_prev] <- sd.ratio
    
    if(holdout) {
      X.s.holdout<-foreach(t = 1:n.t, .combine='cbind')%dopar%{
        X.s.to.Z <- qnorm(1-R[t]/X.s[,t])
        Z.tmp<-condMVNorm::rcmvnorm(1, mean=rep(0,n.s+n.holdout), sigma = Sigma_full, dependent.ind = (n.s+1):(n.s+n.holdout), given.ind = 1:n.s, X.given = X.s.to.Z)
        R[t]/(1-pnorm(Z.tmp))
      }
      X.s.holdout.trace[i-1+i_prev, , ] <- X.s.holdout
    }
    
    # ----------------------------- Echo progress --------------------------------
    cat("Done with", i, "updates,\n")
    if ((i %% echo.interval) == 0) {
      cat(" Acceptance rate for X^* is", mean(X.s.accept.trace[(i + i_prev-echo.interval):(i + i_prev), ]), "\n")
      pdf(file=sprintf("%s_progress.pdf", experiment.name))
      par(mfrow=c(3,3))
      plot(rho.trace[,1], type="l", ylab=expression(rho))
      if (!is.null(true.params)) abline(h=true.params$rho, lty=2, col=2, lwd=3)
      plot(rho.trace[,2], type="l", ylab=expression(nu))
      if (!is.null(true.params)) abline(h=true.params$nu, lty=2, col=2, lwd=3)
      plot(tau.trace, type="l", ylab=expression(tau^2))
      if (!is.null(true.params)) abline(h=true.params$tau, lty=2, col=2, lwd=3)
      plot(delta.trace[,1], type="l", ylab=expression(gamma))
      if (!is.null(true.params)) abline(h=true.params$gamma, lty=2, col=2, lwd=3)
      plot(delta.trace[,2], type="l", ylab=expression(sigma))
      if (!is.null(true.params)) abline(h=true.params$sigma, lty=2, col=2, lwd=3)
      # plot(prob.below.trace, type="l", ylab=expression(pi))
      # if (!is.null(true.params)) abline(h=true.params$prob.below, lty=2, col=2, lwd=3)
      plot(theta.gpd.trace[ ,1], type="l", ylab="a")
      if (!is.null(true.params)) abline(h=true.params$theta.gpd[2], lty=2, col=2, lwd=3)
      plot(theta.gpd.trace[ ,2], type="l", ylab="b")
      if (!is.null(true.params)) abline(h=true.params$theta.gpd[3], lty=2, col=2, lwd=3)
      plot(shape.trace, type="l", ylab="shape")
      if (!is.null(true.params)) abline(h=true.params$shape, lty=2, col=2, lwd=3)
      par(mfrow=c(3, 3))
      for (j in 1:min(18, n.t)) {
        plot(R.trace[ , j], type="l", main=paste("Replication", j), 
             ylab=paste0("R[", j, "]"),  ylim=range(R.trace[ , 1:min(18, n.t)], na.rm=TRUE))
        if (!is.null(true.params)) abline(h=true.params$R[j], lty=2, lwd=3, col="gray80")
      }
      # For each year, plot the location with an exceedance that happens to be first in the matrix
      for (j in 1:min(9, n.t)) { 
        plot.loc <- which(!cen[ ,subset.replicates[j]])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          if (!is.null(true.params)) abline(h=true.params$X.s[plot.loc, j], lty=2, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      # For each year, plot the location with a non-exceedance that happens to be first in the matrix
      for (j in 1:min(9, n.t)) { 
        plot.loc <- which(cen[ ,subset.replicates[j]])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylim=c(min(X.s.trace[ , plot.loc, j], na.rm=TRUE), thresh.X),
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          abline(h=thresh.X, lty=1, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      if(holdout) {
        # ylim<-range(X.s.holdout.trace[ , 1, 3],na.rm=TRUE)
        # if(ylim[2]>1e7) ylim[2]<-1e7
        # if(ylim[1]<-1e7) ylim[1]<--1e7
        ylim=c(-100, 10000)
        plot(X.s.holdout.trace[ , 1, 3], type="l", ylab=paste0("X*holdout[n1,t3]"),ylim=ylim)} 
      dev.off()
      
      
      state <- list(Y=Y, cen=cen, S=S, lon=lon, thresh=thresh,
                    experiment.name="Huser-wadsworth",
                    i=i + i_prev, sigma.m=sigma.m, prop.Sigma=prop.Sigma, sd.ratio.rho=sd.ratio.rho, sd.ratio.trace=sd.ratio.trace,
                    X=X, X.s=X.s, theta.gpd=theta.gpd, shape=shape, prob.below=prob.below,
                    gamma=gamma, sigma=sigma, R=R, tau=tau, rho=rho, nu=nu, missing = missing)
      out.obj <- list(Y=Y, cen=cen, S=S, thresh=thresh,
                      i=i + i_prev, missing = missing,
                      X.s.trace=X.s.trace,
                      X.trace=X.trace,
                      X.s.accept.trace=X.s.accept.trace,
                      rho.trace=rho.trace,
                      tau.trace=tau.trace,
                      delta.trace=delta.trace,
                      theta.gpd.trace=theta.gpd.trace,
                      shape.trace = shape.trace,
                      # prob.below.trace=prob.below.trace,
                      R.trace=R.trace,
                      sd.ratio.trace=sd.ratio.trace)
      if(holdout) out.obj<-c(out.obj,list(X.s.holdout.trace=X.s.holdout.trace))
      save(state, out.obj, file=sprintf("%s_progress_%1.1s.RData", experiment.name, save.bit))
      save.bit <- !save.bit
    }
  }
  
  return(out.obj)
}

