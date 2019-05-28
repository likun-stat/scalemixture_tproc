library(doParallel)
library(foreach)

source("~/Desktop/Research/scalemixture_tproc/scalemix_utils.R")
source("~/Desktop/Research/scalemixture_tproc/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture_tproc/scalemix_priors.R")
source("~/Desktop/Research/scalemixture_tproc/generic_samplers.R")
source("~/Desktop/Research/scalemixture_tproc/scalemix_sampler_02.R")

library(fields)   # For rdist
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

# ---------------- 1. Load the current states --------------------
tmp1<-file.mtime("Huser-wadsworth-DA_progress_F.RData")
tmp2<-file.mtime("Huser-wadsworth-DA_progress_T.RData")
Tdiff<-as.numeric(difftime(tmp1,tmp2), units="secs")
if(Tdiff>0) {load("Huser-wadsworth-DA_progress_F.RData")} else {load("Huser-wadsworth-DA_progress_T.RData")}

Y <- state$Y
cen <- state$cen
S <- state$S
thresh <- state$thresh
lon <- state$lon
missing <- state$missing
initial.values <- list(gamma=state$gamma, sigma=state$sigma,rho=state$rho, nu=state$nu, tau=state$tau, theta.gpd=state$theta.gpd, 
                       shape=state$shape, prob.below=state$prob.below, X.s=state$X.s, R=state$R)
i_prev <- state$i
n.updates<-15000

sigma.m <- state$sigma.m
prop.Sigma <- state$prop.Sigma
sd.ratio <- state$sd.ratio.trace[i_prev]
sd.ratio.rho <- state$sd.ratio.rho

# true.params <- list(delta = out.obj$delta.trace[1], rho=out.obj$rho.trace[1], tau=out.obj$tau.trace[1],
#                     theta.gpd=c(thresh, out.obj$theta.gpd.trace[1,]), shape=out.obj$shape.trace[1], prob.below=state$prob.below, X.s=out.obj$X.s.trace[1,,],
#                     R=out.obj$R.trace[1,])

rm(state)


## --------------- 2. Running Metropolis -------------------
scalemix.sampler.02.cont.DA(Y=Y, S=S, cen=cen, lon=lon, thresh=thresh,
                            initial.values=initial.values , i_prev=i_prev,
                            n.updates=n.updates, thin=10,
                            experiment.name="Huser-wadsworth-DA",
                            echo.interval=50,
                            sigma.m=sigma.m, prop.Sigma=prop.Sigma, missing=missing,
                            sd.ratio=sd.ratio, sd.ratio.rho = sd.ratio.rho, sd.ratio.delta = sd.ratio.delta, lower.prob.lim=0.5, out.obj=out.obj,
                            holdout=TRUE, S_full = ffwi.locs_full, n.holdout = 5)


stopCluster(cl)
