# required packages for implement
# =========================================================================================
library("foreach");library("truncnorm");library("pscl");library("igraph");library("MCMCpack");
library("Rcpp");library("RcppParallel");library("Brobdingnag")

# functions
# =========================================================================================
# The function to adjust sampling variance for target acceptance rate.
# =========================================================================================
# The arguments are:
#  1) accept: current acceptance rate.
#  2) sgm: current sampling variance.
#  3) target: target acceptance rate, default to 0.5.
# =========================================================================================
# The function returns the adjusted sampling variance sgm.

adjust_acceptance  <-  function(accept,sgm,target = 0.5){
  y = 1. + 1000.*(accept-target)*(accept-target)*(accept-target)
  if (y < .95)
    y = .95
  if (y > 1.05)
    y = 1.05
  sgm = sgm* y
  return(sgm)
}
# =========================================================================================
# The function to generate pre-specified kernel-based coordinates.
# =========================================================================================
# The arguments are:
#  1) pattern: the pre-specified spatialization pattern, where 'linear' is for Linear pattern, 'cos' is for Periodic pattern,
#      and 'gau' is for Exponential pattern.
#  2) q: the quantile value for scale parameter, default to 0.5 for the median value of the absolute coordinates.
#  3) loc: an n-by-2 matrix indicating the location of the corresponding spot center.
# =========================================================================================
# The function returns the adjusted sampling variance 
Kernel  <-  function(pattern = c("linear","gau","cos"),q=0.5,loc){
  S <- t(loc)-colMeans(loc)
  S <- t(S/apply(S, 1, sd))
  psi1 <- quantile(abs(S[,1]),q)
  psi2 <- quantile(abs(S[,2]),q)
  if(pattern == "gau"){
    b1 <- exp(-S[,1]^2/2/psi1^2)
    b2 <- exp(-S[,2]^2/2/psi2^2)
  }else if(pattern == "cos"){
    b1 <- cos(2*pi*S[,1]/psi1)
    b2 <- cos(2*pi*S[,2]/psi2)
  }else{
    b1 <- S[,1]
    b2 <- S[,2]
  }
  return(B = cbind(b1,b2))
}

# =========================================================================================
# The function to obtain target significant threshold based on Bayesian false discovery rate (BFDR) control strategy.
# =========================================================================================
# The arguments are:
#  1) PIP_vec: a posterior inclusion probability vector of p genes.
#  2) alpha: desired significance level, default to 0.05.
# =========================================================================================
# The function returns the value of the significant threshold for SV genes.
fdr_cut <- function(PIP_vec, alpha=0.05){
  p <- sort(1 - PIP_vec)
  thres <- cumsum(p)/seq_len(length(p))
  k <- which.max(thres >= alpha)
  cutoff = 1 - p[k-1]
  return(cutoff)
}

# =========================================================================================
# The function to generate newly sampled gamma parameter based on Langevin algorithm 
# =========================================================================================
# The arguments are:
#  1) gamma: current parameter gamma
#  2) tau_gamma: sampling standard deviation of gamma
#  3) sigma_gamma: standard deviation of prior normal distribution of gamma
#  4) L_vec: vector of first derivative of the log-likelihood function
#  5) eigmat_inv_List: list of the inverse of block-wise Laplacian matrix for parallelization 
#  6) vertices_List: list of the vertices of sub-networks for parallelization
# =========================================================================================
# The function returns the newly sampled values of the gamma parameter.
GenLangevin <- function(gamma, tau_gamma, sigma_gamma, L_vec, eigmat_inv_List, vertices_List){
  x  <-  foreach(t = eigmat_inv_List, .combine = 'c') %do% {
    mvnfast::rmvn(n = 1, mu = rep(0, ncol(t)), sigma = sigma_gamma*t)
  }
  y  <-  sqrt(1- tau_gamma^2) * gamma + (1-sqrt(1- tau_gamma^2))*L_vec + tau_gamma * x[order(unlist(vertices_List))] 
  return(y)
}

# =========================================================================================
# The function to generate the initial values for the corresponding component.
# =========================================================================================
# The arguments are:
#  1) Y: an n-by-p raw count matrix, where n is number of spots and p is the number of genes.
#  2) B: an n-by-2 specified matrix of spatial gene expression variation with pre-specified spatialization function,
#     Typical usage is to have the program compute the matrix based on arguments 'pattern', 'q', and 'loc'. Supplying the argument 'B' overrides this.
#  3) W: an n-by-K matrix with K being the number of cell type, where the ith row represents the cellular composition of the ith spot.
#  4) net: input network of the p genes.
#  5) nsim: the total number of the MCMC iterations, default to 2000. 
#  6) ntune: the number of the burn-in iterations, default to 1000. 
#  7) freqTune: the frequency of iterations for sampling variances tuning, default to 100.
#  8) nest: the number of the iterations for estimation, default to 1000.
#  9) nthin: the frequency of iterations for thinning, default to 10.
#  10) pattern: the pre-specified spatialization pattern, where 'linear' is for Linear pattern, 'cos' is for Periodic pattern,
#      and 'gau' is for Exponential pattern.
#  11) q: the quantile value for scale parameter, default to 0.5 for the median value of the absolute coordinates.
#  12) loc: an n-by-2 matrix indicating the location of the corresponding spot center.
#  13) vrbs: a boolean variable which defines whether to print the iteration details, default to TRUE.
# =========================================================================================
# The function returns a list containing the following components: 
#  1) beta1_hat: the initial spatial variation parameter beta1 at the first-axis. 
#  2) beta2_hat: the initial spatial variation parameter beta2 at the second-axis. 
#  3) alpha_hat: the initial alpha parameter for K cell types. 
#  4) phi_hat: the initial dispersion parameter phi. 
#  5) mu0_hat: the initial baseline expression parameter mu0. 

init <- function(Y,B=NULL,W,nsim=2000,ntune=1000,freqTune=100,nest=1000,nthin=10,pattern="linear",q=0.5,loc,vrbs=T){
  
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  K <- dim(W)[2]
  
  if(is.null(B)){
    B <-  Kernel(pattern=pattern,q=q,loc=loc)
  }
  tau.beta <- 0.05
  tau.phi <- 5
  tau.mu0 <- 0.01
  tau.alpha <- 0.05
  
  accept.beta <- accept.phi <- accept.mu0 <- accept.alpha <- 0
  
  mu0.lower <- -5
  mu0.upper <- 5
  mu0mean <- 0
  mu0sd <- 3
  alpha.lower <- -5
  alpha.upper <- 5
  alphamean <- 0
  alphasd <- 3
  
  s <- rep(1,n)
  mu0 <- rtruncnorm(p,a=mu0.lower,b=mu0.upper,mean=mu0mean,sd=mu0sd)
  alpha <- rtruncnorm(K,a=alpha.lower,b=alpha.upper,mean=alphamean,sd=alphasd)
  
  phi.lower <- 0
  phi.upper <- 100
  a.phi <- 10
  b.phi <- 0.1
  
  repeat{
    phi <- rgamma(1, shape = a.phi, scale = 1.0/b.phi)
    if(phi <- phi.upper) {break}
  }
  
  
  a.beta <- 2
  b.beta <- 0.5
  beta.ini <- rep(0,p)
  beta1 <- beta2 <- beta.ini
  burnin <- nsim-nest
  sigma.beta1 <- sigma.beta2 <- 1e-04
  
  a.pi <- b.pi <- 1
  
  numrow = nest/nthin
  
  Phi <- rep(0,numrow)
  Mu0 <- matrix(0,nrow=numrow,ncol=p)
  Beta1 <- matrix(0,nrow=numrow,ncol=p)
  Beta2 <- matrix(0,nrow=numrow,ncol=p)
  Alpha <- matrix(0,nrow=numrow,ncol=K)
  Rsum <- matrix(0,nrow=numrow,ncol=p)
  Sigma.beta1 <- rep(0,numrow)
  Sigma.beta2 <- rep(0,numrow)
  
  iter = seq((nsim-nest+1),nsim,by=nthin)
  
  sim <- 1
  R <- mu <- matrix(0,nrow=n,ncol=p)
  count <- 10
  tic <- Sys.time()
  
  if(vrbs){print('MCMC procedure for initialization starts!')}
  repeat{
    
    if(vrbs){
      if(sim*100/nsim == count) {
        print( paste0(count, "% has been done"))
        count <- count+10
      }
    }
    
    
    R <- update_Rmat(Y,B,W,beta1,beta2,mu0,alpha,s,phi,a.pi,b.pi)
    
    phi_results <- update_phi(R,Y,B,W,beta1,beta2,mu0,alpha,s,phi,phi.lower,phi.upper,
                            tau.phi,a.phi,b.phi,accept.phi)
    accept.phi <- phi_results$accept_phi
    phi <- phi_results$phi
    
    mu0.new  <-  rtruncnorm(p,mu0.lower,mu0.upper,mu0,sqrt(tau.mu0))
    mu0_results <- update_mu0(R,Y,B,W,beta1,beta2,mu0,mu0.new,alpha,s,phi,mu0mean,mu0sd,
                            mu0.lower,mu0.upper,tau.mu0,accept.mu0)
    accept.mu0 <- mu0_results$accept_mu0
    mu0 <- mu0_results$mu0 
    # if(sim<=ntune&sim%%500==0) cat("The acceptance of mu0 is", round(accept.mu0/(freqTune*p),3),"\n")
    
    alpha.new <- rtruncnorm(K,alpha.lower,alpha.upper,alpha,sqrt(tau.alpha))
    alpha_results <- update_alpha(R,Y,B,W,beta1,beta2,mu0,alpha,alpha.new,s,phi,alphamean,alphasd,
                                alpha.lower,alpha.upper,tau.alpha,accept.alpha)
    accept.alpha <- alpha_results$accept_alpha
    alpha <- alpha_results$alpha
    # if(sim<=ntune&sim%%500==0) cat("The acceptance of alpha is", round(accept.alpha/(freqTune*K),3),"\n")
    
    beta_results <- update_beta(R,Y,B,W,beta1,beta2,mu0,alpha,s,phi,sigma.beta1,sigma.beta2,
                              tau.beta,accept.beta)
    beta1  <-  beta_results$beta1
    beta2  <-  beta_results$beta2
    accept.beta  <-  beta_results$accept_beta
    # if(sim<=ntune&sim%%500==0) cat("The acceptance of beta is", accept.beta,"\n")
    
    sigma.beta1=rinvgamma(1,shape = a.beta+length(beta1)/2,as.numeric(b.beta+0.5*sum(beta1^2)))
    sigma.beta2=rinvgamma(1,shape = a.beta+length(beta2)/2,as.numeric(b.beta+0.5*sum(beta2^2)))
    
    
    if(sim<=ntune&sim%%freqTune==0){
      tau.beta <- adjust_acceptance(accept.beta/(freqTune),tau.beta,0.3)
      tau.phi <- adjust_acceptance(accept.phi/(freqTune),tau.phi,0.3) 
      tau.mu0 <- adjust_acceptance(accept.mu0/(freqTune*p),tau.mu0,0.3)
      tau.alpha <- adjust_acceptance(accept.alpha/(freqTune*K),tau.alpha,0.3)
      accept.beta <- 0
      accept.phi <- 0
      accept.mu0 <- 0
      accept.alpha <- 0
    }
    
    if(sim %in% iter){
      idx <- match(sim, iter)
      Beta1[idx,] <- beta1
      Beta2[idx,] <- beta2
      Alpha[idx,] <- alpha
      Phi[idx] <- phi
      Mu0[idx,] <- mu0
      Rsum[idx,] <- apply(R, MARGIN = 2,sum)
      Sigma.beta1[idx] <- sigma.beta1
      Sigma.beta2[idx] <- sigma.beta2
    }
    
    sim <- sim +1
    if(sim>nsim){break}
  }
  
  toc <- Sys.time()
  if(vrbs){
    print('MCMC procedure for initialization has been done :)')
    print(toc-tic)
  }
  
  rownames(Beta1) <- rownames(Beta2) <- paste0("iter", iter)
  colnames(Beta1) <- paste0('beta',1:p,'dim1')
  colnames(Beta2) <- paste0('gene',1:p,'dim2')
  beta1_hat <- colMeans(Beta1)
  beta2_hat <- colMeans(Beta2)
  phi_hat <- mean(Phi)
  mu0_hat <- colMeans(Mu0)
  alpha_hat <- colMeans(Alpha)
  
  return(list("beta1_hat"=beta1_hat,
              "beta2_hat"=beta2_hat,
              "alpha_hat"=alpha_hat,
              "phi_hat"=phi_hat,
              "mu0_hat"=mu0_hat))
}

