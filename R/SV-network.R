# The main function to conduct network-assisted SV genes identification
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
#  10) pars_ini: an optional list containing the initial values for the corresponding component. 
#      Initialization values can be supplied or automatically estimated (the default).
#  11) pattern: the pre-specified spatialization pattern, where 'linear' is for Linear pattern, 'cos' is for Periodic pattern,
#      and 'gau' is for Exponential pattern.
#  12) q: the quantile value for scale parameter, default to 0.5 for the median value of the absolute coordinates.
#  13) loc: an n-by-2 matrix indicating the location of the corresponding spot center.
#  14) a: a variable defines the desired significant level, default to 0.05.
#  15) vrbs: a boolean variable which defines whether to print the iteration details, default to TRUE.
# =========================================================================================
# The SV_network returns a list containing the following components: 
#  1) pars_ini: a list containing the initial values for the corresponding component.
#  2) save_mcmc:  a matrix containing the corresponding samples obtained after burn-in and thinning.
#  3) PIPs: a posterior inclusion probability vector of p genes.
#  4) loglik: the posterior model probability based on estimated parameters.
#  5) Pars: a list containing the estimated parameters for the corresponding component.
#  6) B: the n-by-2 specified matrix of spatial gene expression variation.
#  7) Y: the n-by-p raw count matrix.
#  8) net: the input network of the p genes.
#  9) gene_ids: the set of identified SV genes under the desired significance level 'alpha'.
#  10) tune_set: the vector of length 2 representing the lower and upper bound of the uniform distribution 
#      for the threshold parameter.


source("R/functions.R")
Rcpp::sourceCpp("R/mcmc_algo.cpp")
# =========================================================================================
SV_network  <-  function(Y,B=NULL,W,net=NULL,nsim=2000,ntune=1000,freqTune=100,nest=1000,nthin=10,
                 pars_ini=NULL,pattern,q=0.5,loc,a=0.05,vrbs=T){
  
  n  <-  dim(Y)[1]
  p <- dim(Y)[2]
  K <- dim(W)[2]
  
  a.gamma <- 3.5
  b.gamma <- 0.5
  a.phi <- 10
  b.phi <- 0.1
  a.pi <- b.pi <- 1
  mu0mean <- 0
  mu0sd <- 3
  alphamean <- 0
  alphasd <- 3
  
  phi.lower <- 0
  phi.upper <- 100
  mu0.lower <- -5
  mu0.upper <- 5
  alpha.lower <- -5
  alpha.upper <- 5
  
  tau.gamma <- 5e-04
  tau.lambda <- 0.01
  tau.phi <- 5
  tau.mu0 <- 0.01
  tau.alpha <- 0.05
  
  
  s <- rep(1,n)
  accept.gamma <- accept.lambda <- accept.phi <- accept.mu0 <- accept.alpha <- 0
  
  if(is.null(B)){
    B <- Kernel(pattern=pattern,q=q,loc=loc)
  }
  
  if(is.null(pars_ini)){
    pars_ini  <-  init(Y,B,W,nsim,ntune,freqTune,nest,nthin,vrbs=T)
  }
  beta1_ini <- pars_ini$beta1_hat
  beta2_ini <- pars_ini$beta2_hat
  phi <- pars_ini$phi_hat
  mu0 <- pars_ini$mu0_hat
  alpha <- pars_ini$alpha_hat
  
  lambda.lower <- min(0.5,as.numeric(quantile(abs(c(beta1_ini, beta2_ini)), 0.9)))
  lambda.upper <- min(1, as.numeric(quantile(abs(c(beta1_ini, beta2_ini)), 1)))
  
  if(!is.null(net)){
    degrees  <-  igraph::degree(net)
    vertice_0_name  <-  names(which(degrees==0))
    vertice_net_name  <-  names(which(degrees>0))
    vertice_0_id  <-  which(degrees==0)
    vertice_net_id  <-  which(degrees>0)
    if(length(vertice_0_name)){
      net  <-  net - vertice_0_name
    }
    L0  <-  laplacian_matrix(net,normalized=T,sparse=T)
    mat1  <-  sign(outer(beta1_ini[vertice_net_id], beta1_ini[vertice_net_id], "*"))
    mat2  <-  sign(outer(beta2_ini[vertice_net_id], beta2_ini[vertice_net_id], "*"))
    laplacian1  <-  bdiag(L0 * mat1)
    laplacian2  <-  bdiag(L0 * mat2)
    
    epsilon  <-  -10
    Ip  <-  base::diag(1,nrow(laplacian1))
    eigmat1  <-  laplacian1+exp(epsilon)*Ip
    eigmat1_inverse  <-  solve(eigmat1)
    eigmat2  <-  laplacian2+exp(epsilon)*Ip
    eigmat2_inverse  <-  solve(eigmat2)
    
    
    V(net)$name  <-  seq_len(vcount(net))
    subnets  <-  decompose(net, mode = "weak")
    subnets  <-  subnets[order(sapply(subnets,gorder),decreasing=F)]
    vertices_List  <-  eigmat1_inv_List  <-  eigmat2_inv_List  <-  Lbdcov_List  <-  list()
    for(i in 1:length(subnets)) {
      vertices_List[[i]]  <-  V(subnets[[i]])$name
      eigmat1_inv_List[[i]]  <-  eigmat1_inverse[vertices_List[[i]], vertices_List[[i]]]
      eigmat2_inv_List[[i]]  <-  eigmat2_inverse[vertices_List[[i]], vertices_List[[i]]]
    }
  }
  
  
  burnin <- nsim-nest
  gamma1 <- beta1_ini
  gamma2 <- beta2_ini
  
  if(!is.null(net)){
    wt1  <-  wt2  <-  c()
    temp1  <-  foreach(t = vertices_List, .combine = "c") %do%{
      rep(min((1.0/abs((gamma1[vertice_net_id][t])))^(1.0/2)),length(t))
    }
    
    temp2  <-  foreach(t = vertices_List, .combine = "c") %do%{
      rep(min((1.0/abs((gamma2[vertice_net_id][t])))^(1.0/2)),length(t))
    }
    
    wt1[vertice_net_id]  <-  temp1[order(unlist(vertices_List))] 
    wt2[vertice_net_id]  <-  temp2[order(unlist(vertices_List))]
    wt1[vertice_0_id]  <-  1.0/abs(gamma1[vertice_0_id])^(1.0/2)
    wt2[vertice_0_id]  <-  1.0/abs(gamma2[vertice_0_id])^(1.0/2)
  }
  
  
  
  lbd <- lambda.lower
  lambda1 <- lbd * wt1
  lambda2 <- lbd * wt2
  beta1 <- gamma1*as.numeric(abs(gamma1)>lambda1)
  beta2 <- gamma2*as.numeric(abs(gamma2)>lambda2)
  sigma.gamma1 <- sigma.gamma2 <- 1e-03
  sigma.lambda <- 0.01
  
  
  numrow <- nest/nthin
  
  Phi <- rep(0,numrow)
  Mu0 <- matrix(0, nrow = numrow, ncol = p)
  fgamma1 <- fgamma2 <- matrix(0,nrow=numrow,ncol=p)
  Beta1 <- Beta2 <- matrix(0,nrow=numrow,ncol=p)
  Alpha <- matrix(0,nrow=numrow,ncol=K)
  Gammas1 <- Gammas2 <- matrix(0,nrow=numrow,ncol=p)
  Rsum <- matrix(0,nrow=numrow,ncol=p)
  Lambda1 <- matrix(0, nrow=numrow, ncol = p)
  Lambda2 <- matrix(0, nrow=numrow, ncol = p)
  Sigma.gamma1 <- Sigma.gamma2 <- rep(0,numrow)
  Rs <- matrix(0,nrow=n,ncol=p)
  
  Pars <- list(
    beta1_hat=rep(NA,p),
    beta2_hat=rep(NA,p),
    beta_hat=rep(NA,p*ncol(B)),
    phi_hat=rep(NA,p),
    mu0_hat=rep(NA,p),
    alpha_hat=rep(NA,K),
    pi_hat=rep(NA,p),
    scale_hat=s,
    tags_hat=matrix(NA,n,p))
  
  iter <- seq((nsim-nest+1),nsim,by=nthin)
  
  sim <- 1
  R <- mu <- matrix(0,nrow=n,ncol=p)
  count <- 10
  tic <- Sys.time()
  
  
  repeat{
    if(vrbs){
      if(sim*100/nsim == count) {
        print(paste0(count, "% has been done"))
        count <- count + 10
      }
    }
    
    R <- update_Rmat(Y, B, W, beta1, beta2,mu0, alpha,s,phi,a.pi,b.pi)
    
    mu0.new <- rtruncnorm(p, a=mu0.lower, b=mu0.upper, mean=mu0, sd=sqrt(tau.mu0))
    mu0_results <- update_mu0(R, Y, B, W, beta1, beta2, mu0, mu0.new,
                              alpha, s, phi, mu0mean, mu0sd, 
                              mu0.lower, mu0.upper, tau.mu0, accept.mu0)
    accept.mu0 <- mu0_results$accept_mu0
    mu0 <- mu0_results$mu0 
    
    alpha.new <- rtruncnorm(K, a=alpha.lower, b=alpha.upper, mean=alpha, sd=sqrt(tau.alpha))
    alpha_results <- update_alpha(R, Y, B,W , beta1, beta2, mu0, 
                                  alpha, alpha.new, s, phi, alphamean, alphasd,
                                  alpha.lower, alpha.upper, tau.alpha, accept.alpha)
    accept.alpha <- alpha_results$accept_alpha
    alpha <- alpha_results$alpha
    
    Ls <- deltaUGamma(R, Y, B, W, mu0, 
                      alpha, s, phi, beta1, beta2, 
                      gamma1, gamma2, lambda1, lambda2, 
                      1e-03)
    if(!is.null(net)){
      L_vec1 <- as.vector(sigma.gamma1 * eigmat1_inverse %*% Ls$d1[vertice_net_id])
      L_vec2 <- as.vector(sigma.gamma2 * eigmat2_inverse %*% Ls$d2[vertice_net_id])
      gamma1.new <- gamma1
      gamma2.new <- gamma2
      gamma1.new[vertice_net_id] <- GenLangevin(gamma1[vertice_net_id],tau.gamma,sigma.gamma1,
                                              L_vec1,eigmat1_inv_List,vertices_List)
      gamma2.new[vertice_net_id] <-  GenLangevin(gamma2[vertice_net_id],tau.gamma,sigma.gamma2, 
                                               L_vec2,eigmat2_inv_List,vertices_List)
      gamma1.new[vertice_0_id] <- sqrt(1-tau.gamma^2)*gamma1[vertice_0_id]+
        (1-sqrt(1-tau.gamma^2))*sigma.gamma1*Ls$d1[vertice_0_id]+#/sqrt(n)+ 
        tau.gamma*rnorm(length(vertice_0_id),0, sqrt(sigma.gamma1))
      gamma2.new[vertice_0_id] <- sqrt(1-tau.gamma^2)*gamma2[vertice_0_id]+
        (1-sqrt(1-tau.gamma^2))*sigma.gamma2*Ls$d2[vertice_0_id]+#/sqrt(n)+ 
        tau.gamma*rnorm(length(vertice_0_id),0,sqrt(sigma.gamma2))
    }else{
      gamma1.new <- sqrt(1-tau.gamma^2)*gamma1+(1-sqrt(1-tau.gamma^2))*sigma.gamma1*Ls$d1 + 
        tau.gamma*rnorm(p,0,sqrt(sigma.gamma1))
      gamma2.new <- sqrt(1-tau.gamma^2)*gamma2+(1-sqrt(1-tau.gamma^2))*sigma.gamma2*Ls$d2 + 
        tau.gamma*rnorm(p,0,sqrt(sigma.gamma2))
    }
    
    gamma_results <- update_gamma(R, Y, B, W, beta1, beta2, 
                                  gamma1, gamma2, gamma1.new, gamma2.new, mu0,
                                  alpha, s, phi, lambda1, lambda2, accept.gamma)
    accept.gamma <- gamma_results$accept_gamma
    gamma1 <- gamma_results$gamma1
    gamma2 <- gamma_results$gamma2
    beta1 <- gamma_results$beta1
    beta2 <- gamma_results$beta2

    if(!is.null(net)){
      temp1  <-  foreach(t = vertices_List, .combine = "c") %do%{
        rep(min((1.0/abs((gamma1[vertice_net_id][t])))^(1.0/2)),length(t))
      }
      
      temp2  <-  foreach(t = vertices_List, .combine = "c") %do%{
        rep(min((1.0/abs((gamma2[vertice_net_id][t])))^(1.0/2)),length(t))
      }
      
      wt1[vertice_net_id]  <-  temp1[order(unlist(vertices_List))] 
      wt2[vertice_net_id]  <-  temp2[order(unlist(vertices_List))]
      wt1[vertice_0_id]  <-  1.0/abs(gamma1[vertice_0_id])^(1.0/2)
      wt2[vertice_0_id]  <-  1.0/abs(gamma2[vertice_0_id])^(1.0/2)
    }

    
    if(is.null(net)){
      sigma.gamma1 <- rinvgamma(1,shape = a.gamma+length(gamma1)/2,as.numeric(b.gamma+0.5*sum(gamma1^2)))
      sigma.gamma2 <- rinvgamma(1,shape = a.gamma+length(gamma2)/2,as.numeric(b.gamma+0.5*sum(gamma2^2)))
    }else{
      # sigma.gamma1 <- rinvgamma(1,shape = a.gamma+length(gamma1)/2,as.numeric(b.gamma+0.5*crossprod(gamma1, eigmat1%*%gamma1)))
      # sigma.gamma2 <- rinvgamma(1,shape = a.gamma+length(gamma2)/2,as.numeric(b.gamma+0.5*crossprod(gamma2, eigmat2%*%gamma2)))
      sigma.gamma1 <- rinvgamma(1,shape = a.gamma+length(gamma1)/2,as.numeric(b.gamma+0.5*crossprod(gamma1[vertice_net_id], eigmat1%*%gamma1[vertice_net_id]) + 
                                                                             0.5*sum(gamma1[vertice_0_id]^2)))
      sigma.gamma2 <- rinvgamma(1,shape = a.gamma+length(gamma2)/2,as.numeric(b.gamma+0.5*crossprod(gamma2[vertice_net_id], eigmat2%*%gamma2[vertice_net_id])+
                                                                             0.5*sum(gamma2[vertice_0_id]^2)))
    }
    
    lbd_new <- rnorm_trunc(lbd, sqrt(tau.lambda), lambda.lower, lambda.upper)
    lambda1 <- (lbd*wt1)
    lambda2 <- (lbd*wt2)
    beta1 <- gamma1*as.numeric(abs(gamma1)>(lambda1))
    beta2 <- gamma2*as.numeric(abs(gamma2)>(lambda2))
    beta1_new <- gamma1*as.numeric(abs(gamma1)>(lbd_new*wt1))
    beta2_new <- gamma2*as.numeric(abs(gamma2)>(lbd_new*wt2))
    
    temp <- mh_betaseq(R, Y, B, W, mu0, alpha, s, phi, beta1, beta2, beta1_new, beta2_new)
    log_mh <- temp$log_mh + 
      log(dtruncnorm(lbd,lambda.lower,lambda.upper,lbd_new, sqrt(tau.lambda))) -
      log(dtruncnorm(lbd_new,lambda.lower,lambda.upper,lbd,sqrt(tau.lambda)))

    if(log_mh>log(runif(1))){
      lbd <- lbd_new
      lambda1 <- lbd_new*wt1
      lambda2 <- lbd_new*wt2
      beta1 <- beta1_new
      beta2 <- beta2_new
      accept.lambda <- accept.lambda+1
    }
    
    
    if(sim<=ntune&sim%%freqTune==0){
      tau.gamma <- adjust_acceptance(accept.gamma/(freqTune),tau.gamma,0.3)
      tau.lambda <- adjust_acceptance(accept.lambda/(freqTune), tau.lambda, 0.3)
      tau.mu0 <- adjust_acceptance(accept.mu0/(freqTune*p),tau.mu0,0.3)
      tau.alpha <- adjust_acceptance(accept.alpha/(freqTune*K),tau.alpha,0.3)
      accept.gamma <- 0
      accept.lambda <- 0
      accept.mu0 <- 0
      accept.alpha <- 0
    }
    
    
    if(sim %in% iter){
      idx <- match(sim, iter)
      Gammas1[idx,] <- gamma1
      Gammas2[idx,] <- gamma2
      Sigma.gamma1[idx] <- sigma.gamma1
      Sigma.gamma2[idx] <- sigma.gamma2
      Lambda1[idx,] <- lambda1
      Lambda2[idx,] <- lambda2
      Beta1[idx,] <- beta1
      Beta2[idx,] <- beta2
      Alpha[idx,] <- alpha
      Phi[idx] <- phi
      Mu0[idx,] <- mu0
      Rsum[idx,] <- apply(R,2,sum)
      Rs <- Rs+R
      fgamma1[idx,] <- as.numeric(abs(gamma1)>lambda1)
      fgamma2[idx,] <- as.numeric(abs(gamma2)>lambda2)
    }
    
    sim <- sim+1
    if(sim>nsim){break}
  }
  
  toc <- Sys.time()
  print(toc-tic)
  
  PIPs <- apply(matrix(c(apply(fgamma1,2,mean),apply(fgamma2,2,mean)),nrow=p,ncol=2,byrow=F),1,max)
  prob_select <- fdr_cut(PIPs,alpha=a)
  ids <- which(PIPs > prob_select)
  
  
  beta  <- c(apply(Beta1,2,function(x) mean(x[x!=0])), apply(Beta2,2,function(x) mean(x[x!=0])))
  beta[is.na(beta)] <- 0
  tune_set <- c(lambda.lower,lambda.upper) 
  
  
  Pars$tags_hat <- Rs/numrow
  Pars$beta_hat <- beta
  Pars$beta1_hat <- beta[1:p]
  Pars$beta2_hat <- beta[(p+1):(2*p)]
  Pars$alpha_hat <- colMeans(Alpha)
  Pars$pi_hat <- colMeans(Rsum)/n
  Pars$mu0_hat <- colMeans(Mu0)
  Pars$phi_hat <- mean(Phi)
  
  loglik <- LoglikCpp(Y_mat=Y,B_mat=B, W_mat=W,scale_hat=Pars$scale_hat,mu0_hat=Pars$mu0_hat,beta_hat=Pars$beta_hat,
                    alpha_hat=Pars$alpha_hat,phi_hat=Pars$phi_hat,tags=Pars$tags_hat,thresh=0.8)
  
  save_mcmc <- cbind(Beta1,Beta2,Alpha,Gammas1,Gammas2,Sigma.gamma1,Sigma.gamma2,
                    Lambda1,Lambda2,Phi,Mu0,Rsum)
  rownames(save_mcmc) <- paste0("iter", iter)
  colnames(save_mcmc) <- c(paste0("beta",1:(p*ncol(B))), paste0("alpha",1:K), paste0("gamma",1:(p*ncol(B))),
                         'sigma.gamma1','sigma.gamm2',paste0("Lambda",1:(p*ncol(B))),'phi',paste0('mu0',1:p),paste0('Rsum',1:p))
  
  res <- list("pars_ini"= pars_ini,
              "save_mcmc"=save_mcmc,
              "PIPs"= PIPs,
              "loglik"=loglik,
              "Pars"=Pars,
              "B"=B,
              "count"=Y,
              "net"=net,
              "gene_ids"=sort(unique(ids)),
              "tune_set"=tune_set)
  
  return(res)
}



