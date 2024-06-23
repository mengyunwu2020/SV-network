// [[Rcpp::depends(RcppParallel,RcppArmadillo, BH)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <boost/random.hpp>

 using namespace Rcpp;
 using namespace RcppParallel;
 using namespace arma;

 double norm_rs(double a, double b);
 double half_norm_rs(double a, double b);
 double unif_rs(double a, double b);
 double exp_rs(double a, double b);
 double rnorm_trunc(double mu, double sigma, double lower, double upper);
 double dnorm_trunc(double x, double mu, double sigma, double lower, double upper);
 
 // [[Rcpp::export]]
 double dThresh(double gamma, double lambda, double epsilon) {
   double pi = 3.14159265;
   double temp = 2.0 / pi * gamma * epsilon / (epsilon * epsilon + pow((gamma * gamma - lambda * lambda), 2));
   return temp;
 }
 
 
 // [[Rcpp::export]]
 double Thresh(double gamma, double lambda, double epsilon){
   double pi = 3.14159265;
   double temp = 0.5  + 1.0/pi * atan((gamma * gamma - lambda * lambda)/epsilon);
   return temp;
 }
 
 // [[Rcpp::export]]
 double dnorm_trunc(double x, double mu, double sigma, double lower, double upper)
 {
   double dtemp, p_low, p_upp;
   dtemp = R::dnorm(x, mu, sigma, 0);
   p_low = R::pnorm(lower, mu, sigma, 1, 0);
   p_upp = R::pnorm(upper, mu, sigma, 1, 0);

   double temp = dtemp/(p_upp - p_low);
   return(temp);
 }
 
 
 // [[Rcpp::export]]
 double rnorm_trunc(double mu, double sigma, double lower, double upper)
 {
   int change;
   double a, b;
   double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
   double z, tmp, lograt;

   change = 0;
   a = (lower - mu)/sigma;
   b = (upper - mu)/sigma;

   // First scenario
   if( (a == R_NegInf)||(b == R_PosInf))
   {
     if(a == R_NegInf)
     {
       change = 1;
       a = -b;
       b = R_PosInf;
     }
     // The two possibilities for this scenario
     if(a <= 0.45) z = norm_rs(a, b);
     else z = exp_rs(a, b);
     if(change) z = -z;
   }

   // Second scenario
   else if((a*b) <= 0.0)
   {
     // The two possibilities for this scenario
     if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
     {
       z = norm_rs(a, b);
     }
     else z = unif_rs(a,b);
   }

   // Third scenario
   else
   {
     if(b < 0)
     {
       tmp = b; b = -a; a = -tmp; change = 1;
     }

     lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
     if(lograt <= logt2)
     {
       z = unif_rs(a,b);
     }
     else if((lograt > logt1)&&(a < t3))
     {
       z = half_norm_rs(a,b);
     }
     else
     {
       z = exp_rs(a,b);
     }
     if(change)
     {
       z = -z;
     }
   }
   double output;
   output = sigma*z + mu;
   return (output);
 }

 
 // [[Rcpp::export]]
 NumericVector mvmult( NumericMatrix m, NumericVector v, bool byrow=true)
 {
   if( m.ncol()!= v.length() ) stop("Non-conformable arrays") ;
   NumericVector out(m.nrow()) ;
   for (int i = 0; i < m.nrow(); i++) 
   {
     for (int j = 0; j < v.length(); j++) 
     {
       out[i] += m(i,j) * v[j];
     }
   }
   return out;
 }
 

 double norm_rs(double a, double b)
 {
   double x;
   x = Rf_rnorm(0.0, 1.0);
   while((x < a)||(x > b))
   {
     x = norm_rand();
   }
   return x;
 }

 double half_norm_rs(double a, double b)
 {
   double x;
   x = fabs(norm_rand());
   while((x<a)||(x>b))
   {
     x = fabs(norm_rand());
   }
   return x;
 }

 double exp_rs(double a, double b)
 {
   double  z, u, rate;
   rate = 1/a;

   // Generate a proposal on (0, b-a)
   z = R::rexp(rate);
   while(z > (b-a))
   {
     z = R::rexp(rate);
   }
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
     z = R::rexp(rate);
     while(z > (b-a))
     {
       z = R::rexp(rate);
     }
     u = R::runif(0.0,1.0);
   }
   return(z+a);
 }

 double unif_rs(double a, double b)
 {
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   if(a <= 0.0)
   {
     xstar = 0.0;
   }
   else
   {
     xstar = a;
   }
   logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

   x = R::runif(a, b);
   logu = log(R::runif(0.0, 1.0));
   while(logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
   {
     x = R::runif(a, b);
     logu = log(R::runif(0.0, 1.0));
   }
   return x;
 }


 
 

 struct UpdateRmatWorker : public Worker {
   RMatrix<double> result;
   const RMatrix<double> Y_mat;
   const RMatrix<double> B_mat;
   const RVector<double> beta1;
   const RVector<double> beta2;
   const RVector<double> mu0;
   const RVector<double> alphae;
   const RVector<double> s;
   const double phi;
   const double a_pi;
   const double b_pi;

   UpdateRmatWorker(NumericMatrix result, const NumericMatrix Y_mat, const NumericMatrix B_mat,
                    const NumericVector beta1, const NumericVector beta2, 
                    const NumericVector mu0, const NumericVector alphae, const NumericVector s, const double phi,
                    const double a_pi, const double b_pi)
     : result(result), Y_mat(Y_mat), B_mat(B_mat), beta1(beta1), beta2(beta2), 
       mu0(mu0), alphae(alphae), s(s), phi(phi), a_pi(a_pi), b_pi(b_pi) {}
   

   void operator()(std::size_t begin, std::size_t end) {
     boost::random::mt19937 rng;
     for (std::size_t i = begin; i < end; i++) {
       for (int j = 0; j < result.ncol(); j++) {
         if (Y_mat(i, j) == 0) {
           double phi_j = phi, mu0_j = mu0[j], s_i = s[i];
           double xb = B_mat(i,0)*beta1[j] + B_mat(i,1)*beta2[j] + alphae[i];
           double mu_ij = exp(mu0_j+xb);
           std::vector<double> lg_p_vec(2);
           lg_p_vec[1] = lgamma(a_pi + 1) + lgamma(b_pi) - lgamma(a_pi + 1 + b_pi);
           lg_p_vec[0] = phi_j * (log(phi_j) - log(s_i * mu_ij + phi_j)) +
             lgamma(a_pi) + lgamma(b_pi + 1) - lgamma(a_pi + 1 + b_pi);
           //double max_tmp = max(lg_p_vec);
           double max_tmp = lg_p_vec[0];
           if (lg_p_vec[1] > lg_p_vec[0]){
             max_tmp = lg_p_vec[1];
           }
           double sum_tmp = 0;
           for(int m = 0; m < 2; m++){
             lg_p_vec[m] = lg_p_vec[m] - max_tmp;
             lg_p_vec[m] = exp(lg_p_vec[m]);
             sum_tmp = sum_tmp + lg_p_vec[m];
           }
           for(int m = 0; m < 2; m++){
             lg_p_vec[m] = lg_p_vec[m] / sum_tmp;
           }

           boost::random::binomial_distribution<int> binomial(1, lg_p_vec[1]);
           int r_ij = binomial(rng);
             //rbinom(1, 1, lg_p_vec[1])[0];
           result(i,j) = r_ij;
         }
       }
     }
   }
 };


 // [[Rcpp::export]]
 NumericMatrix update_Rmat(NumericMatrix Y_mat, NumericMatrix B_mat, NumericMatrix W_mat,
                           NumericVector beta1, NumericVector beta2, NumericVector mu0, NumericVector alpha,
                           NumericVector s, double phi, double a_pi, double b_pi) {
   int n = Y_mat.nrow();
   int p = Y_mat.ncol();

   NumericMatrix result(n, p);
   NumericVector alphae =  mvmult(W_mat, alpha);

   UpdateRmatWorker worker(result, Y_mat, B_mat, beta1, beta2, mu0, alphae, s, phi, a_pi, b_pi);

   parallelFor(0, n, worker);

   return result;
 }


 struct MhPhiWorker : public Worker {
   const RMatrix<double> R_mat;
   const RMatrix<double> Y_mat;
   const RMatrix<double> B_mat;
   const RVector<double> beta1;
   const RVector<double> beta2;
   const RVector<double> mu0;
   const RVector<double> alphae;
   const RVector<double> s;
   const double phi;
   const double phi_new;
   double log_mh;

   MhPhiWorker(const NumericMatrix R_mat, const NumericMatrix Y_mat, const NumericMatrix B_mat,
               const NumericVector beta1, const NumericVector beta2, const NumericVector mu0, 
               const NumericVector alphae, const NumericVector s,
               const double phi, const double phi_new)
     : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat), beta1(beta1), beta2(beta2), mu0(mu0), alphae(alphae), s(s),
       phi(phi), phi_new(phi_new), log_mh() {}

   MhPhiWorker(const MhPhiWorker& mpw, Split)
     : R_mat(mpw.R_mat), Y_mat(mpw.Y_mat), B_mat(mpw.B_mat), beta1(mpw.beta1), beta2(mpw.beta2), mu0(mpw.mu0), 
       alphae(mpw.alphae), s(mpw.s), phi(mpw.phi), phi_new(mpw.phi_new), log_mh() {}

   void operator()(std::size_t begin, std::size_t end) {

     for (int j = begin; j < end; j++) {
       for (int i = 0; i < Y_mat.nrow(); i++) {
         if (R_mat(i, j) == 0) {
           double s_i = s[i], y_ij = Y_mat(i, j), mu0_j = mu0[j];
           double xb = B_mat(i, 0) * beta1[j] + B_mat(i, 1) * beta2[j] + alphae[i];
           double mu_ij = exp(mu0_j + xb);

           log_mh +=  lgamma(y_ij + phi_new) + lgamma(phi) - lgamma(y_ij + phi) - lgamma(phi_new) +
             phi_new * (log(phi_new) - log(s_i * mu_ij + phi_new)) -
             phi * (log(phi) - log(s_i * mu_ij + phi)) +
             y_ij * (log(s_i * mu_ij + phi) - log(s_i * mu_ij + phi_new));
         }
       }
     }
   }

   void join(const MhPhiWorker& other) {
     log_mh += other.log_mh;
   }
 };

 // [[Rcpp::export]]
 double mh_phi(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat,
                   NumericVector beta1, NumericVector beta2, NumericVector mu0, NumericVector alphae, 
                   NumericVector s, double phi, double phi_new) {
   int p = Y_mat.ncol();

   MhPhiWorker worker(R_mat, Y_mat, B_mat, beta1, beta2, mu0, alphae, s, phi, phi_new);

   parallelReduce(0, p, worker);

   return worker.log_mh;
 }

 // [[Rcpp::export]]
 List update_phi(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat, NumericMatrix W_mat,
                     NumericVector beta1, NumericVector beta2, NumericVector mu0, NumericVector alpha, NumericVector s, double phi,
                     double phi_lower, double phi_upper, double tau_phi,
                     double a_phi, double b_phi, int accept_phi) {

   NumericVector alphae =  mvmult(W_mat, alpha);
   
   double phi_new = rnorm_trunc(phi, sqrt(tau_phi), phi_lower, phi_upper);
   double log_phi_prior = (a_phi-1) * (log(phi_new) - log(phi)) + b_phi*(phi- phi_new);

   double log_mh = mh_phi(R_mat, Y_mat, B_mat, beta1, beta2, mu0, alphae, s, phi, phi_new);
   double new_dens = dnorm_trunc(phi_new, phi, sqrt(tau_phi), phi_lower, phi_upper );
   double old_dens = dnorm_trunc(phi, phi_new, sqrt(tau_phi), phi_lower, phi_upper );

   log_mh = log_mh + log_phi_prior + log(old_dens) - log(new_dens);

   if(log_mh > log(unif_rand())){
     phi = phi_new;
     accept_phi = accept_phi + 1;
   }

   List result;
   result["phi"] = phi;
   result["accept_phi"] = accept_phi;
   result["log_mh"] = log_mh;
   result["phi_new"] = phi_new;
   result["new_dens"] = new_dens;
   result["log_phi_prior"] = log_phi_prior;
   result["old_dens"] = old_dens;
   result["alphae"] = alphae;

   return result;

 }

 
 struct Mu0Update : public RcppParallel::Worker {
   const  NumericMatrix R_mat;
   const  NumericMatrix Y_mat;
   const  NumericMatrix B_mat;
   const  NumericVector beta1;
   const  NumericVector beta2;
   const  NumericVector alphae;
   const  NumericVector s;
   const double phi;
   const double mu0_mean;
   const double mu0_sd;
   const double mu0_lower;
   const double mu0_upper;
   const double tau_mu0;
   NumericVector mu0;
   const NumericVector mu0_new;
   NumericVector logmh;
   NumericVector acceptmu0;
   
   Mu0Update(const  NumericMatrix R_mat, const  NumericMatrix Y_mat, const  NumericMatrix B_mat,
             const  NumericVector beta1, const  NumericVector beta2,  const  NumericVector alphae, 
             const  NumericVector s, const double phi, const double mu0_mean, const double mu0_sd,const double mu0_lower, 
             const double mu0_upper, const double tau_mu0,  NumericVector mu0, NumericVector mu0_new,
             NumericVector logmh,  NumericVector acceptmu0)
     : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat), beta1(beta1), beta2(beta2), alphae(alphae), s(s), phi(phi),
       mu0_mean(mu0_mean), mu0_sd(mu0_sd), mu0_lower(mu0_lower), mu0_upper(mu0_upper),tau_mu0(tau_mu0), mu0(mu0), mu0_new(mu0_new), 
       logmh(logmh), acceptmu0(acceptmu0) {}
   
   
   void operator()(std::size_t begin, std::size_t end) {
     int n = Y_mat.nrow();
     
     for (int j = begin; j < end; j++) {
       double phi_j = phi;
       double mu0_old = mu0[j];
       double mu0_star = mu0_new[j];
       double log_mu0_prior = 0.5 * ((mu0_old - mu0_mean) *(mu0_old - mu0_mean) - 
                                     (mu0_star - mu0_mean)*( mu0_star - mu0_mean)) / (mu0_sd * mu0_sd); 
       double log_mh = 0;
       for (int i = 0; i < n; i++) {
         if (R_mat(i, j) == 0) {
           double s_i = s[i], y_ij = Y_mat(i, j);
           double xb = B_mat(i, 0) * beta1[j] + B_mat(i, 1) * beta2[j] + alphae[i];
           double mu_star = exp(mu0_star + xb);
           double mu_old = exp(mu0_old + xb);
           log_mh = log_mh + phi_j * (log(s_i * mu_old + phi_j) - log(s_i * mu_star + phi_j)) +
             y_ij * (log(s_i * mu_old + phi_j) + log(mu_star) - log(s_i * mu_star + phi_j) - log(mu_old));
         }
       }
       
       double new_dens = dnorm_trunc(mu0_star, mu0_old, sqrt(tau_mu0),  mu0_lower, mu0_upper);
       double old_dens = dnorm_trunc(mu0_old, mu0_star, sqrt(tau_mu0),  mu0_lower, mu0_upper);
       
       log_mh = log_mh + log_mu0_prior + log(old_dens) - log(new_dens);
       logmh[j] = log_mh;
       
       if (log_mh > log(unif_rand())) {
         mu0[j] = mu0_star;
         acceptmu0[j] = 1;
       }else{
         mu0[j] = mu0_old;
       }
     }
   }
 };
 
 // [[Rcpp::export]]
 List update_mu0( NumericMatrix R_mat,  NumericMatrix Y_mat,  NumericMatrix B_mat, NumericMatrix W_mat,
                      NumericVector beta1,  NumericVector beta2, NumericVector mu0, NumericVector mu0_new,
                      NumericVector alpha, NumericVector s,  double phi, double mu0_mean,
                      double mu0_sd, double  mu0_lower, double mu0_upper, double tau_mu0, int accept_mu0) {
   int p = Y_mat.ncol();
   
   NumericVector alphae =  mvmult(W_mat, alpha);
   NumericVector logmh(p);
   NumericVector acceptmu0(p);
   
   Mu0Update mu0Update(R_mat, Y_mat, B_mat, beta1, beta2, alphae, s, phi, mu0_mean, mu0_sd, 
                       mu0_lower, mu0_upper, tau_mu0, mu0,mu0_new, logmh, acceptmu0);
   
   RcppParallel::parallelFor(0, p, mu0Update);
   
   accept_mu0 = accept_mu0 + sum(acceptmu0);
   List result;
   result["mu0"] = mu0;
   result["mu0_new"] = mu0_new;
   result["accept_mu0"] = accept_mu0;
   result["accept_mu0_vec"] = acceptmu0;
   result["log_mh"] = logmh;
   
   return result;
 }
 
 
 struct MhalphaWorker : public Worker {
   const RMatrix<double> R_mat;
   const RMatrix<double> Y_mat;
   const RMatrix<double> B_mat;
   const RVector<double> beta1;
   const RVector<double> beta2;
   const RVector<double> mu0;
   const RVector<double> alphae;
   const RVector<double> alphae_new;
   const RVector<double> s;
   const double phi;
   double log_mh;
   
   MhalphaWorker(const NumericMatrix R_mat, const NumericMatrix Y_mat, const NumericMatrix B_mat,
               const NumericVector beta1, const NumericVector beta2, const NumericVector mu0, 
               const NumericVector alphae, const NumericVector alphae_new, const NumericVector s,
               const double phi)
     : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat), beta1(beta1), beta2(beta2), mu0(mu0), alphae(alphae), 
       alphae_new(alphae_new), s(s), phi(phi), log_mh() {}
   
   MhalphaWorker(const MhalphaWorker & maw, Split)
     : R_mat(maw.R_mat), Y_mat(maw.Y_mat), B_mat(maw.B_mat), beta1(maw.beta1), beta2(maw.beta2), mu0(maw.mu0), 
       alphae(maw.alphae), alphae_new(maw.alphae_new), s(maw.s), phi(maw.phi), log_mh() {}
   
   void operator()(std::size_t begin, std::size_t end) {
     
     for (int j = begin; j < end; j++) {
       double mu0_j = mu0[j];
       double phi_j = phi;
       for (int i = 0; i < Y_mat.nrow(); i++) {
         if (R_mat(i, j) == 0) {
           double s_i = s[i], y_ij = Y_mat(i, j);
           double xb = B_mat(i, 0) * beta1[j] + B_mat(i, 1) * beta2[j] + alphae[i];
           double xb_star = B_mat(i, 0) * beta1[j] + B_mat(i, 1) * beta2[j] + alphae_new[i];
           double mu_old = exp(mu0_j + xb);
           double mu_star = exp(mu0_j + xb_star);
           
           log_mh +=  phi_j * (log(s_i* mu_old + phi_j) - log(s_i * mu_star + phi_j)) +
             y_ij * (log(s_i * mu_old + phi_j) + log(mu_star) - log(s_i * mu_star + phi_j) - log(mu_old));
           
         }
       }
     }
   }
   
   void join(const MhalphaWorker& other) {
     log_mh += other.log_mh;
   }
 };
 

 // [[Rcpp::export]]
 double mh_alpha(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat,
                   NumericVector beta1, NumericVector beta2, NumericVector mu0, NumericVector alphae, 
                   NumericVector alphae_new, NumericVector s, double phi) {
   int p = Y_mat.ncol();
   
   MhalphaWorker worker(R_mat, Y_mat, B_mat, beta1, beta2, mu0, alphae, alphae_new, s, phi);
   
   parallelReduce(0, p, worker);
   
   return worker.log_mh;
 }
 
 
 // [[Rcpp::export]]
 List update_alpha( NumericMatrix R_mat,  NumericMatrix Y_mat,  NumericMatrix B_mat, NumericMatrix W_mat,
                      NumericVector beta1,  NumericVector beta2, NumericVector mu0, NumericVector alpha, NumericVector alpha_new, 
                      NumericVector s,  double phi, double alpha_mean, double alpha_sd, 
                      double  alpha_lower, double alpha_upper, double tau_alpha, int accept_alpha) {
   int K = W_mat.ncol();

   NumericVector acceptalpha(K);
   NumericVector log_mh_vec(K);
   
   
   for (int k = 0; k < K; k++) {
     
     // Rcpp::Rcout << "The iteration is " << std::endl << k << std::endl;
     // Rcpp::Rcout << "alpha is " << std::endl << alpha << std::endl;
     // Rcpp::Rcout << "alpha_new is " << std::endl << alpha_new << std::endl;
     
     double alpha_old = alpha[k];
     // Rcpp::Rcout << "alpha old is " << std::endl << alpha_old << std::endl;
     double alpha_star = alpha_new[k];
     // Rcpp::Rcout << "alpha star is " << std::endl << alpha_star << std::endl;
     double log_alpha_prior = 0.5 * ((alpha_old - alpha_mean) * (alpha_old - alpha_mean) - 
                                     (alpha_star - alpha_mean) * (alpha_star - alpha_mean)) / (alpha_sd * alpha_sd);
     // Rcpp::Rcout << "alpha prior is " << std::endl <<  log_alpha_prior << std::endl;
     double new_dens = dnorm_trunc(alpha_star, alpha_old, sqrt(tau_alpha),  alpha_lower, alpha_upper);
     // Rcpp::Rcout << "new dens is " << std::endl <<  new_dens << std::endl;
     double old_dens = dnorm_trunc(alpha_old, alpha_star, sqrt(tau_alpha),  alpha_lower, alpha_upper);
     // Rcpp::Rcout << "old dens is " << std::endl <<  old_dens << std::endl;
     
     NumericVector alphae =  mvmult(W_mat, alpha);
     NumericVector alpha_new_vec = alpha;
     alpha_new_vec[k] = alpha_new[k];
     NumericVector alphae_new =  mvmult(W_mat, alpha_new_vec);
     
     double log_mh = mh_alpha(R_mat, Y_mat, B_mat, beta1, beta2, mu0, alphae, alphae_new, s, phi);
     
     log_mh = log_mh + log_alpha_prior + log(old_dens) - log(new_dens);
     log_mh_vec[k] = log_mh;
     
     if(log_mh > log(unif_rand())){
       // Rcpp::Rcout << "Yes" << std::endl << k << std::endl;
       alpha[k] = alpha_star;
       acceptalpha[k] = 1;
     }else{
       alpha[k] = alpha_old;
     }
     
   }
   
   accept_alpha = accept_alpha + sum(acceptalpha);
   List result;
   result["alpha"] = alpha;
   result["alpha_new"] = alpha_new;
   result["accept_alpha"] = accept_alpha;
   result["accept_alpha_vec"] = acceptalpha;
   result["log_mh"] = log_mh_vec;
   
   return result;
 }
 
 
 
 
 struct MhBetaWorker : public Worker {
   const RMatrix<double> R_mat;
   const RMatrix<double> Y_mat;
   const RMatrix<double> B_mat;
   const RVector<double> mu0;
   const RVector<double> alphae;
   const RVector<double> s;
   const double phi;
   RVector<double> beta1;
   RVector<double> beta2;
   const RVector<double> beta1_new;
   const RVector<double> beta2_new;
   double log_mh;
   
   MhBetaWorker(const NumericMatrix R_mat, const NumericMatrix Y_mat, const NumericMatrix B_mat,
                const NumericVector mu0, const NumericVector alphae, const NumericVector s, const double phi,
                NumericVector beta1, NumericVector beta2, const NumericVector beta1_new, const NumericVector beta2_new)
     : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat),  mu0(mu0), alphae(alphae), s(s), phi(phi),
       beta1(beta1), beta2(beta2), beta1_new(beta1_new), beta2_new(beta2_new), log_mh() {}
   
   MhBetaWorker(const MhBetaWorker& other, Split)
     : R_mat(other.R_mat), Y_mat(other.Y_mat), B_mat(other.B_mat), mu0(other.mu0), alphae(other.alphae), 
       s(other.s), phi(other.phi), beta1(other.beta1), beta2(other.beta2), 
       beta1_new(other.beta1_new), beta2_new(other.beta2_new), log_mh() {}
   
   void operator()(std::size_t begin, std::size_t end) {
     
     for (int j = begin; j < end; j++) {
       double phi_j = phi;
       double mu0_j = mu0[j];
       for (int i = 0; i < Y_mat.nrow(); i++) {
         if(R_mat(i,j) == 0){
           double s_i = s[i], y_ij = Y_mat(i,j);
           double xb = B_mat(i,0)* beta1[j] + B_mat(i,1)* beta2[j] + alphae[i];
           double xb_new = B_mat(i,0)* beta1_new[j]+B_mat(i,1)* beta2_new[j] + alphae[i];
           double mu_star = exp(mu0_j + xb_new);
           double mu_old  = exp(mu0_j + xb);
           log_mh = log_mh + phi_j * (log(s_i* mu_old + phi_j) - log(s_i * mu_star + phi_j)) +
             y_ij * (log(s_i * mu_old + phi_j) + log(mu_star) - log(s_i * mu_star + phi_j) - log(mu_old));
         }
       }
     }
   }
   
   void join(const MhBetaWorker& other) {
     log_mh += other.log_mh;
   }
 };

 
 
 // [[Rcpp::export]]
 double mh_beta(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat,
                    NumericVector mu0, NumericVector alphae, NumericVector s, double phi, NumericVector beta1, NumericVector beta2, 
                    NumericVector beta1_new, NumericVector beta2_new) {
   int p = Y_mat.ncol();
   
   MhBetaWorker worker(R_mat, Y_mat, B_mat, mu0, alphae, s, phi, beta1, beta2, beta1_new, beta2_new);
   
   parallelReduce(0, p, worker);
   
   return worker.log_mh;
 }
 
 
 struct MhBetaseqWorker : public Worker {
   const RMatrix<double> R_mat;
   const RMatrix<double> Y_mat;
   const RMatrix<double> B_mat;
   const RVector<double> mu0;
   const RVector<double> alphae;
   const RVector<double> s;
   const double phi;
   RVector<double> beta1;
   RVector<double> beta2;
   const RVector<double> beta1_new;
   const RVector<double> beta2_new;
   RVector<double> log_mh_vec;
   double log_mh;
   
   
   MhBetaseqWorker(const NumericMatrix R_mat, const NumericMatrix Y_mat, const NumericMatrix B_mat,
                   const NumericVector mu0, const NumericVector alphae, const NumericVector s, const double phi,
                   NumericVector beta1, NumericVector beta2, const NumericVector beta1_new, const NumericVector beta2_new, NumericVector log_mh_vec)
     : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat),  mu0(mu0), alphae(alphae), s(s), phi(phi),
       beta1(beta1), beta2(beta2), beta1_new(beta1_new), beta2_new(beta2_new), log_mh_vec(log_mh_vec), log_mh() {}
   
   MhBetaseqWorker(const MhBetaseqWorker& other, Split)
     : R_mat(other.R_mat), Y_mat(other.Y_mat), B_mat(other.B_mat), mu0(other.mu0), alphae(other.alphae), s(other.s), phi(other.phi),
       beta1(other.beta1), beta2(other.beta2), beta1_new(other.beta1_new), beta2_new(other.beta2_new), log_mh_vec(other.log_mh_vec), log_mh() {}
   
   void operator()(std::size_t begin, std::size_t end) {
     
     for (int j = begin; j < end; j++) {
       double phi_j = phi;
       double mu0_j = mu0[j];
       for (int i = 0; i < Y_mat.nrow(); i++) {
         if(R_mat(i,j) == 0){
           double s_i = s[i], y_ij = Y_mat(i,j);
           double xb = B_mat(i,0)* beta1[j] + B_mat(i,1)* beta2[j] + alphae[i];
           double xb_new = B_mat(i,0)* beta1_new[j]+B_mat(i,1)* beta2_new[j] + alphae[i];
           double mu_star = exp(mu0_j + xb_new);
           double mu_old  = exp(mu0_j + xb);
           log_mh_vec[j] = log_mh_vec[j] +  phi_j * (log(s_i* mu_old + phi_j) - log(s_i * mu_star + phi_j)) +
             y_ij * (log(s_i * mu_old + phi_j) + log(mu_star) - log(s_i * mu_star + phi_j) - log(mu_old));
           log_mh = log_mh + phi_j * (log(s_i* mu_old + phi_j) - log(s_i * mu_star + phi_j)) +
             y_ij * (log(s_i * mu_old + phi_j) + log(mu_star) - log(s_i * mu_star + phi_j) - log(mu_old));
         }
       }
     }
   }
   
   void join(const MhBetaseqWorker& other) {
     log_mh += other.log_mh;
   }
 };
 
 // [[Rcpp::export]]
 List mh_betaseq(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat, NumericMatrix W_mat,
                     NumericVector mu0, NumericVector alpha, NumericVector s, double phi, NumericVector beta1, NumericVector beta2, 
                     NumericVector beta1_new, NumericVector beta2_new) {
   
   int p = Y_mat.ncol();
   NumericVector log_mh_vec(p);
   NumericVector alphae =  mvmult(W_mat, alpha);
   
   MhBetaseqWorker worker(R_mat, Y_mat, B_mat, mu0, alphae, s, phi, beta1, beta2, beta1_new, beta2_new, log_mh_vec);
   
   parallelReduce(0, p, worker);
   
   List result;
   result["log_mh"] = worker.log_mh;
   result["log_mh_vec"] = worker.log_mh_vec;
   
   return result;
 }
 
 
 struct DeltaUBeta : public Worker {
   const RMatrix<double> R_mat;
   const RMatrix<double> Y_mat;
   const RMatrix<double> B_mat;
   const RVector<double> mu0;
   const RVector<double> alphae;
   const RVector<double> s;
   const double phi;
   const RVector<double> beta1;
   const RVector<double> beta2;
   RVector<double> d1;
   RVector<double> d2;
   
   DeltaUBeta(const NumericMatrix R_mat, const NumericMatrix Y_mat, const NumericMatrix B_mat,
              const NumericVector mu0, const NumericVector alphae, const NumericVector s, const double phi,
              const NumericVector beta1,  const NumericVector beta2,
              const NumericVector d1, const NumericVector d2)
     : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat), mu0(mu0), alphae(alphae), s(s), phi(phi),
       beta1(beta1), beta2(beta2), d1(d1), d2(d2) {}
   
   void operator()(std::size_t begin, std::size_t end) {
     for (std::size_t j = begin; j < end; j++) {
       double mu0_j = mu0[j];
       double phi_j = phi;
       
       for (int i = 0; i < Y_mat.nrow(); i++) {
         if (R_mat(i, j) == 0) {
           double s_i = s[i];
           double y_ij = Y_mat(i, j);
           double mu = exp(mu0_j + B_mat(i, 0) * beta1[j] + B_mat(i, 1) * beta2[j] + alphae[i]);
           double temp = phi_j * (y_ij - s_i * mu) / (s_i * mu + phi_j);
           d1[j] += temp * B_mat(i, 0);
           d2[j] += temp * B_mat(i, 1);
         }
       }
     }
   }
 };


// [[Rcpp::export]]
List deltaUBeta(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat,
                NumericVector mu0, NumericVector alphae, NumericVector s, double phi,
                NumericVector beta1, NumericVector beta2) {
  int p = Y_mat.ncol();
  
  NumericVector d1(p);
  NumericVector d2(p);
  
  DeltaUBeta deltaUBetaCalc(R_mat, Y_mat, B_mat, mu0, alphae, s, phi, beta1, beta2, d1, d2);
  
  parallelFor(0, p, deltaUBetaCalc);
  
  List result;
  result["d1"] = d1;
  result["d2"] = d2;
  return result;
}


 
 
 // [[Rcpp::export]]
 List update_beta(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat, NumericMatrix W_mat,
                       NumericVector beta1,  NumericVector beta2, NumericVector mu0, NumericVector alpha, 
                       NumericVector s, double phi, double sigma_beta1, double sigma_beta2, double tau_beta, int accept_beta){
   
   int p = beta1.length();
   
   NumericVector rnorm1_new = rnorm(p, 0, sqrt(sigma_beta1));
   NumericVector rnorm2_new = rnorm(p, 0, sqrt(sigma_beta2));

   NumericVector alphae =  mvmult(W_mat, alpha);
   
   List Langevin = deltaUBeta(R_mat, Y_mat, B_mat, mu0, alphae, s, phi, beta1, beta2);
   NumericVector L1 = Langevin["d1"];
   NumericVector L2 = Langevin["d2"];
   NumericVector beta1_new =  pow((1-tau_beta*tau_beta),0.5) * beta1  + (1 - pow((1-tau_beta*tau_beta),0.5)) * sigma_beta1 * L1 + tau_beta * rnorm1_new;
   NumericVector beta2_new =  pow((1-tau_beta*tau_beta),0.5) * beta2  +  (1 - pow((1-tau_beta*tau_beta),0.5)) * sigma_beta2 * L2 + tau_beta * rnorm2_new;
   
   
   double log_mh = mh_beta(R_mat, Y_mat, B_mat, mu0, alphae, s, phi, beta1, beta2, beta1_new, beta2_new);
   
   if(log_mh > log(unif_rand())){
     beta1 = beta1_new;
     beta2 = beta2_new;
     accept_beta = accept_beta + 1;
   }
   
   
   List result;
   result["beta1"] = beta1;
   result["beta2"] = beta2;
   result["log_mh"] = log_mh;
   result["accept_beta"] = accept_beta;
   result["L1"] = L1;
   result["L2"] = L2;
   result["beta1_new"] = beta1_new;
   result["beta2_new"] = beta2_new;
   
   return result;
   
 }




// [[Rcpp::export]]
double LoglikCpp(NumericMatrix Y_mat, NumericMatrix B_mat,  NumericMatrix W_mat, NumericVector scale_hat, 
                        NumericVector mu0_hat, NumericVector beta_hat, NumericVector alpha_hat, double phi_hat, 
                        NumericMatrix tags, double thresh = 0.8) {
  int n = Y_mat.nrow();
  int p = Y_mat.ncol();
  
  double loglik = 0;
  NumericVector alphae_hat =  mvmult(W_mat, alpha_hat);
  NumericMatrix lik(n,p);
  
  for (int j = 0; j < p; j++) {
    for (int i = 0; i < n; i++) {
      double mu = exp(mu0_hat[j] + B_mat(i, 0) * beta_hat[j] + B_mat(i, 1) * beta_hat[j + p] + alphae_hat[i]) * scale_hat[i];
      
      lik(i,j) = dnbinom_mu(Y_mat(i, j), phi_hat, mu, false);
      if (lik(i,j) == 0) {
        lik(i,j) = 1e-08;
      }
      if (tags(i, j) >= thresh) {
        lik(i,j) = 1;
      }
      loglik += log(lik(i,j));
    }
  }
  
  return loglik;
  }



struct DeltaUGamma : public Worker {
  const RMatrix<double> R_mat;
  const RMatrix<double> Y_mat;
  const RMatrix<double> B_mat;
  const RVector<double> mu0;
  const RVector<double> alphae;
  const RVector<double> s;
  const double phi;
  RVector<double> beta1;
  RVector<double> beta2;
  RVector<double> gamma1;
  RVector<double> gamma2;
  RVector<double> lambda1;
  RVector<double> lambda2;
  double epsilon;
  RVector<double> d1;
  RVector<double> d2;

  
  DeltaUGamma(const NumericMatrix R_mat, const NumericMatrix Y_mat, const NumericMatrix B_mat,
              const NumericVector mu0, const NumericVector alphae, const NumericVector s, double phi,
              NumericVector beta1, NumericVector beta2,
              NumericVector gamma1, NumericVector gamma2,
              NumericVector lambda1, NumericVector lambda2, double epsilon,
              NumericVector d1, NumericVector d2)
    : R_mat(R_mat), Y_mat(Y_mat), B_mat(B_mat), mu0(mu0), alphae(alphae), s(s), phi(phi),
      beta1(beta1), beta2(beta2),
      gamma1(gamma1), gamma2(gamma2),
      lambda1(lambda1), lambda2(lambda2), epsilon(epsilon),
      d1(d1), d2(d2) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      double mu0_j = mu0[j];
      double phi_j = phi;
      
      for (int i = 0; i < Y_mat.nrow(); i++) {
        if (R_mat(i, j) == 0) {
          double s_i = s[i];
          double y_ij = Y_mat(i, j);
          double mu = exp(mu0_j + B_mat(i, 0) * beta1[j] + B_mat(i, 1) * beta2[j] + alphae[i]);
          double temp = phi_j * (y_ij - s_i * mu) / (s_i * mu + phi_j);
          // d1[j] += temp * B_mat(i, 0) * (static_cast<double>(beta1[j] != 0) + gamma1[j] * dThresh(gamma1[j], lambda1, epsilon));
          // d2[j] += temp * B_mat(i, 1) * (static_cast<double>(beta2[j] != 0) + gamma2[j] * dThresh(gamma2[j], lambda2, epsilon));
          d1[j] += temp * B_mat(i, 0) * (Thresh(gamma1[j], lambda1[j], epsilon) + gamma1[j] * dThresh(gamma1[j], lambda1[j], epsilon));
          d2[j] += temp * B_mat(i, 1) * (Thresh(gamma2[j], lambda2[j], epsilon) + gamma2[j] * dThresh(gamma2[j], lambda2[j], epsilon));
        }
      }
    }
  }
};


// [[Rcpp::export]]
List deltaUGamma(
    NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat, NumericMatrix W_mat,
    NumericVector mu0, NumericVector alpha, NumericVector s, double phi,
    NumericVector beta1, NumericVector beta2,
    NumericVector gamma1, NumericVector gamma2,
    NumericVector lambda1, NumericVector lambda2, double epsilon
) {
  
  int p = Y_mat.ncol();
  NumericVector d1(p);
  NumericVector d2(p);
  NumericVector alphae =  mvmult(W_mat, alpha);

  DeltaUGamma worker(R_mat, Y_mat, B_mat, mu0, alphae, s, phi, beta1, beta2, gamma1, gamma2, 
                              lambda1, lambda2, epsilon, d1, d2);
  
  parallelFor(0, p, worker);
  
  List result;
  result["d1"] = d1;
  result["d2"] = d2;
  
  
  return result;
}



// [[Rcpp::export]]
List update_gamma(NumericMatrix R_mat, NumericMatrix Y_mat, NumericMatrix B_mat,  NumericMatrix W_mat,
                      NumericVector beta1, NumericVector beta2, NumericVector gamma1, NumericVector gamma2, 
                      NumericVector gamma1_new, NumericVector gamma2_new, NumericVector mu0, NumericVector alpha, 
                      NumericVector s, double phi, NumericVector lambda1, NumericVector lambda2, int accept_gamma){
  
  NumericVector beta1_new(gamma1_new.size());
  NumericVector beta2_new(gamma2_new.size());
  
  for(int j = 0; j < gamma1_new.size(); j++){
    if (abs(gamma1_new[j]) <= lambda1[j]) {
      beta1_new[j] = 0;
    } else {
      beta1_new[j] = gamma1_new[j];
    }
    if (abs(gamma2_new[j]) <= lambda2[j]) {
      beta2_new[j] = 0;
    } else {
      beta2_new[j] = gamma2_new[j];
    }
  }
  
  List temp = mh_betaseq(R_mat, Y_mat, B_mat, W_mat, mu0, alpha, s, phi, beta1, beta2, beta1_new, beta2_new);
  
  double log_mh = temp["log_mh"];
  NumericVector log_mh_vec = temp["log_mh_vec"];
  
  if(log_mh > log(unif_rand())){
    gamma1 = gamma1_new;
    beta1 = beta1_new;
    gamma2 = gamma2_new;
    beta2 = beta2_new;
    accept_gamma = accept_gamma + 1;
  }
  
  
  List result;
  result["gamma1"] = gamma1;
  result["gamma2"] = gamma2;
  result["lambda1"] = lambda1;
  result["lambda2"] = lambda2;
  result["beta1"] = beta1;
  result["beta2"] = beta2;
  result["beta1_new"] = beta1_new;
  result["beta2_new"] = beta2_new;
  result["log_mh"] = log_mh;
  result["log_mh_vec"] = log_mh_vec;
  result["accept_gamma"] = accept_gamma;
  
  return result;
  
}

