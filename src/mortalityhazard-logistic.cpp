/************************************************
 * mortalityhazard-logistic.cpp
 *
 *
 *
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 *
 * dennis, dec 2011
 ************************************************/

#include <Rcpp.h>

#define ALPHA_IDX 0
#define BETA_IDX 1
#define GAMMA_IDX 2
#define DELTA_IDX 3

using namespace Rcpp;

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double logistic_partial_alpha_partI(double alpha, double beta, double delta, double gamma, double x) {

   double logistic_partial_alpha_partI_result;
      logistic_partial_alpha_partI_result = -(-log(beta*delta*exp(beta*x) + beta) + log(beta*delta*exp(beta*x + beta) + beta))*exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta))/(beta*delta*(-1 + exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta))));
         return logistic_partial_alpha_partI_result;

}

double logistic_partial_alpha_partII(double beta, double delta, double x) {

   double logistic_partial_alpha_partII_result;
      logistic_partial_alpha_partII_result = -(-log(beta*delta*exp(beta*x) + beta) + log(beta*delta*exp(beta*x + beta) + beta))/(beta*delta);
         return logistic_partial_alpha_partII_result;

}

double logistic_partial_beta_partI(double alpha, double beta, double delta, double gamma, double x) {

   double logistic_partial_beta_partI_result;
      logistic_partial_beta_partI_result = -((alpha*(beta*delta*(x + 1)*exp(beta*x + beta) + delta*exp(beta*x + beta) + 1)/(beta*delta*exp(beta*x + beta) + beta) - alpha*(beta*delta*x*exp(beta*x) + delta*exp(beta*x) + 1)/(beta*delta*exp(beta*x) + beta) + delta*gamma)/(beta*delta) - (-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(pow(beta, 2)*delta))*exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta))/(-1 + exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta)));
         return logistic_partial_beta_partI_result;

}

double logistic_partial_beta_partII(double alpha, double beta, double delta, double gamma, double x) {

   double logistic_partial_beta_partII_result;
      logistic_partial_beta_partII_result = -(alpha*(beta*delta*(x + 1)*exp(beta*x + beta) + delta*exp(beta*x + beta) + 1)/(beta*delta*exp(beta*x + beta) + beta) - alpha*(beta*delta*x*exp(beta*x) + delta*exp(beta*x) + 1)/(beta*delta*exp(beta*x) + beta) + delta*gamma)/(beta*delta) + (-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(pow(beta, 2)*delta);
         return logistic_partial_beta_partII_result;

}

double logistic_partial_gamma_partI(double alpha, double beta, double delta, double gamma, double x) {

   double logistic_partial_gamma_partI_result;
      logistic_partial_gamma_partI_result = -exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta))/(-1 + exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta)));
         return logistic_partial_gamma_partI_result;

}

double logistic_partial_gamma_partII() {

   double logistic_partial_gamma_partII_result;
      logistic_partial_gamma_partII_result = -1;
         return logistic_partial_gamma_partII_result;

}

double logistic_partial_delta_partI(double alpha, double beta, double delta, double gamma, double x) {

   double logistic_partial_delta_partI_result;
      logistic_partial_delta_partI_result = -((alpha*beta*exp(beta*x + beta)/(beta*delta*exp(beta*x + beta) + beta) - alpha*beta*exp(beta*x)/(beta*delta*exp(beta*x) + beta) + beta*gamma)/(beta*delta) - (-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*pow(delta, 2)))*exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta))/(-1 + exp(-(-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*delta)));
         return logistic_partial_delta_partI_result;

}

double logistic_partial_delta_partII(double alpha, double beta, double delta, double gamma, double x) {

   double logistic_partial_delta_partII_result;

      logistic_partial_delta_partII_result = -(alpha*beta*exp(beta*x + beta)/(beta*delta*exp(beta*x + beta) + beta) - alpha*beta*exp(beta*x)/(beta*delta*exp(beta*x) + beta) + beta*gamma)/(beta*delta) + (-alpha*log(beta*delta*exp(beta*x) + beta) + alpha*log(beta*delta*exp(beta*x + beta) + beta) + beta*delta*gamma)/(beta*pow(delta, 2));

      return logistic_partial_delta_partII_result;

}

/*
 * END CODE PASTED IN FROM SAGE
 */

// [[Rcpp::export]]
NumericVector mortalityhazard_logistic_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = exp(theta[3]);

  double k = 0.0;
  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    k = exp(beta * z[i]);
    num = alpha * k;
    denom = 1 + (delta*k);

    /*Rprintf("denom is %f\n", denom);*/

    res[i] = (num/denom) + gamma;
  }

  bool anyNA = is_true( any( is_na(res) ) );

  if (anyNA) {
    for(int i=0; i<len; i++) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_logistic_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = exp(theta[3]);

  double temp = 0.0;


  for(int i=0; i < len; i++) {
    temp = (delta*exp(beta*z[i]) + 1.0) / (delta*exp(beta*(z[i]+1)) + 1.0);
    res[i] = 1.0 - exp(-1.0*gamma) * pow(temp, alpha/(beta*delta));
  }

  //for(int i=0; i < len; i++) {
  //  k0 = exp(beta * z[i]);
  //  k1 = exp(beta * (z[i]+1));
  //  temp = gamma + (alpha/(beta*delta))*(log(beta*delta*k1 + beta) -
  //                                       log(beta*delta*k0 + beta));
  //  res[i] = 1 - exp(-1.0*temp);
  //}

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_logistic_binomial_grad_cpp(NumericVector theta,
                                                         NumericVector ages,
                                                         NumericVector Dx,
                                                         NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = exp(theta[3]);

  double partI = 0.0, partII =0.0;

  double cur_age = 0.0;

  // initialize result vector to 0
  for(int i=0; i < theta_len; i++) {
    res[i] = 0.0;
  }

  /*
   * alpha
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * logistic_partial_alpha_partI(alpha, beta, delta, gamma, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * logistic_partial_alpha_partII(beta, delta, cur_age);

  }

  //Rcout << "partial wrt alpha: " << partI + partII << std::endl;

  res[ALPHA_IDX] = partI + partII;

  /*
   * beta
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * logistic_partial_beta_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * logistic_partial_beta_partII(alpha,
                                                                               beta,
                                                                               delta,
                                                                               gamma,
                                                                               cur_age);

  }

  //Rcout << "partial wrt beta: " << partI + partII << std::endl;

  res[BETA_IDX] = partI + partII;

  /*
   * gamma
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * logistic_partial_gamma_partI(alpha, beta, delta, gamma, cur_age);
      // partII is always -1, but we keep the function call here for clarity and maintainability
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * logistic_partial_gamma_partII();

  }

  //Rcout << "partial wrt gamma: " << partI + partII << std::endl;

  res[GAMMA_IDX] = partI + partII;

  /*
   * delta
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * logistic_partial_delta_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * logistic_partial_delta_partII(alpha,
                                                                               beta,
                                                                               delta,
                                                                               gamma,
                                                                               cur_age);

  }

  //Rcout << "partial wrt delta: " << partI + partII << std::endl;

  res[DELTA_IDX] = partI + partII;

  return(res);

}

