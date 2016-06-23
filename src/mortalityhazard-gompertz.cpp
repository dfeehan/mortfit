/************************************************
 * mortalityhazard-gompertz.cpp
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


using namespace Rcpp;

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double gompertz_partial_alpha_partI(double alpha, double beta, double x) {

   double gompertz_partial_alpha_partI_result;
   gompertz_partial_alpha_partI_result = -(-exp(beta*x)/beta + exp(beta*x + beta)/beta)*exp(-alpha*(-exp(beta*x)/beta + exp(beta*x + beta)/beta))/(-1 + exp(-alpha*(-exp(beta*x)/beta + exp(beta*x + beta)/beta)));
   return gompertz_partial_alpha_partI_result;

}

double gompertz_partial_alpha_partII(double beta, double x) {

   double gompertz_partial_alpha_partII_result;
   gompertz_partial_alpha_partII_result = exp(beta*x)/beta - exp(beta*x + beta)/beta;
   return gompertz_partial_alpha_partII_result;

}

double gompertz_partial_beta_partI(double alpha, double beta, double x) {

   double gompertz_partial_beta_partI_result;
   gompertz_partial_beta_partI_result = -alpha*(-x*exp(beta*x)/beta + (x + 1)*exp(beta*x + beta)/beta + exp(beta*x)/pow(beta, 2) - exp(beta*x + beta)/pow(beta, 2))*exp(-alpha*(-exp(beta*x)/beta + exp(beta*x + beta)/beta))/(-1 + exp(-alpha*(-exp(beta*x)/beta + exp(beta*x + beta)/beta)));
   return gompertz_partial_beta_partI_result;

}

double gompertz_partial_beta_partII(double alpha, double beta, double x) {

   double gompertz_partial_beta_partII_result;
   gompertz_partial_beta_partII_result = -alpha*(-x*exp(beta*x)/beta + (x + 1)*exp(beta*x + beta)/beta + exp(beta*x)/pow(beta, 2) - exp(beta*x + beta)/pow(beta, 2));
   return gompertz_partial_beta_partII_result;

}

/*
 * END CODE PASTED IN FROM SAGE
 */

// [[Rcpp::export]]
NumericVector mortalityhazard_gompertz_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];

  for(int i=0; i < len; i++) {
    res[i] = alpha*exp(beta*z[i]);
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_gompertz_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];
  double temp1 = 0.0, temp2=0.0;

  /*
   * we want
   *     1 - exp(-(alpha/beta) * [exp(beta*(z+1)) - exp(beta*z)])
   * but we'll perform this as
   *     1 - [exp(-(alpha/beta)*exp(beta*(z+1))) * exp((alpha/beta)*exp(beta*z))]
   * to try to minimize possible float cancellation problems with the difference of the exps
   */

  for(int i=0; i < len; i++) {
    temp1 = -1.0 * (alpha/beta)*exp(beta*(z[i]+1));
    temp2 = (alpha/beta)*(exp(beta*z[i]));
    res[i] = 1 - (exp(temp1) * exp(temp2));
  }

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_gompertz_binomial_grad_cpp(NumericVector theta,
                                                         NumericVector ages,
                                                         NumericVector Dx,
                                                         NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = exp(theta[0]);
  double beta = theta[1];

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
      partI += (double)Dx[ageidx] * gompertz_partial_alpha_partI(alpha, beta, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * gompertz_partial_alpha_partII(beta, cur_age);

  }

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
      partI += (double)Dx[ageidx] * gompertz_partial_beta_partI(alpha, beta, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * gompertz_partial_beta_partII(alpha, beta, cur_age);

  }

  res[BETA_IDX] = partI + partII;

  return(res);

}

