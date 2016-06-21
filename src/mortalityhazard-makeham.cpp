/************************************************
 * mortalityhazard-makeham.cpp
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

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double makeham_partial_alpha_partI(double alpha, double beta, double gamma, double x) {

   double makeham_partial_alpha_partI_result;
   makeham_partial_alpha_partI_result = -(-exp(beta*x)/beta + exp(beta*x + beta)/beta)*exp((alpha*exp(beta*x) + beta*gamma*x)/beta - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/beta)/(exp((alpha*exp(beta*x) + beta*gamma*x)/beta - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/beta) - 1);
   return makeham_partial_alpha_partI_result;

}

double makeham_partial_alpha_partII(double beta, double x) {

   double makeham_partial_alpha_partII_result;
   makeham_partial_alpha_partII_result = exp(beta*x)/beta - exp(beta*x + beta)/beta;
   return makeham_partial_alpha_partII_result;

}

double makeham_partial_beta_partI(double alpha, double beta, double gamma, double x) {

   double makeham_partial_beta_partI_result;
   makeham_partial_beta_partI_result = -(-(alpha*x*exp(beta*x) + gamma*x)/beta + (alpha*(x + 1)*exp(beta*x + beta) + gamma*x + gamma)/beta + (alpha*exp(beta*x) + beta*gamma*x)/pow(beta, 2) - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/pow(beta, 2))*exp((alpha*exp(beta*x) + beta*gamma*x)/beta - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/beta)/(exp((alpha*exp(beta*x) + beta*gamma*x)/beta - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/beta) - 1);
   return makeham_partial_beta_partI_result;

}

double makeham_partial_beta_partII(double alpha, double beta, double gamma, double x) {

   double makeham_partial_beta_partII_result;
   makeham_partial_beta_partII_result = (alpha*x*exp(beta*x) + gamma*x)/beta - (alpha*(x + 1)*exp(beta*x + beta) + gamma*x + gamma)/beta - (alpha*exp(beta*x) + beta*gamma*x)/pow(beta, 2) + (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/pow(beta, 2);
   return makeham_partial_beta_partII_result;

}

double makeham_partial_gamma_partI(double alpha, double beta, double gamma, double x) {

   double makeham_partial_gamma_partI_result;
   makeham_partial_gamma_partI_result = (x - (beta*x + beta)/beta)*exp((alpha*exp(beta*x) + beta*gamma*x)/beta - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/beta)/(exp((alpha*exp(beta*x) + beta*gamma*x)/beta - (alpha*exp(beta*x + beta) + beta*gamma*x + beta*gamma)/beta) - 1);
   return makeham_partial_gamma_partI_result;

}

double makeham_partial_gamma_partII(double beta, double x) {

   double makeham_partial_gamma_partII_result;
   makeham_partial_gamma_partII_result = x - (beta*x + beta)/beta;
   return makeham_partial_gamma_partII_result;

}

/*
 * CODE ABOVE PASTED IN FROM SAGE
 */

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mortalityhazard_makeham_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  int parlen = theta.size();

  if (parlen != 3) {

  }

  NumericVector res(len);

  //double alpha = exp(theta[0]);
  //double beta = theta[1];
  //double gamma = exp(theta[2]);

  double alpha = std::abs(theta[0]);
  double beta = theta[1];
  double gamma = std::abs(theta[2]);

  for(int i=0; i < len; i++) {

    res[i] = gamma + alpha*exp(beta*z[i]);

    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_makeham_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  //double alpha = exp(theta[0]);
  //double beta = theta[1];
  //double gamma = exp(theta[2]);

  double alpha = std::abs(theta[0]);
  double beta = theta[1];
  double gamma = std::abs(theta[2]);


  double temp1 = 0.0, temp2=0.0;

  /*
   * we want
   *     1 - exp(-gamma -(alpha/beta) * [exp(beta*(z+1)) - exp(beta*z)])
   * but we'll perform this as
   *     1 - [exp(gamma) * exp(-(alpha/beta)*exp(beta*(z+1))) * exp((alpha/beta)*exp(beta*z))]
   * to try to minimize possible float cancellation problems with the difference of the exps
   */

  for(int i=0; i < len; i++) {
    temp1 = -1.0 * (alpha/beta)*exp(beta*(z[i]+1));
    temp2 = (alpha/beta)*(exp(beta*z[i]));
    res[i] = 1 - (exp(-1.0*gamma) * exp(temp1) * exp(temp2));
  }

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_makeham_binomial_grad_cpp(NumericVector theta,
                                                         NumericVector ages,
                                                         NumericVector Dx,
                                                         NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  //double alpha = exp(theta[0]);
  //double beta = theta[1];
  //double gamma = exp(theta[2]);

  double alpha = std::abs(theta[0]);
  double beta = theta[1];
  double gamma = std::abs(theta[2]);

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
      partI += (double)Dx[ageidx] * makeham_partial_alpha_partI(alpha, beta, gamma, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * makeham_partial_alpha_partII(beta, cur_age);

  }

  Rcout << "partial wrt alpha: " << partI + partII << std::endl;

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
      partI += (double)Dx[ageidx] * makeham_partial_beta_partI(alpha, beta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * makeham_partial_beta_partII(alpha,
                                                                               beta,
                                                                               gamma,
                                                                               cur_age);

  }

  Rcout << "partial wrt beta: " << partI + partII << std::endl;

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
      partI += (double)Dx[ageidx] * makeham_partial_gamma_partI(alpha, beta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * makeham_partial_gamma_partII(beta, cur_age);

  }

  Rcout << "partial wrt gamma: " << partI + partII << std::endl;

  res[GAMMA_IDX] = partI + partII;

  return(res);

}

