/************************************************
 * mortalityhazard-beard.cpp
 *
 *
 *
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 *
 * dennis, dec 2011
 ************************************************/

#include <Rcpp.h>

// NB: DELTA_IDX is 2 here, since gamma is fixed at 0 for
//     the Beard hazard
#define ALPHA_IDX 0
#define BETA_IDX 1
#define DELTA_IDX 2

using namespace Rcpp;

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double beard_partial_alpha_partI(double alpha, double beta, double delta, double x) {

   double beard_partial_alpha_partI_result;
   beard_partial_alpha_partI_result = (-exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta))*log(beta*delta*exp(beta*x) + beta) + exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta))*log(beta*delta*exp(beta*x + beta) + beta))/(-beta*delta*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta)) + beta*delta*exp(alpha*log(beta*delta*exp(beta*x + beta) + beta)/(beta*delta)));
   return beard_partial_alpha_partI_result;

}

double beard_partial_alpha_partII(double beta, double delta, double x) {

   double beard_partial_alpha_partII_result;
   beard_partial_alpha_partII_result = -(-log(beta*delta*exp(beta*x) + beta) + log(beta*delta*exp(beta*x + beta) + beta))/(beta*delta);
   return beard_partial_alpha_partII_result;

}

double beard_partial_beta_partI(double alpha, double beta, double delta, double x) {

   double beard_partial_beta_partI_result;
   beard_partial_beta_partI_result = -(-(alpha*beta*pow(delta, 2)*exp(2*beta*x + beta) + (alpha*beta*delta*exp(beta) + delta*x*(alpha*beta*exp(beta) - alpha*beta))*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta)) - (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta))*log(beta*delta*exp(beta*x) + beta) + (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta))*log(beta*delta*exp(beta*x + beta) + beta))/(-(pow(beta, 2)*pow(delta, 3)*exp(2*beta*x + beta) + pow(beta, 2)*delta + pow(delta, 2)*(pow(beta, 2)*exp(beta) + pow(beta, 2))*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta)) + (pow(beta, 2)*pow(delta, 3)*exp(2*beta*x + beta) + pow(beta, 2)*delta + pow(delta, 2)*(pow(beta, 2)*exp(beta) + pow(beta, 2))*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x + beta) + beta)/(beta*delta)));
   return beard_partial_beta_partI_result;

}

double beard_partial_beta_partII(double alpha, double beta, double delta, double x) {

   double beard_partial_beta_partII_result;
   beard_partial_beta_partII_result = -(alpha*beta*pow(delta, 2)*exp(2*beta*x + beta) + (alpha*beta*delta*exp(beta) + delta*x*(alpha*beta*exp(beta) - alpha*beta))*exp(beta*x) + (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*log(beta*delta*exp(beta*x) + beta) - (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*log(beta*delta*exp(beta*x + beta) + beta))/(pow(beta, 2)*pow(delta, 3)*exp(2*beta*x + beta) + pow(beta, 2)*delta + pow(delta, 2)*(pow(beta, 2)*exp(beta) + pow(beta, 2))*exp(beta*x));
   return beard_partial_beta_partII_result;

}

double beard_partial_delta_partI(double alpha, double beta, double delta, double x) {

   double beard_partial_delta_partI_result;
   beard_partial_delta_partI_result = (delta*(alpha*exp(beta) - alpha)*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta) + beta*x) + (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta))*log(beta*delta*exp(beta*x) + beta) - (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta))*log(beta*delta*exp(beta*x + beta) + beta))/(-(beta*pow(delta, 4)*exp(2*beta*x + beta) + beta*pow(delta, 2) + pow(delta, 3)*(beta*exp(beta) + beta)*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x) + beta)/(beta*delta)) + (beta*pow(delta, 4)*exp(2*beta*x + beta) + beta*pow(delta, 2) + pow(delta, 3)*(beta*exp(beta) + beta)*exp(beta*x))*exp(alpha*log(beta*delta*exp(beta*x + beta) + beta)/(beta*delta)));
   return beard_partial_delta_partI_result;

}

double beard_partial_delta_partII(double alpha, double beta, double delta, double x) {

   double beard_partial_delta_partII_result;
   beard_partial_delta_partII_result = -(delta*(alpha*exp(beta) - alpha)*exp(beta*x) + (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*log(beta*delta*exp(beta*x) + beta) - (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*log(beta*delta*exp(beta*x + beta) + beta))/(beta*pow(delta, 4)*exp(2*beta*x + beta) + beta*pow(delta, 2) + pow(delta, 3)*(beta*exp(beta) + beta)*exp(beta*x));
   return beard_partial_delta_partII_result;

}


/*
 * END CODE PASTED IN FROM SAGE
 */

// [[Rcpp::export]]
NumericVector mortalityhazard_beard_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[ALPHA_IDX]);
  double beta = exp(theta[BETA_IDX]);
  double delta = exp(theta[DELTA_IDX]);

  double k = 0.0;
  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    k = exp(beta * z[i]);
    num = alpha * k;
    denom = 1 + (delta*k);

    /*Rprintf("denom is %f\n", denom);*/

    res[i] = (num/denom);
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
NumericVector mortalityhazard_to_prob_beard_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[ALPHA_IDX]);
  double beta = exp(theta[BETA_IDX]);
  double delta = exp(theta[DELTA_IDX]);

  double temp = 0.0;


  for(int i=0; i < len; i++) {
    temp = (delta*exp(beta*z[i]) + 1.0) / (delta*exp(beta*(z[i]+1)) + 1.0);
    res[i] = 1.0 - pow(temp, alpha/(beta*delta));
  }

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_beard_binomial_grad_cpp(NumericVector theta,
                                                      NumericVector ages,
                                                      NumericVector Dx,
                                                      NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = exp(theta[ALPHA_IDX]);
  double beta = exp(theta[BETA_IDX]);
  double delta = exp(theta[DELTA_IDX]);

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
      partI += (double)Dx[ageidx] * beard_partial_alpha_partI(alpha, 
                                                              beta, 
                                                              delta, 
                                                              cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * beard_partial_alpha_partII(beta, 
                                                                               delta, 
                                                                               cur_age);

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
      partI += (double)Dx[ageidx] * beard_partial_beta_partI(alpha, 
                                                             beta, 
                                                             delta, 
                                                             cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * beard_partial_beta_partII(alpha,
                                                                             beta,
                                                                             delta,
                                                                             cur_age);

  }

  //Rcout << "partial wrt beta: " << partI + partII << std::endl;

  res[BETA_IDX] = partI + partII;

  /*
   * delta
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * beard_partial_delta_partI(alpha, 
                                                              beta, 
                                                              delta, 
                                                              cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * beard_partial_delta_partII(alpha,
                                                                               beta,
                                                                               delta,
                                                                               cur_age);

  }

  //Rcout << "partial wrt delta: " << partI + partII << std::endl;

  res[DELTA_IDX] = partI + partII;

  return(res);

}

