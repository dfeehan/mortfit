/************************************************
 * mortalityhazard-perks.cpp
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

double perks_partial_alpha_partI(double alpha, double beta, double delta, double gamma, double x) {

   double perks_partial_alpha_partI_result;
   perks_partial_alpha_partI_result = (-exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x) + 1) + exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x + beta) + 1))/(-beta*delta*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta) + beta*delta*exp(alpha*log(delta*exp(beta*x + beta) + 1)/(beta*delta) + gamma + gamma*log(delta*exp(beta*x) + 1)/beta));
   return perks_partial_alpha_partI_result;

}

double perks_partial_alpha_partII(double beta, double delta, double x) {

   double perks_partial_alpha_partII_result;
   perks_partial_alpha_partII_result = -(-log(delta*exp(beta*x) + 1) + log(delta*exp(beta*x + beta) + 1))/(beta*delta);
   return perks_partial_alpha_partII_result;

}

double perks_partial_beta_partI(double alpha, double beta, double delta, double gamma, double x) {

   double perks_partial_beta_partI_result;
   perks_partial_beta_partI_result = -(-((-alpha*beta*pow(delta, 2)*exp(beta) + beta*pow(delta, 3)*gamma*exp(beta))*exp(2*beta*x) + (-alpha*beta*delta*exp(beta) + beta*pow(delta, 2)*gamma*exp(beta) + x*(pow(delta, 2)*gamma*(beta*exp(beta) - beta) - delta*(alpha*beta*exp(beta) - alpha*beta)))*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta) - (-alpha + delta*gamma + (-alpha*pow(delta, 2)*exp(beta) + pow(delta, 3)*gamma*exp(beta))*exp(2*beta*x) + (pow(delta, 2)*gamma*(exp(beta) + 1) - delta*(alpha*exp(beta) + alpha))*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x) + 1) + (-alpha + delta*gamma + (-alpha*pow(delta, 2)*exp(beta) + pow(delta, 3)*gamma*exp(beta))*exp(2*beta*x) + (pow(delta, 2)*gamma*(exp(beta) + 1) - delta*(alpha*exp(beta) + alpha))*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x + beta) + 1))/((pow(beta, 2)*pow(delta, 3)*exp(2*beta*x + beta) + pow(beta, 2)*delta + pow(delta, 2)*(pow(beta, 2)*exp(beta) + pow(beta, 2))*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta) - (pow(beta, 2)*pow(delta, 3)*exp(2*beta*x + beta + gamma) + pow(beta, 2)*delta*exp(gamma) + pow(delta, 2)*(pow(beta, 2)*exp(beta) + pow(beta, 2))*exp(beta*x + gamma))*exp(alpha*log(delta*exp(beta*x + beta) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x) + 1)/beta));
   return perks_partial_beta_partI_result;

}

double perks_partial_beta_partII(double alpha, double beta, double delta, double gamma, double x) {

   double perks_partial_beta_partII_result;
   perks_partial_beta_partII_result = ((-alpha*beta*pow(delta, 2)*exp(beta) + beta*pow(delta, 3)*gamma*exp(beta))*exp(2*beta*x) + (-alpha*beta*delta*exp(beta) + beta*pow(delta, 2)*gamma*exp(beta) + x*(pow(delta, 2)*gamma*(beta*exp(beta) - beta) - delta*(alpha*beta*exp(beta) - alpha*beta)))*exp(beta*x) + (-alpha + delta*gamma + (-alpha*pow(delta, 2)*exp(beta) + pow(delta, 3)*gamma*exp(beta))*exp(2*beta*x) + (pow(delta, 2)*gamma*(exp(beta) + 1) - delta*(alpha*exp(beta) + alpha))*exp(beta*x))*log(delta*exp(beta*x) + 1) - (-alpha + delta*gamma + (-alpha*pow(delta, 2)*exp(beta) + pow(delta, 3)*gamma*exp(beta))*exp(2*beta*x) + (pow(delta, 2)*gamma*(exp(beta) + 1) - delta*(alpha*exp(beta) + alpha))*exp(beta*x))*log(delta*exp(beta*x + beta) + 1))/(pow(beta, 2)*pow(delta, 3)*exp(2*beta*x + beta) + pow(beta, 2)*delta + pow(delta, 2)*(pow(beta, 2)*exp(beta) + pow(beta, 2))*exp(beta*x));
   return perks_partial_beta_partII_result;

}

double perks_partial_gamma_partI(double alpha, double beta, double delta, double gamma, double x) {

   double perks_partial_gamma_partI_result;
   perks_partial_gamma_partI_result = (beta*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta) + exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x) + 1) - exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x + beta) + 1))/(-beta*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta) + beta*exp(alpha*log(delta*exp(beta*x + beta) + 1)/(beta*delta) + gamma + gamma*log(delta*exp(beta*x) + 1)/beta));
   return perks_partial_gamma_partI_result;

}

double perks_partial_gamma_partII(double beta, double delta, double x) {

   double perks_partial_gamma_partII_result;
   perks_partial_gamma_partII_result = -(beta + log(delta*exp(beta*x) + 1) - log(delta*exp(beta*x + beta) + 1))/beta;
   return perks_partial_gamma_partII_result;

}

double perks_partial_delta_partI(double alpha, double beta, double delta, double gamma, double x) {

   double perks_partial_delta_partI_result;
   perks_partial_delta_partI_result = ((pow(delta, 2)*gamma*(exp(beta) - 1) - delta*(alpha*exp(beta) - alpha))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + beta*x + gamma*log(delta*exp(beta*x + beta) + 1)/beta) - (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x) + 1) + (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta)*log(delta*exp(beta*x + beta) + 1))/((beta*pow(delta, 4)*exp(2*beta*x + beta) + beta*pow(delta, 2) + pow(delta, 3)*(beta*exp(beta) + beta)*exp(beta*x))*exp(alpha*log(delta*exp(beta*x) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x + beta) + 1)/beta) - (beta*pow(delta, 4)*exp(2*beta*x + beta + gamma) + beta*pow(delta, 2)*exp(gamma) + pow(delta, 3)*(beta*exp(beta) + beta)*exp(beta*x + gamma))*exp(alpha*log(delta*exp(beta*x + beta) + 1)/(beta*delta) + gamma*log(delta*exp(beta*x) + 1)/beta));
   return perks_partial_delta_partI_result;

}

double perks_partial_delta_partII(double alpha, double beta, double delta, double gamma, double x) {

   double perks_partial_delta_partII_result;
   perks_partial_delta_partII_result = ((pow(delta, 2)*gamma*(exp(beta) - 1) - delta*(alpha*exp(beta) - alpha))*exp(beta*x) - (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*log(delta*exp(beta*x) + 1) + (alpha*pow(delta, 2)*exp(2*beta*x + beta) + alpha + delta*(alpha*exp(beta) + alpha)*exp(beta*x))*log(delta*exp(beta*x + beta) + 1))/(beta*pow(delta, 4)*exp(2*beta*x + beta) + beta*pow(delta, 2) + pow(delta, 3)*(beta*exp(beta) + beta)*exp(beta*x));
   return perks_partial_delta_partII_result;

}


/*
 * END CODE PASTED IN FROM SAGE
 */

/*
 * TODO - LEFT OFF HERE
 *        NEED TO FILL IN ALL OF THE FUNCTIONS BELOW
 * TODO ALSO -
 *        NEED METHOD FOR PICKING STARTING VALUES
 *        AND NEED ROUGH PARAMETER RANGES
 * TODO - unit tests for hazard to prob
 *
 */

// [[Rcpp::export]]
NumericVector mortalityhazard_perks_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[ALPHA_IDX]);
  double beta = exp(theta[BETA_IDX]);
  double gamma = exp(theta[GAMMA_IDX]);
  double delta = exp(theta[DELTA_IDX]);

  double k = 0.0;
  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    k = exp(beta * z[i]);
    num = gamma + (alpha*k);
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
NumericVector mortalityhazard_to_prob_perks_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[ALPHA_IDX]);
  double beta = exp(theta[BETA_IDX]);
  double gamma = exp(theta[GAMMA_IDX]);
  double delta = exp(theta[DELTA_IDX]);

  double temp = 0.0;


  for(int i=0; i < len; i++) {
    temp = (delta*exp(beta*z[i]) + 1.0) / (delta*exp(beta*(z[i]+1)) + 1.0);

    res[i] = 1.0 - exp(-1.0*gamma) * 
               pow(temp, alpha/(beta*delta)) * 
               pow(temp, gamma/beta);
  }

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_perks_binomial_grad_cpp(NumericVector theta,
                                                         NumericVector ages,
                                                         NumericVector Dx,
                                                         NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = exp(theta[ALPHA_IDX]);
  double beta = exp(theta[BETA_IDX]);
  double gamma = exp(theta[GAMMA_IDX]);
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
      partI += (double)Dx[ageidx] * perks_partial_alpha_partI(alpha, beta, delta, gamma, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * perks_partial_alpha_partII(beta, delta, cur_age);

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
      partI += (double)Dx[ageidx] * perks_partial_beta_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * perks_partial_beta_partII(alpha,
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
      partI += (double)Dx[ageidx] * perks_partial_gamma_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * perks_partial_gamma_partII(beta, delta, cur_age);

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
      partI += (double)Dx[ageidx] * perks_partial_delta_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * perks_partial_delta_partII(alpha,
                                                                               beta,
                                                                               delta,
                                                                               gamma,
                                                                               cur_age);

  }

  //Rcout << "partial wrt delta: " << partI + partII << std::endl;

  res[DELTA_IDX] = partI + partII;

  return(res);

}

