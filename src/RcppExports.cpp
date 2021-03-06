// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mortalityhazard_beard_cpp
NumericVector mortalityhazard_beard_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_beard_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_beard_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_beard_cpp
NumericVector mortalityhazard_to_prob_beard_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_beard_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_beard_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_beard_binomial_grad_cpp
NumericVector mortalityhazard_beard_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_beard_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_beard_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_gompertz_cpp
NumericVector mortalityhazard_gompertz_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_gompertz_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_gompertz_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_gompertz_cpp
NumericVector mortalityhazard_to_prob_gompertz_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_gompertz_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_gompertz_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_gompertz_binomial_grad_cpp
NumericVector mortalityhazard_gompertz_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_gompertz_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_gompertz_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_kannisto_cpp
NumericVector mortalityhazard_kannisto_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_kannisto_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_kannisto_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_kannisto_cpp
NumericVector mortalityhazard_to_prob_kannisto_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_kannisto_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_kannisto_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_kannisto_binomial_grad_cpp
NumericVector mortalityhazard_kannisto_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_kannisto_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_kannisto_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_lb_cpp
NumericVector mortalityhazard_lb_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_lb_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_lb_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_lb_cpp
NumericVector mortalityhazard_to_prob_lb_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_lb_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_lb_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_lb_binomial_grad_cpp
NumericVector mortalityhazard_lb_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_lb_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_lb_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_logistic_cpp
NumericVector mortalityhazard_logistic_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_logistic_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_logistic_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_logistic_cpp
NumericVector mortalityhazard_to_prob_logistic_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_logistic_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_logistic_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_logistic_binomial_grad_cpp
NumericVector mortalityhazard_logistic_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_logistic_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_logistic_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_makeham_cpp
NumericVector mortalityhazard_makeham_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_makeham_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_makeham_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_makeham_cpp
NumericVector mortalityhazard_to_prob_makeham_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_makeham_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_makeham_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_makeham_binomial_grad_cpp
NumericVector mortalityhazard_makeham_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_makeham_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_makeham_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_perks_cpp
NumericVector mortalityhazard_perks_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_perks_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_perks_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_perks_cpp
NumericVector mortalityhazard_to_prob_perks_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_perks_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_perks_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_perks_binomial_grad_cpp
NumericVector mortalityhazard_perks_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_perks_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_perks_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_quadratic_cpp
NumericVector mortalityhazard_quadratic_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_quadratic_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_quadratic_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_quadratic_cpp
NumericVector mortalityhazard_to_prob_quadratic_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_quadratic_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_quadratic_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_quadratic_binomial_grad_cpp
NumericVector mortalityhazard_quadratic_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_quadratic_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_quadratic_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_weibull_cpp
NumericVector mortalityhazard_weibull_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_weibull_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_weibull_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_to_prob_weibull_cpp
NumericVector mortalityhazard_to_prob_weibull_cpp(NumericVector theta, NumericVector z);
RcppExport SEXP mortfit_mortalityhazard_to_prob_weibull_cpp(SEXP thetaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_to_prob_weibull_cpp(theta, z));
    return rcpp_result_gen;
END_RCPP
}
// mortalityhazard_weibull_binomial_grad_cpp
NumericVector mortalityhazard_weibull_binomial_grad_cpp(NumericVector theta, NumericVector ages, NumericVector Dx, NumericVector Nx);
RcppExport SEXP mortfit_mortalityhazard_weibull_binomial_grad_cpp(SEXP thetaSEXP, SEXP agesSEXP, SEXP DxSEXP, SEXP NxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ages(agesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dx(DxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nx(NxSEXP);
    rcpp_result_gen = Rcpp::wrap(mortalityhazard_weibull_binomial_grad_cpp(theta, ages, Dx, Nx));
    return rcpp_result_gen;
END_RCPP
}
// rmnom_cpp
NumericVector rmnom_cpp(SEXP n, SEXP p);
RcppExport SEXP mortfit_rmnom_cpp(SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type n(nSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(rmnom_cpp(n, p));
    return rcpp_result_gen;
END_RCPP
}
