/************************************************
 * rmnom.h
 * 
 * header file for the rmnom funcion
 * (uses Rcpp). see rmnom.cpp for more details.
 * 
 * dennis, dec 2011
 ************************************************/

#ifndef _rmnom_RCPP_RMNOM_H
#define _rmnom_RCPP_RMNOM_H

#include <Rcpp.h>

RcppExport SEXP rmnom(SEXP n, SEXP p);

#endif


