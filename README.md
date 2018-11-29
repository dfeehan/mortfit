**mortfit**
-------

Fit mortality hazard models, with special attention to humans at old ages.

NOTES
-----

* 20160526 - this depends on the version of RcppFaddeeva that is on github and
  *not* the CRAN version. (The github version exports the Rcpp headers properly.)
  So it's necessary to

    install_github("baptiste/RcppFaddeeva@5fb5c32")

  in order to build this version

(NB: The issue opened by Tim Riffe on github has a solution)


To add a hazard function
------------------------

Example steps used to create the Beard hazard object:

* create `mortalityhazard-beard.R`
	- include code to create new `mortalityHazard` object
	
* calculate
	1. hazard to probability function
	2. gradients

* Using Sage, I produced C++ code to calculate (1) and (2) above.
  I put the results in `mortalityhazard-beard.cpp`. The function names are:
	* `mortalityhazard_beard_cpp` - calculates the mortality hazard, given params
	* `mortalityhazard_to_prob_beard_cpp` - given params and ages, calculates cond prob of death
	* `mortalityhazard_beard_binomial_grad_cpp` - calculates the gradient of the binomial likelihood
	* `beard_partial_alpha_partI` - calculates first part of gradient wrt alpha
	* `beard_partial_alpha_partII` - calculates second part of gradient wrt alpha
	* similar functions to calculate gradient wrt other params



