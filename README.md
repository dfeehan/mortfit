**mortfit**
-------

Fit mortality hazard models, with special attention to humans at old ages.

Package currently under development. Stay tuned.

NOTE
----

* 20160526 - this depends on the version of RcppFaddeeva that is on github and
  *not* the CRAN version. (The github version exports the Rcpp headers properly.)
  So it's necessary to

    install_github("baptiste/RcppFaddeeva@5fb5c32")

  in order to build this version


ideas for tests
---------------

* be sure that the predicted values from all of the pre-defined functional forms always make sense (probs between 0 and 1, numbers of deaths never more than number of people, etc)

* be sure that optim always converges on test datasets

* test for haz.to.prob: one way of doing this is to take fairly easy cases, with known closed-form solutions, and to test the haz.to.prob output against them

To add a hazard function
------------------------

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
	* similarl functions to calculate gradient wrt other params

TODO
----

	* `mortalitymodels.R` 
		- add createBinomialModel() call for new hazard
		- add the new model to the list binomial.models
	* `mortalityhazard-[NAME].R` - main purpose of this file is to create a
	  new `mortalityhazard` object. Most of the work here involves creating
		- the hazard function itself (probably in C++)
		* the gradient function (probably in C++)
		* a function to turn hazards over an age range into conditional probabilities
		  (probably in C++)
		- a function to calculate starting values for optimization
	* `mortalityhazard-[NAME].cpp` will have
		- the hazard function itself
		- the gradient function (may involve a few functions)
		- the hazard to probability function

TODO (legacy)
----

- go through and remove subset(); this should only be used interactively. see https://github.com/hadley/devtools/wiki/Evaluation for useful info



Components of a mortality model
-------------------------------

* a probability model; in all of this work, we use the Binomial model
  (see `mortalitymodels.R`)



