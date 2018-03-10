optim.fits <- list()

#####################################
## go through each mortality model
## (eventually)
for(mort.model in binomial.models) {

  context(paste("testing optimFitMethod for",
                mort.model@name))

  ## these are the true values of the parameters we'll
  ## test
  #sim.thetas <- list(c(mort.model@theta.default))
  sim.thetas <- list(c(mort.model@theta.default),
                     plyr::laply(mort.model@theta.range, function(x) x[1]),
                     plyr::laply(mort.model@theta.range, function(x) x[2]))

  theta.hats <- list()
  sim.data <- list()
  folded.data <- list()

  for(i in seq_along(sim.thetas)) {

      ## say we have 100 people at each age from 80 to 99
      these.ages <- 80:99

      ## NB: had to increase the population sizes to be sure that simulated
      ##     data had maxlik param vals sufficiently close to 'true' ones for
      ##     logistic model
      sim.data <- try(mort.model@simulate.fn(theta=sim.thetas[[i]],
                                             Nx=rep(1e7,length(these.ages)),
                                             age=these.ages,
                                             means=TRUE))

      test_that(paste0(mort.model@name, " param set ", i, " - creating simulated data"),
                expect_is(sim.data, "mortalityData"))

      # debugonce(mort.model@loglik.fn)
      true.ll <- try(mort.model@loglik.fn(theta=sim.thetas[[i]],
                                          Dx=sim.data@data$Dx,
                                          Nx=sim.data@data$Nx,
                                          ages=sim.data@data$age))

      test_that(paste0(mort.model@name, " param set ", i, " - true ll is numeric"),
                {
                  expect_is(true.ll, "numeric")
                  expect_false(is.nan(true.ll))
                })


      ## uncomment to actually plot the hazard and the drawn data
      ## plot(sim.data, mort.model@hazard)

      set.seed(1000)

      this.fit <- try(mort.fit(mort.model,
                               sim.data,
                               optimMultipleFit,
                               verbose=FALSE))

      ## uncomment to actually plot the fit
      ## plot(this.fit)

      test_that("fit object has correct type", expect_is(this.fit, "mortalityFit"))

      ## check to see if this has the right class, ie didn't return an error
      ## before trying to grab theta.hat; this prevents all of the tests from
      ## stopping if one model/param combination fails
      if(is(this.fit, "mortalityFit")) {
        theta.hat <- setNames(this.fit@theta.hat, NULL) # strip off names, which confuse testthat
        est.ll <- this.fit@log.likelihood
      } else {
        theta.hat <- rep(NA, mort.model@num.param)
        est.ll <- NA
      }

      ## for now, we're not testing the actual parameter estimates, since what we care
      ## about is the model fit (ie, the likelihood)
      #test_that(paste0(paste(mort.model@name, " param set ", i, sep=""), " - recovered correct parameter estimates"),
      #          expect_equal(theta.hat, sim.thetas[[i]], tolerance=1e-2))

      test_that(paste0(paste(mort.model@name, " param set ", i, sep=""), " - maximized log-likelihood is the same as the one for true params"),
                expect_equal(true.ll, est.ll, tolerance=1e-4))


    }

}



