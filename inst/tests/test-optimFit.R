optim.fits <- list()

#####################################
## go through each mortality model
## (eventually)
for(mort.model in all.models) {

  context(paste("testing optimFitMethod for",
                mort.model@name))
  
  ## these are the true values of the parameters we'll
  ## test
  sim.thetas <- list(c(mort.model@theta.default),
                     c(.75*mort.model@theta.default),
                     c(1.25*mort.model@theta.default))
  theta.hats <- list()
  sim.data <- list()
  folded.data <- list()
  
  ##cat("\nStarting tests for ", mort.model@name, "\n")
  
  for(i in seq_along(sim.thetas)) {

      ## say we have 100 people at each age from 80 to 99
      these.ages <- 80:99
      
      ## NB: had to increase the population sizes to be sure that simulated
      ##     data had maxlik param vals sufficiently close to 'true' ones for
      ##     logistic model
      sim.data <- mort.model@simulate.fn(theta=sim.thetas[[i]],
                                         Nx=rep(100000,length(these.ages)),
                                         age=these.ages,
                                         means=TRUE)

      true.ll <- mort.model@loglik.fn(theta=sim.thetas[[i]],
                                      Dx=sim.data@data$Dx,
                                      Nx=sim.data@data$Nx,
                                      ages=(these.ages-min(these.ages)))
      
      ## TODO -- eventually test on folded data...
      ##folded.data <- partition.into.folds(5, sim.data)

      ## uncomment to actually plot the hazard and the drawn data
      ##sim.data.plots[[mort.model@name]] <- plot(these.sim.data, weibull.haz)

      this.fit <- try(mort.fit(mort.model,
                               ##folded.data,
                               sim.data,
                               optimFit,
                               verbose=FALSE))

      expect_that(this.fit, is_a("mortalityFit"))

      ## check to see if this has the right class, ie didn't return an error
      ## before trying to grab theta.hat; this prevents all of the tests from
      ## stopping if one model/param combination fails
      if(is(this.fit, "mortalityFit")) {
        theta.hat <- this.fit@theta.hat
        est.ll <- this.fit@log.likelihood
      } else {
        theta.hat <- rep(NA, mort.model@num.params)
        est.ll <- NA
      }

      relerr <- abs((theta.hat-sim.thetas[[i]])/sim.thetas[[i]])
      thistest <- all(relerr < 0.05)
      ##expect_true(thistest, paste(mort.model@name, " param set ", i, sep=""))
      expect_that(thistest,
                  is_true(),
                  paste(mort.model@name, " param set ", i, sep=""),                  
                  "fitted estimates of theta within 0.05 relative error of true ones")

      expect_that(true.ll,
                  equals(est.ll),
                  paste(mort.model@name, " param set ", i, sep=""),
                  "maximized log-likelihood is the same as the one for true params")

    }

}



