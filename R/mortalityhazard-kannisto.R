########################
# kannisto
########################

# hazard fn
kannisto.haz.fn <- function(theta, z) {

    alpha <- exp(theta[1])
    beta <- theta[2]

    res.num <- alpha*exp(beta*z)
    res.denom <- 1 + alpha*exp(beta*z)
    res <- (res.num/res.denom)

    ## ensure only non-negative results are returned
    ## (hazards can't be negative)
    res[res<0] <- NA

    return(res)

}

kannisto.haz.to.prob <- function(haz.fn, theta, z) {

    alpha <- exp(theta[1])
    beta <- theta[2]

    ## NB: if values of alpha are negative, then this may not work
    ## (in that case, you might end up trying to take the log of
    ##  a negative number)
    res <- (1/beta)*(log1p(alpha*exp(beta*(z+1)))-log1p(alpha*exp(beta*z)))
    res <- 1 - exp(-res)
    
    return(res)
  
}

## kannisto hazard fn implemented in c++
kannisto.haz.fn.cpp <- function(theta, z) {
  res <- .Call("mortalityhazard_kannisto_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## kannisto hazard fn implemented in c++
kannisto.haz.to.prob.cpp <- function(haz.fn, theta, z) {
  ## note that we don't need the first argument, haz.fn
  res <- .Call("mortalityhazard_to_prob_kannisto_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## these starting values have been updated based on preliminary analysis
kannisto.haz   <- new("mortalityHazard",
                      name="Kannisto",
                      num.param=2L,
                      theta.default=c(.05, .09),
                      theta.range=list(c(log(0.03), log(0.07)),
                                       c(log(0.08), log(0.11))),
                      theta.start.fn=function(data.obj) {
                        ## choose starting values by getting b from
                        ## a logistic regression...
                        dat <- data.obj@data
                        crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                        crude[crude > 1] <- 1
                        prelim.theta <- coef(glm(crude ~ age,data=dat,
                                                 family=quasibinomial(link="logit"),
                                                 weights=dat$Nx))

                        res <- c(prelim.theta[1], prelim.theta[2])
                        if(res[2] > 0) {
                            res[2] <- log(res[2])
                        }

                        return(res)
                      },                      
                      ## NB: method was Nelder-Mead
                      optim.default=list(method="BFGS",
                                         control=list(parscale=c(log(.05), 
                                                                 log(.09)),
                                                      reltol=1e-10,
                                                      maxit=10000)),
                      haz.fn=mortalityhazard_kannisto_cpp,
                      #binomial.grad.fn=NULL,
                      binomial.grad.fn=mortalityhazard_kannisto_binomial_grad_cpp,
                      haz.to.prob.fn=mortalityhazard_to_prob_kannisto_cpp)
