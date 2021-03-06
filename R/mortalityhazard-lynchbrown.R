########################
# lynch and brown hazard

## parameter transformation
lb.param.trans <- function(theta) {
    return(c(theta[1],
             exp(theta[2]),
             exp(theta[3]),
             theta[4]))
    #return(c(exp(theta[1]),
    #         exp(theta[2]),
    #         exp(theta[3]),
    #         theta[4]))
}

## hazard function
lb.haz.fn <- function(theta, z) {

    theta.trans <- lb.param.trans(theta)
    alpha <- theta.trans[1]
    beta <- theta.trans[2]
    gamma <- theta.trans[3]
    delta <- theta.trans[4]

    res <- alpha + beta * atan(gamma*(z-delta))

    ## ensure only non-negative results are returned
    ## (hazards can't be negative)
    res[res<0] <- NA

    return(res)
}

lb.haz.to.prob <- function(haz.fn, theta, z) {

    theta.trans <- lb.param.trans(theta)
    alpha <- theta.trans[1]
    beta <- theta.trans[2]
    gamma <- theta.trans[3]
    delta <- theta.trans[4]

    res1 <- alpha
    res2 <- (1/(2*gamma))*beta
    res3 <- 2*gamma*(z-delta)*atan(gamma*(delta-z))
    res4 <- 2*gamma*(delta-z-1)*atan(gamma*(delta-z-1))
    res5 <- log1p((gamma^2)*((delta-z)^2))
    res6 <- log1p((gamma^2)*((delta-z-1)^2))

    res <- res1 + res2*(res3 + res4 + res5 - res6)

    res <- 1 - exp(-res)

    return(res)
}

## lb hazard fn implemented in c++
lb.haz.fn.cpp <- function(theta, z) {
  res <- .Call("mortalityhazard_lb_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## lb hazard fn implemented in c++
lb.haz.to.prob.cpp <- function(haz.fn, theta, z) {
  ## note that we don't need the first argument, haz.fn
  res <- .Call("mortalityhazard_to_prob_lb_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

# these values come from Table 1 of Lynch and Brown (2001).
# I took the point estimates of the parameters they present
# for the baseline year, plus and minus 3 standard deviations
# (I subtracted 79 from the parameter delta, because I'm indexing
#  age from 1, not 80). I log them because of the parameters are
# all constrained to be positive
#d> log(alphas)
#[1] -1.177331 -1.140998 -1.105939
#d> log(betas)
#[1] -1.528319 -1.485010 -1.443500
#d> log(gammas)
#[1] -2.294617 -2.231195 -2.171557
#d> log(deltas)
#[1] 2.743186 2.785011 2.825157

## TODO -- consider adding theta.scale?

## these starting values have been updated based on preliminary analysis
lb.haz   <- new("mortalityHazard",
                name="Lynch-Brown",
                num.param=4L,
                #theta.default=c(-1.57, -1.485, -2.23, 16),
                theta.default=c(1.0, -1.485, -2.23, 16),
                ## NB: these are the ranges that random thetas
                ##     are drawn from; they are not typical ranges
                ##     but rather safe starting values
                #theta.range=list(c(-.1, 1.57),
                theta.range=list(c(0.9, 4.8),
                                 c(-1.53, -1.43),
                                 c(-3.0, -2.0),
                                 c(-5, 22)),

                theta.start.fn=function(data.obj) {

                  ## choose starting values by getting b from
                  ## a regression...
                  dat <- data.obj@data
                  crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)

                  crude[crude > 1] <- .9999
                  crude[crude == 0] <- .000001

                  roughmod <- glm(crude ~ age,
                                  data=dat,
                                  weight=dat$Nx,
                                  family=quasibinomial(link="logit"))
                  fv <- fitted.values(roughmod)

                  mid <- mean(range(fv))
                  age.mid.idx <- which.min(abs(fv-mid))
                  ## age at which fitted probs are closest to the mean
                  ## (one way of thinking about middle age wrt death rates)
                  age.mid <- dat$age[age.mid.idx]

                  ## beta is the maximum value that the hazard rates attain,
                  ## scaled by pi
                  beta <- max(crude)/pi

                  ## beta*gamma is slope at inflection point
                  ## compute an approximation to this numerically

                  ## hack to handle edge cases...
                  if(age.mid.idx==1) {
                    age.mid.idx <- 2
                  }
                  if(age.mid.idx==length(dat$age)) {
                    age.mid.idx <- age.mid.idx-1
                  }

                  gamma <- (fv[age.mid.idx-1]+fv[age.mid.idx+1])/
                           (dat$age[age.mid.idx+1]-dat$age[age.mid.idx-1])

                  if (gamma < 0) {
                    gamma <- 1e-5
                  }

                  ## rough starting value for delta is age at which
                  ## exposure-weighted linear probability model on central
                  ## death rates is halfway to max...
                  ## substitute this in from rough code / logistic reg
                  delta <- age.mid

                  ## pick alpha so that we know all of the hazards
                  ## will be positive at their initial values
                  etmp <- c(beta,gamma,delta)
                  rawval <- etmp[1]*atan(etmp[2]*(1-etmp[3]))

                  rawval.atinf <- etmp[1]*atan(etmp[2]*(age.mid-etmp[3]))

                  if (rawval < 0) {
                    alpha <- -rawval + .01
                  } else {
                    alpha <- rawval
                  }


                  return(c(alpha,log(beta),log(gamma),delta))

                },
                optim.default=list(method="BFGS",
                                   #control=list(parscale=c(0.01,
                                   control=list(parscale=c(0.1,
                                                           0.1,
                                                           0.1,
                                                           0.1),
                                                reltol=1e-10,
                                                maxit=50000)),
                haz.fn=mortalityhazard_lb_cpp,
                #binomial.grad.fn=NULL,
                binomial.grad.fn=mortalityhazard_lb_binomial_grad_cpp,
                haz.to.prob.fn=mortalityhazard_to_prob_lb_cpp)

