########################
# quadratic hazard

## hazard fn
quad.haz.fn <- function(theta, z) {
     alpha <- theta[1]
     beta <-  theta[2]
     gamma <- theta[3]

     res <- exp(alpha + beta*z + gamma*z^2)

     ## ensure only non-negative results are returned
     ## (hazards can't be negative)
     res[res<0] <- NA

     return(res)
}



## from mathematica (pade[5,5] expansion about x=2)
erfi.approx <- function(x) {

  res <- (-5.733966903607422e9 + x*(-1.0422029063023047e11 + 
        x*(1.166096291419065e11 + x*
            (-6.238574268109686e10 + (1.6683214106948364e10 - 2.282425831401211e9*x)*x)
           )))/(-1.1624511006316472e11 + 
     x*(1.5419057578654575e11 + x*(-8.212176055901636e10 + 
           x*(2.1900269307177174e10 + 
              x*(-2.9145905967631354e9 + 1.541031251468339e8*x)))))

  return(res)
}

## from mathematica (pade[5,5] expansion about x=0.3) of erf(x)
erf.approx <- function(x) {
  res <- (-1.794456570098843e-10 + x*(1.0390620375805752 + 
        x*(0.18186305861283816 + x*
            (0.001502720433999738 + (0.056049232331913376 + 0.018889241003572095*x)*x))
        ))/(0.9208447490906735 + x*
      (0.16117203632722585 + x*(0.3082790329869982 + 
           x*(0.10340277455907947 + (0.027385632367503975 + 0.018447945290923147*x)*x))
        ))

  return(res)
}

## from mathematica, pade[5,5] expansion about x=2 of erfi(x)
erfi.approx.nr2 <- function(x) {

  erfiof2 <- 18.5648
  expof4 <- exp(4)
  sqrtofpi <- sqrt(pi)
  
  res <- (39421363760*expof4 - 65584280236*sqrtofpi*erfiof2 + 
     x*(-54337699400*expof4 + 86992716740*sqrtofpi*erfiof2 + 
        x*(30059331440*expof4 - 46332241890*sqrtofpi*erfiof2 + 
           x*(-8589300460*expof4 + 12355903820*sqrtofpi*erfiof2 + 
              x*(1296600940*expof4 - 1644381655*sqrtofpi*erfiof2 + 
                 x*(-94203190*expof4 + 86943378*sqrtofpi*erfiof2))))))/
   (-65584280236*sqrtofpi + x*(86992716740*sqrtofpi + 
        x*(-46332241890*sqrtofpi + 
           x*(12355903820*sqrtofpi + x*(-1644381655*sqrtofpi + 86943378*sqrtofpi*x)))))

  return(res)

}

## from mathematica, pade[7,7] expansion about x=1 of erfi(x)
erfi.approx.nr1 <- function(x) {

  erfiof1 <- 1.65043
  expof1 <- exp(1)
  sqrtofpi <- sqrt(pi)

  res <- (173465337818777746628*expof1 - 161189343277473823463*sqrtofpi*erfiof1 + 
     x*(-181236135883415078228*expof1 + 58207050515281020997*sqrtofpi*erfiof1 + 
        x*(-32990023275939643620*expof1 + 70448356034529178041*sqrtofpi*erfiof1 + 
           x*(53037366273114757940*expof1 - 37842691311251513855*sqrtofpi*erfiof1 + 
              x*(-7666360342517050580*expof1 - 5516314131408582425*sqrtofpi*erfiof1 + 
                 x*(-6665334075150129708*expof1 + 7528886304748135803*sqrtofpi*erfiof1 + 
                    x*(2371517339152075268*expof1 - 1823214457143623537*sqrtofpi*erfiof1 + 
                       x*(-316367854022677700*expof1 + 146437170211591039*sqrtofpi*erfiof1))
                    ))))))/
   (-161189343277473823463*sqrtofpi + 
     x*(58207050515281020997*sqrtofpi + 
        x*(70448356034529178041*sqrtofpi + 
           x*(-37842691311251513855*sqrtofpi + 
              x*(-5516314131408582425*sqrtofpi + 
                 x*(7528886304748135803*sqrtofpi + 
                    x*(-1823214457143623537*sqrtofpi + 146437170211591039*sqrtofpi*x)))
              ))))

  return(res)

}

quad.haz.to.prob <- function(haz.fn, theta, z) {
  
  alpha <- theta[1]
  beta <-  theta[2]
  gamma <- as.complex(theta[3])

  k <- exp(alpha - ((beta^2)/(4*gamma)))*sqrt(pi)

  if (! is.finite(k)) {
    return(rep(NA, length(z)))
  }
  
  res <- (k/(2*sqrt(gamma)))*(erfi.approx.nr1((beta+ (2*gamma*(z+1)))/
                                              (2*sqrt(gamma)))-
                              erfi.approx.nr1((beta+(2*gamma*(z)))/
                                              (2*sqrt(gamma))))

  res <- 1 - exp(-res)
  res <- Re(res)

  ##browser()
  
  return(res)

}

## quad hazard fn implemented in c++
quad.haz.fn.cpp <- function(theta, z) {
  res <- .Call("mortalityhazard_quadratic_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## quad hazard fn implemented in c++
quad.haz.to.prob.cpp <- function(haz.fn, theta, z) {
  ## note that we don't need the first argument, haz.fn
  res <- .Call("mortalityhazard_to_prob_quadratic_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## these starting values have been updated based on preliminary analysis
quad.haz <- new("mortalityHazard",
                name="Log-Quadratic",
                num.param=3L,
                theta.default=c(-3.04, .097, -8e-4),
                theta.range=list(c(-3.46, -2.76),
                                 c(0.071, 0.118),
                                 c(-.0018, 0.0003)),

                theta.start.fn=function(data.obj) {
                  ## choose starting values via a
                  ## regression on the log of the crude rates
                  dat <- data.obj@data
                  crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                  crude[crude<=0] <- 1e-9
                  dat$age2 <- dat$age^2
                  prelim.b <- coef(lm(log(crude) ~ age + age2,
                                      weights=dat$Nx,
                                      data=dat))

                  return(prelim.b)
                },                
                ## NB: was Nelder-Mead
                optim.default=list(method="BFGS",
                                   control=list(parscale=c(-1, .01, 1e-4),
                                                reltol=1e-12)),
                haz.fn=mortalityhazard_quadratic_cpp,
                #binomial.grad.fn=NULL,
                binomial.grad.fn=mortalityhazard_quadratic_binomial_grad_cpp,
                haz.to.prob.fn=mortalityhazard_to_prob_quadratic_cpp)



