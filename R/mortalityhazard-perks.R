########################
# perks hazard object


## these starting values have been updated based on preliminary analysis
perks.haz   <- new("mortalityHazard",
                   name="Perks",
                   num.param=4L,
                   theta.default=c(-2.95, 
                                   -1.9, 
                                   0, 
                                   -1.7),
                   theta.range=list(c(-3.43, -2.76),
                                    c(-2.5, -1.8),
                                    c(-160, -4.7),
                                    c(-7.4, -1.75)),
                   theta.start.fn=function(data.obj) {

                     ## choose starting values by getting b from
                     ## a logistic regression...
                     dat <- data.obj@data
                     crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                     
                     crude[crude > 1] <- .99999
                     crude[crude == 0] <- .000001

                     ## TODO -- we really should insist that the second
                     ##   parameter start off positive; otherwise, we have
                     ##   initially downward-sloping hazards

                     roughmod <- glm(crude ~ age,
                                     data=dat,
                                     family=quasibinomial(link="logit"),
                                     weights=dat$Nx)
                     fv <- fitted.values(roughmod)

                     mid <- mean(range(fv))
                     age.mid.idx <- which.min(abs(fv-mid))
                     age.mid <- dat$age[age.mid.idx]

                     ## hack to handle edge cases...
                     if(age.mid.idx==1) {
                       age.mid.idx <- 2
                     }
                     if(age.mid.idx==length(dat$age)) {
                       age.mid.idx <- age.mid.idx-1
                     }

                     ## this gives us a rough estimate of the slope
                     ## at the inflection point
                     inflect.slope <- (fv[age.mid.idx-1]+fv[age.mid.idx+1])/
                       (dat$age[age.mid.idx+1]-dat$age[age.mid.idx-1])

                     ## NB: below, in the expressions for beta and
                     ## alpha, we should use max(crude) - gamma,
                     ## but we don't have an estimate of gamma yet, so
                     ## we'll just use max(crude)
                     
                     ## starting value for beta
                     beta.init <- 4*inflect.slope/(max(crude))
                     if (beta.init < 0) {
                       beta.init <- 1e-5
                     } else if (beta.init >= .3) {
                       beta.init <- .3
                     }

                     ## starting value for delta
                     delta.init <- 1 / exp(beta.init*age.mid)

                     ## in some cases, where central death rates
                     ## look much more linear than sigmoid, initial delta parameter
                     ## is waaaay too small. so put a limit on how small it can be.
                     if (delta.init < exp(-3)) {
                       delta.init <- exp(-3)
                     }
                     
                     ## starting value for alpha
                     alpha.init <- (max(crude))*delta.init

                     ## starting value for gamma
                     gamma.init <- min(crude[1:min(5,length(crude))]) - 
                                       (alpha.init / (1 + delta.init))
                     #gamma.init <- min(crude[1:min(5,length(crude))]) - 
                     #                  (alpha.init / (1 + alpha.init))
                     
                     if (gamma.init <= 0) {
                       gamma.init <- 1e-9
                     }

                     return(c(log(alpha.init), log(beta.init),
                              log(gamma.init), log(delta.init)))
                   },
                   optim.default=list(method="BFGS",
                                      control=list(
                                                   parscale=c(-3.9, 
                                                              -1.8,
                                                              -5,
                                                              -2.6),
                                                   reltol=1e-10,
                                                   maxit=50000)),
                                                   ##trace=6)),                      
                    haz.fn=mortalityhazard_perks_cpp,
                    binomial.grad.fn=mortalityhazard_perks_binomial_grad_cpp,
                    haz.to.prob.fn=mortalityhazard_to_prob_perks_cpp)
