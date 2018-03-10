#########################################################
## weibull hazard object

## (these starting values have been updated based on preliminary analysis)
weibull.haz <- new("mortalityHazard",
                   name="Weibull",
                   num.param=2L,
                   theta.default=c(-3, .1),
                   #theta.range=list(c(-5, -1e-5),
                   theta.range=list(c(-5, .1),
                                    c(-2, .5)),
                   theta.start.fn=function(data.obj) {

                     dat <- data.obj@data
                     crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                     crude[crude==0] <- 1e-5
                     prelim.t <- coef(lm(log(crude) ~ log(age),
                                         weights=dat$Nx,
                                         data=dat))
                     #return(c(log(prelim.t[1]), log(prelim.t[2])-1))
                     #return(c(log(prelim.t[1]), log(prelim.t[2]+1)))
                     return(c(prelim.t[1], log(prelim.t[2]+1)))

                   },
                   ##NB: was Nelder-Mead
                   optim.default=list(method="BFGS",
                                      control=list(reltol=1e-10,
                                                   parscale=c(-3.2,
                                                              -0.4),
                                                   maxit=10000)),
                   haz.fn=mortalityhazard_weibull_cpp,
                   #binomial.grad.fn=NULL,
                   binomial.grad.fn=mortalityhazard_weibull_binomial_grad_cpp,
                   haz.to.prob.fn=mortalityhazard_to_prob_weibull_cpp)
