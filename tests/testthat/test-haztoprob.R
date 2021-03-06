#########################################################
## test the conversion of hazards to probabilities
##

## tolerance to use for unit tests that produce floating point values
## (this is more stringent than the default, which is about 1.5e-8
unit.test.tol <- 1e-10

#############################
## gompertz hazard

context("testing haz.to.prob - gompertz hazard")

## these are for the parameter values for a fit to switz-f-cohort-1894
## theta = {'alpha':-3.4282, 'beta':0.0939}
gompertz.sage.val <- c(0.0334453292547037,0.0366769172159332,0.0402142036030128,0.0440847674241909,0.0483184005696209,0.0529472275419508,0.0580058184791177,0.0635312915407369,0.0695633998360972,0.0761445970356908,0.0833200746181912,0.0911377623527231,0.0996482820992584,0.108904843332960,0.118963066977598,0.129880722200984,0.141717358835341,0.154533816119608,0.168391586637264,0.183352012806526)

test.theta <- c(-3.4282, 0.0939)

test.gompertz <- gomp.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.gompertz,equals(gompertz.sage.val, tol=unit.test.tol))

#############################
## makeham hazard

context("testing haz.to.prob - makeham hazard")

## these are for the parameter values for a fit to italy-f-cohort-1883
## theta={alpha: 4.307938e-2,
##        beta:  8.749655e-2,
##        gamma: 2.93565e-6},

makeham.sage.val <- c(0.0440246537105089,0.0479519149056272,0.0522198673149450,0.0568562407901143,0.0618906752409653,0.0673547841167738,0.0732822048668431,0.0797086322791085,0.0866718298540695,0.0942116135425198,0.102369801260095,0.111190120594920,0.120718066058182,0.131000696111261,0.142086359068175,0.154024335863695,0.166864386658608,0.180656187409227,0.195448641969770,0.211289055166128)

test.theta <- c(4.307938e-2, 8.749655e-2, 2.93565e-6)
#
test.makeham <- mak.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.makeham,equals(makeham.sage.val, tol=unit.test.tol))


#############################
## logistic hazard

context("testing haz.to.prob - logistic hazard")

## these are for the parameter values for a fit to denma-m-cohort-1841
#theta = {alpha : exp(-2.7117), beta : exp(-1.9032),
#         gamma : exp(-7.8955), delta : exp(-1.5982)}
logistic.sage.val <- c(0.0574450766491218,0.0645257736205574,0.0722024747661948,0.0804575483962248,0.0892568129487079,0.0985485729595442,0.108263659317502,0.118316649492393,0.128608326616152,0.139029287447311,0.149464452748055,0.159798100224883,0.169918959073370,0.179724895138786,0.189126779213444,0.198051253300509,0.206442263582168,0.214261382983907,0.221487073709059,0.228113124338542)

test.theta <- c(-2.7117, -1.9032, -7.8955, -1.5982)

test.logistic <- logistic.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.logistic,equals(logistic.sage.val, tol=unit.test.tol))

#############################
## perks hazard

context("testing haz.to.prob - perks hazard")

## these are for the parameter values for a fit to italy-f-cohort-1883
#alpha : exp(-3.112581),
#beta : exp(-2.152711),
#gamma : exp(-25.493997),
#delta : exp(-2.261383)}

perks.sage.val <- c(0.0415853222361324,0.0460355371130853,0.0508833980416314,0.0561474057527608,0.0618432594317888,0.0679830151217361,0.0745742021426637,0.0816189326082915,0.0891130488530065,0.0970453620692201,0.105397041310721,0.114141213761286,0.123242833450751,0.132658865494268,0.142338816223513,0.152225617016513,0.162256843015464,0.172366219996191,0.182485346774984,0.192545540143129)

test.theta <- c(-3.112581, -2.152711, -25.493997, -2.261383)

test.perks <- perks.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.perks,equals(perks.sage.val, tol=unit.test.tol))

#############################
## beard hazard

context("testing haz.to.prob - beard hazard")

## these are for the parameter values for a fit to italy-f-cohort-1883
##  theta = {alpha : exp(-3.112054),
##           beta : exp(-2.151455),
##           delta : exp(-2.255184)},

beard.sage.val <- c(0.0415842902714916,0.0460371462887652,0.0508877074108572,0.0561543450553039,0.0618525896694496,0.0679942851970466,0.0745867029433015,0.0816316505523181,0.0891246215574681,0.0970540393926370,0.105400655454541,0.114137162292573,0.123228078929974,0.132629954785788,0.142291921508127,0.152156599057160,0.162161335465477,0.172239731669484,0.182323377049485,0.192343701281297)

test.theta <- c(-3.112054, -2.151455, -2.255184)

test.beard <- beard.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.beard,equals(beard.sage.val, tol=unit.test.tol))

#############################
## kannisto hazard

context("testing haz.to.prob - kannisto hazard")

## these are for the parameter values for a fit to denma-m-cohort-1841
## theta = {alpha : exp(-2.7543), beta : exp(-2.26240330336121)},
kannisto.sage.val <- c(0.0609287928160195,0.0669370907411333,0.0734648659485122,0.0805422353113762,0.0881980790913786,0.0964593309007302,0.105350174950482,0.114891159291111,0.125098240534676,0.135981783349115,0.147545546491612,0.159785695752682,0.172689892155351,0.186236510141782,0.200394044202759,0.215120762362367,0.230364660161930,0.246063758655769,0.262146774321965,0.278534168261983)

test.theta <- c(-2.7543, -2.26240330336121)

test.kannisto <- kannisto.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.kannisto,equals(kannisto.sage.val, tol=unit.test.tol))


#############################
## weibull hazard

context("testing haz.to.prob - weibull hazard")

## these are for the parameter values for a fit to denma-m-cohort-1845
test.theta <- log(c(0.0487, 0.4585))

weibull.sage.val <- c(0.100769536063004,0.0389575434064293,0.0293791843778882,0.0244786993700468,0.0213737490955179,0.0191833381054328,0.0175331244951558,0.0162331651809542,0.0151754513234406,0.0142934394392114,0.0135436084426492,0.0128961432921207,0.0123298464772963,0.0118291822516003,0.0113824719088509,0.0109807454027909,0.0106169847384775,0.0102856104360656,0.00998212389437392,0.00970285266402993)

test.weibull <- weibull.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.weibull,equals(weibull.sage.val, tol=unit.test.tol))

#############################
## log-quadratic hazard

context("testing haz.to.prob - log-quadratic hazard")

## these are for the parameter values for a fit to denma-m-cohort-1841
## theta = {'alpha':-2.8905, 'beta':0.1262, 'gamma':-0.0023}
lq.sage.val <- c(0.0574430728656160,0.0646226401604976,0.0723434382034114,
              0.0805901232194092,0.0893378715517511,0.0985520398947749,
              0.108188073700139,0.118191683497049,0.128499297052722,
              0.139038782324000,0.149730422948784,0.160488115613563,
              0.171220748002960,0.181833708048967,0.192230470480138,
              0.202314205550720,0.211989357294424,0.221163144352396,
              0.229746944738201,0.237657535982441)

test.theta <- c(-2.8905, 0.1262, -0.0023)

test.log_quadratic <- quad.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.log_quadratic,equals(lq.sage.val))

## these are for the parameter values for a fit to denma-m-cohort-1841
## theta={alpha:-3.4649, beta:0.1039, gamma:-5e-4}
lq.sage.val <- c(0.0324160203916671,0.0358652529680389,0.0396350632396515,
  0.0437492469642166,0.0482325824635806,0.0531107588024771,
  0.0584102830182595,0.0641583642832294,0.0703827729262864,
  0.0771116723419858,0.0843734219838408,0.0921963498928657,
  0.100608493551380,0.109637308298288,0.119309343093989,
  0.129649884094233,0.140682567283930,0.152428962336345,
  0.164908130889175,0.178136163571156)

test.theta <- c(-3.4649, 0.1039, -5e-4)

test.log_quadratic <- quad.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.log_quadratic,equals(lq.sage.val))

#############################
## lynch-brown hazard

context("testing haz.to.prob - lynch-brown hazard")

## these are for the parameter values for a fit to denma-m-cohort-1875
##theta={alpha:0.528823791482327, beta:exp(-0.4379),
##       gamma:exp(-3.7814), delta:41.5003},
lb.sage.val <- c(0.0428656268798101,0.0504376025487915,0.0581213849351238,0.0659167193500115,0.0738231353800686,0.0818399301760535,0.0899661513584399,0.0982005796520187,0.106541711383618,0.114987740999886,0.123536543786777,0.132185658997708,0.140932273622873,0.149773207058438,0.158704896958571,0.167723386576834,0.176824313923802,0.186002903084583,0.195253958051511,0.204571859433122)

test.theta <- c(0.528823791482327, -0.4379, -3.7814, 41.5003)

test.lb <- lb.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.lb,equals(lb.sage.val))


## NB: using the approximation exp(-Mx) does not
##     approximate numerical results especially well,
##     according to this test

## TODO -- run for a few scenarios and be sure that
##         it generally converges, and also that the values
##         are probabilities (ie, 0 <= phat <= 1)

## TODO -- take some simple examples with easy closed-form answers
##         and check them (Gompertz, Weibull, etc)

## TODO -- eventually, use Mathematica to get some external
##         computations for some examples; these can be good
##         validators
