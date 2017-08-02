rm(list = ls())
suppressMessages(library("EpiModelHIV"))

load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.m.rda")
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.p.rda")
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/sj/fit.i.rda")
est <- list(fit.m, fit.p, fit.i)


load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/st.rda")

#these are default est/st files to simulate from in epimodel but not my est/st files for age*PT
#data(est)
#data(st)

param <- param_msm(nwstats = st,
                   ai.scale = 1.0, #1.4 keeps age prev at 10%
                   prep.coverage = 0)
init <- init_msm(nwstats = st,
                 prev.Y = 0.2192,
                 prev.O = 0.3452)
control <- control_msm(simno = 0.253,
                       nsteps = 1300, #1500
                       nsims = 1, #1, 5 is alot - hyac;
                       ncores = 1,
                       save.nwstats = TRUE,
                       verbose.int = 1)
sim <- netsim(est, param, init, control)
# to see code for netsim, just type netsim into console
# debug(stergm_prep)

# df <- as.data.frame(sim)
# df$incid
# df$trans.YY
# plot(sim)
# df$nwstats

param
init
control

#set.seed(86)
#set.seed(Sys.time())

at <- 1
dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)
# mf <- dat$p[[1]]$model.form
# mf$terms[[4]]


at <- at + 1

dat <- aging_msm(dat, at)
dat <- deaths_msm(dat, at)
dat <- births_msm(dat, at)
dat <- test_msm(dat, at)
dat <- tx_msm(dat, at)
dat <- prep_msm(dat, at)
dat <- progress_msm(dat, at)
dat <- vl_msm(dat, at)
# dat <- update_aiclass_msm(dat, at)
# dat <- update_roleclass_msm(dat, at)
dat <- simnet_msm(dat, at)
dat <- disclose_msm(dat, at)
dat <- acts_msm(dat, at)
dat <- condoms_msm(dat, at)
dat <- riskhist_msm(dat, at)
dat <- position_msm(dat, at)
dat <- trans_msm(dat, at)
dat <- prevalence_msm(dat, at)

