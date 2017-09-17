rm(list = ls())
suppressMessages(library("EpiModelHIV"))

#Laptop;
load("C:/Users/Johanna Chapin/Documents/EpiModelHIV/est/fit.m.rda")
load("C:/Users/Johanna Chapin/Documents/EpiModelHIV/est/fit.p.rda")
load("C:/Users/Johanna Chapin/Documents/EpiModelHIV/est/sj/fit.i.rda")
load("C:/Users/Johanna Chapin/Documents/EpiModelHIV/est/st.rda")

#School;
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.m.rda")
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.p.rda")
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/sj/fit.i.rda")
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/st.rda")


est <- list(fit.m, fit.p, fit.i)

#these are default est/st files to simulate from in epimodel but not my est/st files for age*PT
#data(est)
#data(st)

param <- param_msm(nwstats = st,
                   ai.scale.YY=2.5,
                   ai.scale.OY=2.5,
                   ai.scale.OO=0.5,

                   cond.rr.YY = 0.658,
                   cond.rr.OY = 0.723,
                   cond.rr.OO = 3,


                   #ai.scale = 1.8, #1.4 keeps age prev at 10%
                   #ai.scale.OO=0.829,
                   #cond.rr.YY = 0.1,
                   #cond.rr.OY = 0.763,
                   #cond.rr.OO = 1.261,
                   prep.coverage = 0)
init <- init_msm(nwstats = st,
                 prev.Y = 0.2192,
                 prev.O = 0.3452)
control <- control_msm(simno = 0.253,
                       nsteps = 2600, #1500
                       nsims = 1, #1, 5 is alot - hyac;
                       ncores = 1,
                       save.nwstats = TRUE,
                       verbose.int = 1)
sim28<- netsim(est, param, init, control)
# to see code for netsim, just type netsim into console
# debug(stergm_prep)


#12
par(mfrow = c(1, 1))
plot(sim28)
par(mfrow = c(2, 1))
plot(sim28, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim28, y = "i.prev.O", main = "HIV Prev - Older")




#1
#ai.scale=1, nsteps=2600
plot(sim)
par(mfrow = c(2, 1))
plot(sim, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim, y = "i.prev.O", main = "HIV Prev - Older")

df <- as.data.frame(sim)
df$num[2501:2600]
df$num.Y[2501:2600]
df$num.O[2501:2600]
df$incid[2501:2600]
df$incid.infd.Y[2501:2600]
df$incid.infd.O[2501:2600]
#can't get nwstats

#2
#ai.scale=1.4, nsteps=2600
par(mfrow = c(1, 1))
plot(sim2)
par(mfrow = c(2,1))
plot(sim2, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim2, y = "i.prev.O", main = "HIV Prev - Older")

df2 <- as.data.frame(sim2)
df2$num[2501:2600]
df2$num.Y[2501:2600]
df2$num.O[2501:2600]
df2$incid[2501:2600]
df2$incid.infd.Y[2501:2600]
df2$incid.infd.O[2501:2600]
#can't get nwstats

#3
#ai.scale=1.8, nsteps=2600
par(mfrow = c(1, 1))
plot(sim3)
par(mfrow = c(2, 1))
plot(sim3, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim3, y = "i.prev.O", main = "HIV Prev - Older")

df3 <- as.data.frame(sim3)
df3$num[2501:2600]
df3$num.Y[2501:2600]
df3$num.O[2501:2600]
df3$incid[2501:2600]
df3$incid.infd.Y[2501:2600]
df3$incid.infd.O[2501:2600]
df3$incidrate.Y[2501:2600]
df3$incidrate.O[2501:2600]
#can't get nwstats

#3
#ai.scale=1.8, nsteps=2600
par(mfrow = c(1, 1))
plot(sim3)
par(mfrow = c(2, 1))
plot(sim3, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim3, y = "i.prev.O", main = "HIV Prev - Older")

df3 <- as.data.frame(sim3)
df3$num[2501:2600]
df3$num.Y[2501:2600]
df3$num.O[2501:2600]
df3$incid[2501:2600]
df3$incid.infd.Y[2501:2600]
df3$incid.infd.O[2501:2600]
df3$incidrate.Y[2501:2600]
df3$incidrate.O[2501:2600]
#can't get nwstats

#4
par(mfrow = c(1, 1))
plot(sim4)
par(mfrow = c(2, 1))
plot(sim4, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim4, y = "i.prev.O", main = "HIV Prev - Older")

df4 <- as.data.frame(sim4)
df4$num[2501:2600]
df4$num.Y[2501:2600]
df4$num.O[2501:2600]
df4$incid[2501:2600]
df4$incid.infd.Y[2501:2600]
df4$incid.infd.O[2501:2600]
mean(df4$incidrate.Y[2501:2600])
mean(df4$incidrate.O[2501:2600])
#can't get nwstats

#5
par(mfrow = c(1, 1))
plot(sim5)
par(mfrow = c(2, 1))
plot(sim5, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim5, y = "i.prev.O", main = "HIV Prev - Older")

df5 <- as.data.frame(sim5)
df5$num[2501:2600]
df5$num.Y[2501:2600]
df5$num.O[2501:2600]
df5$incid[2501:2600]
df5$incid.infd.Y[2501:2600]
df5$incid.infd.O[2501:2600]

#6
par(mfrow = c(1, 1))
plot(sim6)
par(mfrow = c(2, 1))
plot(sim6, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim6, y = "i.prev.O", main = "HIV Prev - Older")

df6 <- as.data.frame(sim6)
df6$num[2501:2600]
df6$num.Y[2501:2600]
df6$num.O[2501:2600]
df6$incid[2501:2600]
df6$incid.infd.Y[2501:2600]
df6$incid.infd.O[2501:2600]

#7
par(mfrow = c(1, 1))
plot(sim7)
par(mfrow = c(2, 1))
plot(sim7, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim7, y = "i.prev.O", main = "HIV Prev - Older")

df7 <- as.data.frame(sim7)
df7$num[2501:2600]
df7$num.Y[2501:2600]
df7$num.O[2501:2600]
df7$incid[2501:2600]
df7$incid.infd.Y[2501:2600]
df7$incid.infd.O[2501:2600]

#8
par(mfrow = c(1, 1))
plot(sim8)
par(mfrow = c(2, 1))
plot(sim8, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim8, y = "i.prev.O", main = "HIV Prev - Older")


#9
par(mfrow = c(1, 1))
plot(sim9)
par(mfrow = c(2, 1))
plot(sim9, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim9, y = "i.prev.O", main = "HIV Prev - Older")

#10
par(mfrow = c(1, 1))
plot(sim10)
par(mfrow = c(2, 1))
plot(sim10, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim10, y = "i.prev.O", main = "HIV Prev - Older")


nws<-sim10$stats$nwstats
plot(sim10, type="formation", stats=nws)


#7b - fixed NAs
par(mfrow = c(1, 1))
plot(sim11)
par(mfrow = c(2, 1))
plot(sim11, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim11, y = "i.prev.O", main = "HIV Prev - Older")

df7b <- as.data.frame(sim11)
df7b$num[2501:2600]
df7b$num.Y[2501:2600]
df7b$num.O[2501:2600]
df7b$incid[2501:2600]
df7b$incid.infd.Y[2501:2600]
df7b$incid.infd.O[2501:2600]


# plot(sim)
# df <- as.data.frame(sim)
# df$incid
# df$trans.YY
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
dat <- tx_msm(dat, at)s
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

