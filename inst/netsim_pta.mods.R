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
                   ai.scale.YY=1.437,
                   ai.scale.OY=1.279,
                   ai.scale.OO=0.829,

                   cond.rr.YY = 0.658,
                   cond.rr.OY = 0.723,
                   cond.rr.OO = 1.322,


                   #ai.scale = 1.8, #1.4 keeps age prev at 10%
                   #ai.scale.OO=0.829,
                   #cond.rr.YY = 0.1,
                   #cond.rr.OY = 0.763,
                   #cond.rr.OO = 1.261,
                   prep.coverage = 0)
init <- init_msm(nwstats = st,
                 prev.Y = 0.2192,
                 prev.O = 0.3452)

control <- control_msm(simno = 1,
                        nsteps = 2600, #1500
                        nsims = 1, # 5 is alot - hyac;
                        ncores = 1,
                        save.nwstats = TRUE,
                        verbose.int = 1)
sim_prac<- netsim(est, param, init, control)

sim31<- netsim(est, param, init, control)
# to see code for netsim, just type netsim into console
# debug(stergm_prep)

#simno to define what sim scenario

#31
par(mfrow = c(1, 1))
plot(sim31)
par(mfrow = c(2, 1))
plot(sim31, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim31, y = "i.prev.O", main = "HIV Prev - Older")





## additional simulation options below -- using scenario 31


#### A ####
#checking on the NAs for incidence vectors - fixed, now 0's

control <- control_msm(simno = 1,
                       nsteps = 500, #1500
                       nsims = 2, #1, 5 is alot - hyac;
                       ncores = 1,
                       save.nwstats = TRUE,
                       verbose.int = 1)
sim31a<- netsim(est, param, init, control)
# to see code for netsim, just type netsim into console
# debug(stergm_prep)


#31a
par(mfrow = c(1, 1))
plot(sim31a)
par(mfrow = c(2, 1))
plot(sim31a, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim31a, y = "i.prev.O", main = "HIV Prev - Older")

df <- as.data.frame(sim31a)
df$num
df$num.Y
df$num.O
df$incid
df$incid.infd.Y
df$incid.infd.O
mean(df$incidrate.Y)
mean(df$incidrate.O)

sim31a$stats
#this works to get nwstats with conc and agedisc
#across all timesteps
class(sim31a$stats)
dfstat<-as.data.frame(sim31a$stats)

dfstat$nwstats.sim1.agedisc.2
mean(dfstat$nwstats.sim1.agedisc.2)
mean(dfstat$nwstats.sim1.agedisc.1)
mean(dfstat$nwstats.sim1.agedisc)
#gives a new variable name to all columns based on sim#.var.network
#where sim#.var has no number on end and is network 1 (main)
#and sim#.var.1 is network 2 (casual) and sim#.var.2 is network 3 (one-off)

#how to i get mean of age disc for instantaneous
#across all timesteps of 1 sim?
#across all timesteps of all sims, when multiple?
#matrix class
#class(object)
#convert to df

#in general, here's what you're getting in epi after multiple sims:
#columns = sims (500)
#rows = timesteps (100)
#these each make up a df for each variable


#plotting aspects of nwstats
plot(sim31a, type="formation", stats=sim31a$stats)
#Error: One or more requested stats not contained in netsim object


#### B ####
#checking on the NAs for incidence vectors - fixed, now 0's

control <- control_msm(simno = 0.253,
                        nsteps = 2609, #1500
                        nsims = 3, #1, 5 is alot - hyac;
                        ncores = 3,
                        save.nwstats = TRUE,
                        verbose.int = 1)
sim31b<- netsim(est, param, init, control)
# to see code for netsim, just type netsim into console
# debug(stergm_prep)

#in verbose, how to include sim number in monitoring output when running?
#just print(s)
#what is diff again between simno and simulation number? simno is the
#version of simulation you are running -- so like scenario 31;

#
# complete stable version of EpiModelHIV
# pushed to github
# output of session info after have loaded up all software
# sessionInfo() make sure all other pacjages are compatible
# script used to generate the individual simulations
# network stats file
# model fit file


#31b
par(mfrow = c(1, 1))
plot(sim31b)
par(mfrow = c(2, 1))
plot(sim31b, y = "i.prev.Y", main = "HIV Prev - Young")
plot(sim31b, y = "i.prev.O", main = "HIV Prev - Older")

#do we use plot to choose which 10-year interval to take?
# take 10-year when running interventions and need time for intervention to
# have impact on incidence, but if just looking at incidence then can just
# look at the last 2 years and best to look at the last timeframe -- definitely
# stable and is the most 'recent' if were real time
# if incidence is low as it may be when you have many categories like age,
# then you need to up the number of timesteps so you dont have so many 0 and
# so you dont get so much variation/instability/wonky CIs
# so then you take the mean of the var at each timestep across all sims
# then take mean of the mean val for the var across all timesteps


df <- as.data.frame(sim31b)
df$num[2506:2609]
df$num.Y[2506:2609]
df$num.O[2506:2609]
df$incid[2506:2609]
df$incid.infd.Y[2506:2609]
df$incid.infd.O[2501:2609]
mean(df$incidrate.Y[2506:2609])
mean(df$incidrate.O[2506:2609])

sim31b$stats
#this works to get nwstats with conc and agedisc
#across all timesteps
dfstat_b<-as.data.frame(sim31b$stats)

dfstat_b$nwstats.sim5.agedisc.2
mean(dfstat_b$nwstats.sim5.agedisc.2)
mean(dfstat_b$nwstats.sim1.agedisc.1)
mean(dfstat_b$nwstats.sim1.agedisc)
#at all stages, simulated partnerships networks were checked
#to ensure that they indeed retained the expected cross-sectional
#structure and relational durations throughout the simulations

# to get across sims, i think need to create a df that has timesteps x sims
# for each var to then get mean across sims and then across timesteps

df$trans.YY
mean(df$trans.YY)
mean(df$trans.OY)
mean(df$trans.OO)

#need to look back at croi code, but how to get sum stats for analysis...
#croi code below:

sim <- truncate_sim(sim31b, at=2506)
       #left truncates and leaves vector from timestep at:nsteps

sim$epi #df for each variable with timesteps x sims

#remove nan's from data
for (j in 1:length(sim$epi)) {
  for (i in 1:ncol(sim$epi[[j]])) {
    to.fix <- which(is.nan(sim$epi[[j]][,i]))
    if (length(to.fix) > 0) {
      sim$epi[[j]][,i][to.fix] <- 0
    }
  }
}

for (j in 1:length(sim$epi)) {
  for (i in 1:ncol(sim$epi[[j]])) {
    to.fix <- which(is.na(sim$epi[[j]][,i]))
    if (length(to.fix) > 0) {
      sim$epi[[j]][,i][to.fix] <- 0
    }
  }
}


#recall rows=timesteps and col=sims, each value in frame is trans.YY for sim# at timestep XXXX.
#summing over columns, sums over the timesteps of each simulation...?
sim$epi$trans.YY
colMeans(sim$epi$trans.YY)
colMeans(sim$epi$trans.OY)
colMeans(sim$epi$trans.OO)
rowMeans(sim$epi$trans.YY)
mean(rowMeans(sim$epi$trans.YY)) #if you want one overall answer - but we want to know variation over all the simulations

f.trans.YY<-unname(colMeans(sim$epi$trans.YY))
round(quantile(f.trans.YY, probs=c(0.025,0.5,0.975)),3)

median(f.trans.YY)
mean(f.trans.YY)


f.trans.Y<-unname(colMeans(sim$epi$trans.Y))
round(quantile(f.trans.Y, probs=c(0.025,0.5,0.975)),3)
#21% with 5 sims -- CDC estimates 27-29% for 13-24;
#recall that estimating hiv prevalence of 7% pre-18 so if these prev are new
#infections since young guys (likely) then that's about right 21+7=28


#okay, one issue next is thinking about denominator and may need to go back to trans/prev modules
#if want to express certain #s from a different denom
#for example, we have % of all transmissions by YYmain, YYcas, etc.
#what i want to say...among all transmissions that occurred to Y people, % that came from main partnerships
#denominator isn't all that are infected but rather all young that were infected. Maybe I did this before
#but check since this would require slightly different denom.

#need to go through all output vars and tables to make sure that everything there before handing off to sam to run
f.trans.main<-unname(colMeans(sim$epi$trans.main))
round(quantile(f.trans.main, probs=c(0.025,0.5,0.975)),3)
f.trans.casl<-unname(colMeans(sim$epi$trans.casl))
round(quantile(f.trans.casl, probs=c(0.025,0.5,0.975)),3)
f.trans.inst<-unname(colMeans(sim$epi$trans.inst))
round(quantile(f.trans.inst, probs=c(0.025,0.5,0.975)),3)

f.trans.recpt.sus<-unname(colMeans(sim$epi$trans.recpt.sus))
round(quantile(f.trans.recpt.sus, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus<-unname(colMeans(sim$epi$trans.inst.sus))
round(quantile(f.trans.inst.sus, probs=c(0.025,0.5,0.975)),3)

f.trans.stage.act<-unname(colMeans(sim$epi$trans.stage.act))
round(quantile(f.trans.stage.act, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.chr<-unname(colMeans(sim$epi$trans.stage.chr))
round(quantile(f.trans.stage.chr, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.aids<-unname(colMeans(sim$epi$trans.stage.aids))
round(quantile(f.trans.stage.aids, probs=c(0.025,0.5,0.975)),3)

f.trans.undx<-unname(colMeans(sim$epi$trans.undx))
round(quantile(f.trans.undx, probs=c(0.025,0.5,0.975)),3)
f.trans.notinitiated<-unname(colMeans(sim$epi$trans.notinitiated))
round(quantile(f.trans.notinitiated, probs=c(0.025,0.5,0.975)),3)
f.trans.notretained<-unname(colMeans(sim$epi$trans.notretained))
round(quantile(f.trans.notretained, probs=c(0.025,0.5,0.975)),3)
f.trans.partsup<-unname(colMeans(sim$epi$trans.partsup))
round(quantile(f.trans.partsup, probs=c(0.025,0.5,0.975)),3)
f.trans.fullsup<-unname(colMeans(sim$epi$trans.fullsup))
round(quantile(f.trans.fullsup, probs=c(0.025,0.5,0.975)),3)

f.trans.Y<-unname(colMeans(sim$epi$trans.Y))  #recall age of infector not infected;
round(quantile(f.trans.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.O<-unname(colMeans(sim$epi$trans.O))
round(quantile(f.trans.O, probs=c(0.025,0.5,0.975)),3)

f.trans.YY<-unname(colMeans(sim$epi$trans.YY))
round(quantile(f.trans.YY, probs=c(0.025,0.5,0.975)),3)
f.trans.OY<-unname(colMeans(sim$epi$trans.OY))
round(quantile(f.trans.OY, probs=c(0.025,0.5,0.975)),3)
f.trans.OO<-unname(colMeans(sim$epi$trans.OO))
round(quantile(f.trans.OO, probs=c(0.025,0.5,0.975)),3)
f.trans.OYd<-unname(colMeans(sim$epi$trans.OYd)) #recall age of infector to infected
round(quantile(f.trans.OYd, probs=c(0.025,0.5,0.975)),3)
f.trans.YOd<-unname(colMeans(sim$epi$trans.YOd)) #and denom still all infections
round(quantile(f.trans.YOd, probs=c(0.025,0.5,0.975)),3)

sim$epi$incid.infd.Y
unname(colMeans(sim$epi$incid.infd.Y))
sim$epi$incid.OYd
unname(colMeans(sim$epi$incid.OYd))
f.trans.OYd.Y<- unname(colMeans(sim$epi$incid.OYd))/unname(colMeans(sim$epi$incid.infd.Y))
                #means over timesteps: sum(inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
round(quantile(f.trans.OYd.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.YY.Y<- unname(colMeans(sim$epi$incid.YY))/unname(colMeans(sim$epi$incid.infd.Y))
round(quantile(f.trans.YY.Y, probs=c(0.025,0.5,0.975)),3)
##yes, checks out 70.6% of Y infections come from O and 29.4% come from Y
##added these variables to trans and prev modules for future;

#age of infector + PT
f.trans.Ymain<-unname(colMeans(sim$epi$trans.Ymain))
round(quantile(f.trans.Ymain, probs=c(0.025,0.5,0.975)),3)
f.trans.Omain<-unname(colMeans(sim$epi$trans.Omain))
round(quantile(f.trans.Omain, probs=c(0.025,0.5,0.975)),3)
f.trans.Ycasl<-unname(colMeans(sim$epi$trans.Ycasl))
round(quantile(f.trans.Ycasl, probs=c(0.025,0.5,0.975)),3)
f.trans.Ocasl<-unname(colMeans(sim$epi$trans.Ocasl))
round(quantile(f.trans.Ocasl, probs=c(0.025,0.5,0.975)),3)
f.trans.Yinst<-unname(colMeans(sim$epi$trans.Yinst))
round(quantile(f.trans.Yinst, probs=c(0.025,0.5,0.975)),3)
f.trans.Oinst<-unname(colMeans(sim$epi$trans.Yinst))
round(quantile(f.trans.Oinst, probs=c(0.025,0.5,0.975)),3)


#agecombo + PT
f.trans.YYmain<-unname(colMeans(sim$epi$trans.YYmain))
round(quantile(f.trans.YYmain, probs=c(0.025,0.5,0.975)),3)
f.trans.OYmain<-unname(colMeans(sim$epi$trans.OYmain))
round(quantile(f.trans.OYmain, probs=c(0.025,0.5,0.975)),3)
f.trans.OOmain<-unname(colMeans(sim$epi$trans.OOmain))
round(quantile(f.trans.OOmain, probs=c(0.025,0.5,0.975)),3)

f.trans.YYcasl<-unname(colMeans(sim$epi$trans.YYcasl))
round(quantile(f.trans.YYcasl, probs=c(0.025,0.5,0.975)),3)
f.trans.OYcasl<-unname(colMeans(sim$epi$trans.OYcasl))
round(quantile(f.trans.OYcasl, probs=c(0.025,0.5,0.975)),3)
f.trans.OOcasl<-unname(colMeans(sim$epi$trans.OOcasl))
round(quantile(f.trans.OOcasl, probs=c(0.025,0.5,0.975)),3)

f.trans.YYinst<-unname(colMeans(sim$epi$trans.YYinst))
round(quantile(f.trans.YYinst, probs=c(0.025,0.5,0.975)),3)
f.trans.OYinst<-unname(colMeans(sim$epi$trans.OYinst))
round(quantile(f.trans.OYinst, probs=c(0.025,0.5,0.975)),3)
f.trans.OOinst<-unname(colMeans(sim$epi$trans.OOinst))
round(quantile(f.trans.OOinst, probs=c(0.025,0.5,0.975)),3)

#directional;
f.trans.OYdmain<-unname(colMeans(sim$epi$trans.OYdmain))
round(quantile(f.trans.OYdmain, probs=c(0.025,0.5,0.975)),3)
f.trans.YOdmain<-unname(colMeans(sim$epi$trans.YOdmain))
round(quantile(f.trans.YOdmain, probs=c(0.025,0.5,0.975)),3)

f.trans.OYdcasl<-unname(colMeans(sim$epi$trans.OYdcasl))
round(quantile(f.trans.OYdcasl, probs=c(0.025,0.5,0.975)),3)
f.trans.YOdcasl<-unname(colMeans(sim$epi$trans.YOdcasl))
round(quantile(f.trans.YOdcasl, probs=c(0.025,0.5,0.975)),3)

f.trans.OYdinst<-unname(colMeans(sim$epi$trans.OYdinst))
round(quantile(f.trans.OYdinst, probs=c(0.025,0.5,0.975)),3)
f.trans.YOdinst<-unname(colMeans(sim$epi$trans.YOdinst))
round(quantile(f.trans.YOdinst, probs=c(0.025,0.5,0.975)),3)


##after adding in output params to measure PAFs among infections to young msm;
control <- control_msm(simno = 1,
                       nsteps = 2609, #1500
                       nsims = 10, #1, 5 is alot - hyac;
                       ncores = 4,
                       save.nwstats = TRUE,
                       verbose.int = 1)
sim31c<- netsim(est, param, init, control)

control <- control_msm(simno = 1,
                       nsteps = 2609, #1500
                       nsims = 20, #1, 5 is alot - hyac;
                       ncores = 4,
                       save.nwstats = TRUE,
                       verbose.int = 1)
sim31d<- netsim(est, param, init, control)

sim <- truncate_sim(sim31d, at=2506) #last 104 weeks (2 years)
#left truncates and leaves vector from timestep at:nsteps

sim31d$epi$trans.YY.Y
sim<-sim31d

#remove nan's from data
for (j in 1:length(sim$epi)) {
  for (i in 1:ncol(sim$epi[[j]])) {
    to.fix <- which(is.nan(sim$epi[[j]][,i]))
    if (length(to.fix) > 0) {
      sim$epi[[j]][,i][to.fix] <- 0
    }
  }
}

for (j in 1:length(sim$epi)) {
  for (i in 1:ncol(sim$epi[[j]])) {
    to.fix <- which(is.na(sim$epi[[j]][,i]))
    if (length(to.fix) > 0) {
      sim$epi[[j]][,i][to.fix] <- 0
    }
  }
}

sim$epi$trans.YY.Y

f.trans.Y<-unname(colMeans(sim$epi$trans.Y))
round(quantile(f.trans.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.O<-unname(colMeans(sim$epi$trans.O))
round(quantile(f.trans.O, probs=c(0.025,0.5,0.975)),3)
f.trans.main<-unname(colMeans(sim$epi$trans.main))
round(quantile(f.trans.main, probs=c(0.025,0.5,0.975)),3)
f.trans.casl<-unname(colMeans(sim$epi$trans.casl))
round(quantile(f.trans.casl, probs=c(0.025,0.5,0.975)),3)
f.trans.inst<-unname(colMeans(sim$epi$trans.inst))
round(quantile(f.trans.inst, probs=c(0.025,0.5,0.975)),3)


#checking new params for among young infections;
f.trans.OYd.Y<-unname(colMeans(sim$epi$trans.OYd.Y))
round(quantile(f.trans.OYd.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.YY.Y<-unname(colMeans(sim$epi$trans.YY.Y))
round(quantile(f.trans.YY.Y, probs=c(0.025,0.5,0.975)),3)
colMeans(sim$epi$trans.OYd.Y)+colMeans(sim$epi$trans.YY.Y)
colMeans(sim$epi$trans.O)+colMeans(sim$epi$trans.Y)
colMeans(sim$epi$trans.OYd)+colMeans(sim$epi$trans.YOd)+colMeans(sim$epi$trans.YY)+colMeans(sim$epi$trans.OO)
colMeans(sim$epi$trans.OY)+colMeans(sim$epi$trans.YY)+colMeans(sim$epi$trans.OO)

f.trans.main.Y<-unname(colMeans(sim$epi$trans.main.Y))
round(quantile(f.trans.main.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.casl.Y<-unname(colMeans(sim$epi$trans.casl.Y))
round(quantile(f.trans.casl.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.Y<-unname(colMeans(sim$epi$trans.inst.Y))
round(quantile(f.trans.inst.Y, probs=c(0.025,0.5,0.975)),3)

f.trans.OYdmain.Y<-unname(colMeans(sim$epi$trans.OYdmain.Y))
round(quantile(f.trans.OYdmain.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.YYmain.Y<-unname(colMeans(sim$epi$trans.YYmain.Y))
round(quantile(f.trans.YYmain.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.OYdcasl.Y<-unname(colMeans(sim$epi$trans.OYdcasl.Y))
round(quantile(f.trans.OYdcasl.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.YYcasl.Y<-unname(colMeans(sim$epi$trans.YYcasl.Y))
round(quantile(f.trans.YYcasl.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.OYdinst.Y<-unname(colMeans(sim$epi$trans.OYdinst.Y))
round(quantile(f.trans.OYdinst.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.YYinst.Y<-unname(colMeans(sim$epi$trans.YYinst.Y))
round(quantile(f.trans.YYinst.Y, probs=c(0.025,0.5,0.975)),3)
##all of these currently null--> that's bc hadn't set these particular ones
##using [at], just these ones tho--fixed now.


f.trans.stage.act.Y.Y<-unname(colMeans(sim$epi$trans.stage.act.Y.Y))
round(quantile(f.trans.stage.act.Y.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.chr.Y.Y<-unname(colMeans(sim$epi$trans.stage.chr.Y.Y))
round(quantile(f.trans.stage.chr.Y.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.aids.Y.Y<-unname(colMeans(sim$epi$trans.stage.aids.Y.Y))
round(quantile(f.trans.stage.aids.Y.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.act.O.Y<-unname(colMeans(sim$epi$trans.stage.act.O.Y))
round(quantile(f.trans.stage.act.O.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.chr.O.Y<-unname(colMeans(sim$epi$trans.stage.chr.O.Y))
round(quantile(f.trans.stage.chr.O.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.stage.aids.O.Y<-unname(colMeans(sim$epi$trans.stage.aids.O.Y))
round(quantile(f.trans.stage.aids.O.Y, probs=c(0.025,0.5,0.975)),3)

f.trans.recpt.sus<-unname(colMeans(sim$epi$trans.recpt.sus))
round(quantile(f.trans.recpt.sus, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus<-unname(colMeans(sim$epi$trans.inst.sus))
round(quantile(f.trans.inst.sus, probs=c(0.025,0.5,0.975)),3)


f.trans.recpt.sus.Y<-unname(colMeans(sim$epi$trans.recpt.sus.Y))
round(quantile(f.trans.recpt.sus.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.recpt.sus.O<-unname(colMeans(sim$epi$trans.recpt.sus.O))
round(quantile(f.trans.recpt.sus.O, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus.Y<-unname(colMeans(sim$epi$trans.inst.sus.Y))
round(quantile(f.trans.inst.sus.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus.O<-unname(colMeans(sim$epi$trans.inst.sus.O))
round(quantile(f.trans.inst.sus.O, probs=c(0.025,0.5,0.975)),3)


f.trans.recpt.sus.amongY<-unname(colMeans(sim$epi$trans.recpt.sus.amongY))
round(quantile(f.trans.recpt.sus.amongY, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus.amongY<-unname(colMeans(sim$epi$trans.inst.sus.amongY))
round(quantile(f.trans.inst.sus.amongY, probs=c(0.025,0.5,0.975)),3)

dat$epi$trans.recpt.sus.Y.Y
dat$epi$trans.inst.sus.Y.Y
dat$epi$trans.recpt.sus.O.Y
dat$epi$trans.inst.sus.O.Y

f.trans.recpt.sus.Y.Y<-unname(colMeans(sim$epi$trans.recpt.sus.Y.Y))
round(quantile(f.trans.recpt.sus.Y.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.recpt.sus.O.Y<-unname(colMeans(sim$epi$trans.recpt.sus.O.Y))
round(quantile(f.trans.recpt.sus.O.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus.O.Y<-unname(colMeans(sim$epi$trans.inst.sus.O.Y))
round(quantile(f.trans.inst.sus.O.Y, probs=c(0.025,0.5,0.975)),3)
f.trans.inst.sus.Y.Y<-unname(colMeans(sim$epi$trans.inst.sus.Y.Y))
round(quantile(f.trans.inst.sus.Y.Y, probs=c(0.025,0.5,0.975)),3)



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
at<-at+1
for (at in 2:100) {
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
}


