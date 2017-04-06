## Partner Type By Age ESTIM file ##

suppressMessages(library(EpiModelHIV))
rm(list = ls())

# load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/st.rda")
# load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nw.main.rda")
# load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nw.pers.rda")



###############################################################
## CONCURRENT OVERALL (NOT BY AGE), AGE*PT DISSOLUTION RATES ##
###############################################################

#Ok, assuming tergm lite issue can be resolved, adding in
#differential dissolution rates by pt and age
#this requires recoding the estimation file for ‘st’
#with not only an ~edges term but also a nodemix term for agecat2
#requires coding main and casual models with different rates
#by age and PT


# 1. Main Model -----------------------------------------------------------


# Initialize network
nw.main <- base_nw_msm(st)


# Assign degree -- assign casual degree to main network (because global constraint)
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)
#nw.main
#save(nw.main, file = "C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nw.main.rda")

# Formulas
formation.m <- ~edges + nodemix("agecat2", base = 1) +
                nodefactor("deg.pers") +
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

#think what are the factors that go into the network formation -- the formation of main edges;
#the number of casual partners they have
#age homophily
#age - main degrees differ by age;


# Fit model
fit.m <- netest(nw.main,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.m,
                coef.diss = st$coef.diss.m,
                constraints = ~bd(maxout = 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e10, MCMLE.maxit = 250))

fit.m
summary(fit.m)

#4/5/2017
save(fit.m, file="C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.m.rda")

# Main Diagnostics -------------------------------------------------------------

#new with homophily;
dx <- netdx(fit.m, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
            nwstats.formula = ~edges + nodemix("agecat2", base = 1) +
                              nodefactor("deg.pers") +
                              absdiff("sqrt.age") +
                              offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
dx
plot(dx)




# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw.pers <- nw.main

# Assign degree -- assign main degree to casual network (because global constraint)
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st)
#nw.pers
#save(nw.pers, file = "C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nw.pers.rda")

# Formulas
formation.p <- ~edges + nodemix("agecat2", base = 1) +
  nodefactor("deg.main") +
  concurrent +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, keep = 1:2))


# Fit model --with approximation
fit.p <- netest(nw.pers,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.p,
                coef.diss = st$coef.diss.p,
                constraints = ~bd(maxout = 2),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9, MCMLE.maxit = 250))

fit.p
summary(fit.p)

#4/5/2017
save(fit.p, file="C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.p.rda")

# Casual Diagnostics -------------------------------------------------------------

#new, with homophily;
dx.p <- netdx(fit.p, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
              nwstats.formula = ~edges + nodemix("agecat2", base = 1) +
                nodefactor("deg.main") +
                concurrent+
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
dx.p

plot(dx.p)



# 3. Fit inst model ----------------------------------------------------------

#SAM ran:
load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/sj/fit.i.rda")



# Initialize network
nw.inst <- nw.main


# Assign degree -- assign main/casual degrees to inst network (because global constraint)
nw.inst <- set.vertex.attribute(nw.inst, "deg.main", nw.pers %v% "deg.main")
nw.inst <- set.vertex.attribute(nw.inst, "deg.pers", nw.main %v% "deg.pers")
table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers")


# Formulas
#not having age-homophily in inst model -- mainly because needed to run
#on core and only needed homophily to run the differential age*PT dissolutions

formation.i <- ~edges +
  nodefactor(c("deg.main", "deg.pers")) +
  nodefactor("agecat2") +
  nodefactor("riskg", base = 3) +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.i <- netest(nw.inst,
                formation = formation.i,
                target.stats = st$stats.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9, MCMLE.maxit = 250))

#Current error:
#Error : cannot allocate vector of size 4.5 Gb
#Error in fit$covar : $ operator is invalid for atomic vectors
#Solution: Sam ran on UW core

# Instantaneous Diagnostics -------------------------------------------------------------

dx.i <- netdx(fit.i, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
              nwstats.formula = ~edges +
                nodefactor(c("deg.main", "deg.pers")) +
                nodefactor("agecat2") +
                nodefactor("riskg", base = 3) +
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2)))

dx.i
plot(dx.i)




# Save data
est <- list(fit.m, fit.p, fit.i)
save(est, file = "C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/fit.rda")















#####################################################
## CONCURRENT BY AGE, STILL EGDES-ONLY DISSOLUTION ##
#####################################################

#
#
# # 1. Main Model -----------------------------------------------------------
#
#
# # Initialize network
# nw.main <- base_nw_msm(st)
#
#
# # Assign degree -- assign casual degree to main network (because global constraint)
# nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)
# #nw.main
# #save(nw.main, file = "C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nw.main.rda")
#
# # Formulas
# formation.m <- ~edges +
#   nodefactor("deg.pers") +
#   nodefactor("agecat2") +
#   absdiff("sqrt.age") +
#   offset(nodematch("role.class", diff = TRUE, keep = 1:2))
#
#
# # Fit model
# fit.m <- netest(nw.main,
#                 formation = formation.m,
#                 coef.form = c(-Inf, -Inf),
#                 target.stats = st$stats.m,
#                 coef.diss = st$coef.diss.m,
#                 constraints = ~bd(maxout = 1),
#                 set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e10, MCMLE.maxit = 250))
#
# fit.m
# summary(fit.m)
#
#
# #pre-age homophily/differential diss;
# dx <- netdx(fit.m, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
#             nwstats.formula = ~edges +
#               nodefactor("deg.pers") +
#               nodefactor("agecat2") +
#               absdiff("sqrt.age") +
#               offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
# dx
#
#
#
# # 2. Casual Model ---------------------------------------------------------
#
# # Initialize network
# nw.pers <- nw.main
#
# # Assign degree -- assign main degree to casual network (because global constraint)
# nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st)
# #nw.pers
# #save(nw.pers, file = "C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nw.pers.rda")
#
# # Formulas
# formation.p <- ~edges +
#               nodefactor("deg.main") +
#               nodefactor("agecat2") +
#               concurrent(by="agecat2") +
#               absdiff("sqrt.age") +
#               offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
#
#
# # Fit model -- with approximation
# fit.p <- netest(nw.pers,
#                 formation = formation.p,
#                 coef.form = c(-Inf, -Inf),
#                 target.stats = st$stats.p,
#                 coef.diss = st$coef.diss.p,
#                 constraints = ~bd(maxout = 2),
#                 set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9, MCMLE.maxit = 250))
#
#
# # Fit model -- without approximation
# fit.p <- netest(nw.pers,
#                 formation = formation.p,
#                 coef.form = c(-Inf, -Inf),
#                 target.stats = st$stats.p,
#                 coef.diss = st$coef.diss.p,
#                 constraints = ~bd(maxout = 2),
#                 set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9, MCMLE.maxit = 250),
#                 edapprox = FALSE)
#
# fit.p
# summary(fit.p)
#
#
#
# # Casual Diagnostics -------------------------------------------------------------
#
#
# ##Concurrent by age checks;
#
# dx.p <- netdx(fit.p, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
#               nwstats.formula = ~edges +
#                 nodefactor("deg.main") +
#                 nodefactor("agecat2") +
#                 concurrent(by="agecat2") +
#                 absdiff("sqrt.age") +
#                 offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
# dx.p
#
# plot(dx.p)
#
#
# #Plots and sim means are off from target stat for most of these;
# #Could be issue with tergm lite;
# #Try diagnostics on cross-section
#
# #concurrent check;
# #note test stats for concurrent:
# #conc.p.overall: 920.289
# #note test stats for concurrent(by="agecat2"):
# #conc.p.by.age[1:2]:  276.975        643.314
# #                     #Y w/ 2ConCasp #O w/ 2ConCasp
#
# dx.p2 <- netdx(fit.p, nsims = 10000, ncores = 1, dynamic = FALSE, #Try diagnostics on cross-section
#                nwstats.formula = ~edges +
#                  nodefactor("deg.main") +
#                  nodefactor("agecat2") +
#                  concurrent(by="agecat2") +
#                  absdiff("sqrt.age") +
#                  offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
# dx.p2
#
# plot(dx.p2)
#
#
#



