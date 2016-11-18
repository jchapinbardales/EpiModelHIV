## Partner Type By Age ESTIM file ##

suppressMessages(library(EpiModelHIV)) 
rm(list = ls()) 

load("C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nwstats.rda") 

# 1. Main Model ----------------------------------------------------------- 

st$stats.m
st$coef.diss.m

# Initialize network 
nw.main <- base_nw_msm(st) 


# Assign degree -- assign casual degree to main network (because global constraint)
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st) 
 nw.main
?assign_degree()
# Formulas 
formation.m <- ~edges + 
                nodefactor("deg.pers") + 
                #nodefactor("deg.pers","agecat2") + #deg.pers.age
                absdiff("sqrt.age") + 
                offset(nodematch("role.class", diff = TRUE, keep = 1:2)) 

#think what are the factors that go into the network formation -- the formation of main edges;
#the number of casual partners they have
#age homophily
#race - but not considering this really;
#age! - main degrees differ by age;
st$stats.m
st$deg.pers
nw.main$role.class

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
#
#Note: The constraint on the sample space is not dyad-independent. 
#Null model likelihood is only implemented for dyad-independent constraints at this time. 
#Number of observations is similarly ill-defined.
#OKAY? This is a dyad-dependent model, since the probability of edge formation between 
#any two nodes depends on the existence of edges between those nodes and other nodes.

# EpiModel Network Estimation
# =======================
#   Model class: netest
# Estimation Method: ERGM with Edges Approximation
# 
# Model Form
# -----------------------
#   Formation: ~edges + nodefactor("deg.pers") + absdiff("sqrt.age") + offset(nodematch("role.class", 
#                                                                                       diff = TRUE, keep = 1:2))
# Target Statistics: 1420 480 210 658.4067
# Constraints: ~bd(maxout = 1)
# 
# Dissolution: ~offset(edges)
# Target Statistics: 34.12442
# ==========================
#   Summary of model fit
# ==========================
#   
#   Formula:   nw ~ edges + nodefactor("deg.pers") + absdiff("sqrt.age") + offset(nodematch("role.class", 
#                                                                                           diff = TRUE, keep = 1:2))
# <environment: 0x0000000036c9b088>
#   
#   Iterations:  113 out of 250 
# 
# Monte Carlo MLE Results:
#   Estimate Std. Error MCMC % p-value    
# edges                  -8.69862    0.05636      0  <1e-04 ***
#   nodefactor.deg.pers.1  -0.44990    0.06039      0  <1e-04 ***
#   nodefactor.deg.pers.2  -0.37501    0.08386      0  <1e-04 ***
#   absdiff.sqrt.age       -1.20365    0.06866      0  <1e-04 ***
#   nodematch.role.class.I     -Inf    0.00000      0  <1e-04 ***
#   nodematch.role.class.R     -Inf    0.00000      0  <1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Log-likelihood was not estimated for this fit.
# To get deviances, AIC, and/or BIC from fit `object$fit` run 
# > object$fit<-logLik(object$fit, add=TRUE)
# to add it to the object or rerun this function with eval.loglik=TRUE.
# 
# The following terms are fixed by offset and are not estimated:
#   nodematch.role.class.I nodematch.role.class.R 
# 
# 
# Dissolution Coefficients
# =======================
#   Dissolution Model: ~offset(edges)
# Target Statistics: 34.12442
# Crude Coefficient: 3.500271
# Mortality/Exit Rate: 3.588345e-05
# Adjusted Coefficient: 3.502723


# Main Diagnostics ------------------------------------------------------------- 
# 
# dx <- netdx(fit.m, nsims = 10000, ncores = 1, dynamic = FALSE,
#             nwstats.formula = ~edges +
#               nodefactor("deg.pers") + 
#               #nodefactor("deg.pers","agecat2") + 
#               absdiff("sqrt.age") + 
#               offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
# dx 

# EpiModel Network Diagnostics
# =======================
#   Diagnostic Method: Static
# Simulations: 10000
# 
# Formation Diagnostics
# ----------------------- 
#   Target Sim Mean Pct Diff Sim SD
# edges                  1420.000 1417.410   -0.002 27.824
# nodefactor.deg.pers.1   480.000  478.232   -0.004 20.045
# nodefactor.deg.pers.2   210.000  209.146   -0.004 13.055
# absdiff.sqrt.age        658.407  658.475    0.000 19.684
# nodematch.role.class.I       NA    0.000       NA  0.000
# nodematch.role.class.R       NA    0.000       NA  0.000

# dx <- netdx(fit.m, nsims = 10000, ncores = 1, dynamic = FALSE,
#             nwstats.formula = ~edges +
#               nodefactor("deg.pers") + 
#               nodefactor("deg.pers","agecat2") + 
#               absdiff("sqrt.age") + 
#               offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
# dx 
#Error: argument number 2 to  model term is not of the expected numeric type
# I think can't have nodefactor with numeric and character var? but nodefactor for deg.pers
# earlier anyhow...?

dx <- netdx(fit.m, nsims = 10000, ncores = 1, dynamic = FALSE,
            nwstats.formula = ~edges +
              nodefactor("deg.pers") + 
              nodefactor("deg.pers.age") + 
              absdiff("sqrt.age") + 
              offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
dx 

# EpiModel Network Diagnostics
# =======================
#   Diagnostic Method: Static
# Simulations: 10000
# 
# Formation Diagnostics
# ----------------------- 
#                              Target Sim Mean Pct Diff Sim SD
# edges                      1420.000 1419.020   -0.001 27.212
# nodefactor.deg.pers.1       480.000  480.723    0.002 20.514
# nodefactor.deg.pers.2       210.000  210.214    0.001 12.903
# absdiff.sqrt.age            658.407  658.513    0.000 20.269
# nodematch.role.class.I           NA    0.000       NA  0.010
# nodematch.role.class.R           NA    0.000       NA  0.020
# nodefactor.deg.pers.age.O1       NA  302.487       NA 16.089
# nodefactor.deg.pers.age.O2       NA  131.067       NA 10.032
# nodefactor.deg.pers.age.Y0       NA  768.101       NA 26.108
# nodefactor.deg.pers.age.Y1       NA  178.236       NA 12.312
# nodefactor.deg.pers.age.Y2       NA   79.147       NA  7.852

#               745.986   166.185   77.553  1437.996  321.657  132.447
#               Y,0C      Y,1C      Y,2C    O,0C      O,1C     O,2C
#number of people with a main partner among those who are Y/O with X# casual partners
#target stats:  768.101   178.236   79.147  ref       302.487  131.067
#               Y,0C      Y,1C      Y,2C    O,0C      O,1C     O,2C

plot(dx)


# 2. Casual Model --------------------------------------------------------- 
 

# Initialize network 
nw.pers <- nw.main 


# Assign degree -- assign main degree to casual network (because global constraint)
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st) 
 

# Formulas 
formation.p <- ~edges + 
                #nodefactor("deg.main","agecat2") + 
                concurrent + 
                absdiff("sqrt.age") + 
                offset(nodematch("role.class", diff = TRUE, keep = 1:2)) 


# Fit model 
fit.p <- netest(nw.pers, 
                formation = formation.p, 
                coef.form = c(-Inf, -Inf), 
                target.stats = st$stats.p, 
                coef.diss = st$coef.diss.p, 
                constraints = ~bd(maxout = 2), 
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9, MCMLE.maxit = 250)) 


# Fit inst model ---------------------------------------------------------- 


# Initialize network 
nw.inst <- nw.main 
 

# Assign degree -- assign main/casual degrees to inst network (because global constraint)
nw.inst <- set.vertex.attribute(nw.inst, "deg.main", nw.pers %v% "deg.main") 
nw.inst <- set.vertex.attribute(nw.inst, "deg.pers", nw.main %v% "deg.pers") 
table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers") 


# Formulas 
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

# Save data 
est <- list(fit.m, fit.p, fit.i) 
save(est, file = "H:/a_PhD Second Year/EpiModel/EpiModelHIV/est/fit_pta.rda") 


# Diagnostics ------------------------------------------------------------- 
# 
 dx <- netdx(est[[3]], nsims = 10000, ncores = 1, dynamic = FALSE,
             nwstats.formula = ~edges +
               nodefactor(c("deg.main", "deg.pers")) + 
               absdiff("sqrt.age") +
               offset(nodematch("role.class", diff = TRUE, keep = 1:2)) +
               concurrent +
               nodefactor("riskg", base = 3) + 
               concurrent +
               nodefactor("deg.main","agecat2") +
               nodefactor("deg.pers","agecat2") +
               nodefactor("deg.inst","agecat2") )
 dx 
 
 st$stats.m
 st$stats.p
 st$stats.i
