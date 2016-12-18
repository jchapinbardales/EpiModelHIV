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
  nodefactor("agecat2") + 
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



# Main Diagnostics ------------------------------------------------------------- 


dx <- netdx(fit.m, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
            nwstats.formula = ~edges +
              nodefactor("deg.pers") + 
              nodefactor("agecat2") + 
              absdiff("sqrt.age") + 
              offset(nodematch("role.class", diff = TRUE, keep = 1:2)))
dx 


plot(dx)


# 2. Casual Model --------------------------------------------------------- 


# Initialize network 
nw.pers <- nw.main 


# Assign degree -- assign main degree to casual network (because global constraint)
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st) 


# Formulas 
formation.p <- ~edges + 
  nodefactor("deg.main") +
  nodefactor("agecat2") + 
  concurrent(by="agecat2") + 
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

fit.p


summary(fit.p)



# Casual Diagnostics ------------------------------------------------------------- 

dx.p <- netdx(fit.p, nsims = 10, nsteps=1000, ncores = 1, dynamic = TRUE,
              nwstats.formula = ~edges + 
                nodefactor("deg.main") +
                nodefactor("agecat2") +
                concurrent(by="agecat2") + 
                absdiff("sqrt.age") + 
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))) 
dx.p

plot(dx.p)

dx.p2 <- netdx(fit.p, nsims = 10000, ncores = 1, dynamic = FALSE,
              nwstats.formula = ~edges + 
                nodefactor("deg.main") +
                nodefactor("agecat2") +
                concurrent(by="agecat2") + 
                absdiff("sqrt.age") + 
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))) 
dx.p2

plot(dx.p2)


# 3. Fit inst model ---------------------------------------------------------- 


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

#### CURRENT ERROR ####
#   Error in vector.namesmatch(target.stats, names(nw.stats)) : 
#   Length of "target.stats" is 11 but should be 12.

## I’m thinking that this has to do with using nodefactor for deg,main
## in the formation model to create an interaction variable. 
## For some reason, it omits 1 deg of freedom per variable
## so 1df missing for deg.pers and 1df missing for deg.main
## Question is where do I create this variable -- in assign degree I assume?
## and create a deg.inst.age variable there while I'm at it for dx?


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
              nodefactor("deg.main","agecat2") +  #replace now with deg.main.age
              nodefactor("deg.pers","agecat2") +  #replace now with deg.pers.age
              nodefactor("deg.inst","agecat2") )  #need a deg.inst.age in assign_deg
dx 

st$stats.m
st$stats.p
st$stats.i