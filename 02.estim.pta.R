## Partner Type By Age ESTIM file ##

suppressMessages(library(EpiModelHIV)) 
rm(list = ls()) 

load("H:/a_PhD Second Year/EpiModel/EpiModelHIV-master/est/nwstats_pta.rda") 

# 1. Main Model ----------------------------------------------------------- 


# Initialize network 
nw.main <- base_nw_msm(st) 


# Assign degree -- assign casual degree to main network (because global constraint)
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st) 
 
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


# Fit model 
fit.m <- netest(nw.main, 
                formation = formation.m, 
                coef.form = c(-Inf, -Inf), 
                target.stats = st$stats.m, 
                coef.diss = st$coef.diss.m, 
                constraints = ~bd(maxout = 1), 
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e10, MCMLE.maxit = 250)) #normally 250


# 2. Casual Model --------------------------------------------------------- 
 

# Initialize network 
nw.pers <- nw.main 


# Assign degree -- assign main degree to casual network (because global constraint)
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st) 
 

# Formulas 
formation.p <- ~edges + 
                nodefactor("deg.main","agecat2") + 
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
                #nodefactor("deg.main","agecat2") +  ????
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
save(est, file = "H:/a_PhD Second Year/EpiModel/EpiModelHIV-master/est/fit_pta.rda") 


# Diagnostics ------------------------------------------------------------- 
# 
# dx <- netdx(est[[3]], nsims = 10000, ncores = 1, dynamic = FALSE) 
# dx 
