
devtools::install_github("statnet/tergmLite")
devtools::install_github("statnet/EpiModelHIV")

system("scp *.* libra:~/jcb")
setwd("~/Downloads")

library(EpiModelHIV)
load("st.rda")
load("nw.main.rda")
load("nw.pers.rda")

# Initialize network
nw.inst <- nw.main

# Assign degree
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

save(fit.i, file = "fit.i.rda")
system("scp libra:~/jcb/fit.i.rda .")

dx.i <- netdx(fit.i, nsims = 1e5, dynamic = FALSE)
dx.i

