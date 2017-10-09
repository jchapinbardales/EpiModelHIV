devtools::install_github("statnet/EpiModel")
library(EpiModel)

# param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
# init <- init.icm(s.num = 500, i.num = 1)
# control <- control.icm(type = "SI", nsteps = 500, nsims = 20)
# mod1 <- icm(param, init, control)
# mdf <- as.data.frame(mod1, out = "vals")

library(ggplot2)
theme_set(theme_classic())
mdf <- as.data.frame(sim31d, out = "vals")

ggplot(mdf, aes(x = time, y = i.num, group = sim)) +
  geom_line(alpha = 0.5)

ggplot(mdf, aes(x = time)) +
  geom_line(aes(y = i.num, group = sim), alpha = 0.5,
            lwd = 0.25, color = "firebrick") +
  geom_bands(aes(y = i.num), lower = 0.025, upper = 0.975, fill = "firebrick") +
  geom_smooth(aes(y = i.num), lwd = 1, col = "firebrick")
