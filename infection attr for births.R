time.sex.active <- 1  #pmax(1, round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *  min(dat$init$ages), 0))
#always going to be only 1 week since entry to sex pool -- could set longer ago...?
vlar.int <- dat$param$vl.acute.rise.int
vlap <- dat$param$vl.acute.peak
vlaf.int <- dat$param$vl.acute.fall.int
vlsp <- dat$param$vl.set.point
vldo.int <- dat$param$vl.aids.onset.int
vl.aids.int <- dat$param$vl.aids.int
vlf  <- dat$param$vl.fatal
vlds <- (vlf - vlsp) / vl.aids.int
vl.acute.int <- vlar.int + vlaf.int


### Non-treater type: tester and non-tester
selected <- which(dat$attr$status[newIds] == 1 & dat$attr$tt.traj[newIds] %in% c(1, 2)) #6, 6th position of newIds is TRUE
max.inf.time <- pmin(time.sex.active[selected], vldo.int + vl.aids.int) # always going to be 1;
time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
dat$attr$inf.time[selected] <- 1 - time.since.inf
dat$attr$tx.status[selected] <- 0
dat$attr$cum.time.on.tx[selected] <- 0
dat$attr$cum.time.off.tx[selected] <- time.since.inf #0

dat$attr$stage[selected[time.since.inf <= vlar.int]] <- 1
dat$attr$stage[selected[time.since.inf > vlar.int & time.since.inf <= vl.acute.int]] <- 2
dat$attr$stage[selected[time.since.inf > vl.acute.int & time.since.inf <= vldo.int]] <- 3
dat$attr$stage[selected[time.since.inf > vldo.int]] <- 4

dat$attr$stage.time[selected][dat$attr$stage[selected] == 1] <- time.since.inf[dat$attr$stage[selected] == 1]
dat$attr$stage.time[selected][dat$attr$stage[selected] == 2] <- time.since.inf[dat$attr$stage[selected] == 2] -
  vlar.int
dat$attr$stage.time[selected][dat$attr$stage[selected] == 3] <- time.since.inf[dat$attr$stage[selected] == 3] -
  vl.acute.int
dat$attr$stage.time[selected][dat$attr$stage[selected] == 4] <- time.since.inf[dat$attr$stage[selected] == 4] -
  vldo.int

dat$attr$vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
  ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
  (time.since.inf > vlar.int + vlaf.int) * (time.since.inf <= vldo.int) * (vlsp) +
  (time.since.inf > vldo.int) * (vlsp + (time.since.inf - vldo.int) * vlds)

selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 1)
dat$attr$diag.status[selected] <- 0

selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 2)

# Time to next test
if (dat$param$testing.pattern == "interval") {
  ttntest <- ceiling(runif(length(selected),
                           min = 0,
                           max = dat$param$mean.test.B.int * (race[selected] == "B") +
                             dat$param$mean.test.W.int * (race[selected] == "W")))
}
if (dat$param$testing.pattern == "memoryless") {
  ttntest <- rgeom(length(selected),
                   1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
}

twind.int <- dat$param$test.window.int
dat$attr$diag.status[selected][ttntest > dat$attr$cum.time.off.tx[selected] - twind.int] <- 0
dat$attr$last.neg.test[selected][ttntest > dat$attr$cum.time.off.tx[selected] - twind.int] <-
  -ttntest[ttntest > dat$attr$cum.time.off.tx[selected] - twind.int]

dat$attr$diag.status[selected][ttntest <= dat$attr$cum.time.off.tx[selected] - twind.int] <- 1


### Full adherent type

# Create set of expected values for (cum.time.off.tx, cum.time.on.tx)

tx.init.time.B <- twind.int + dat$param$last.neg.test.B.int + 1 / dat$param$tx.init.B.prob
tx.init.time.W <- twind.int + dat$param$last.neg.test.W.int + 1 / dat$param$tx.init.W.prob

# Stage for Blacks
prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
  (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                  nrow = tx.init.time.B)
numsteps.B <- (dat$param$max.time.off.tx.full.int - tx.init.time.B) /
  (1 - prop.time.on.tx.B)
offon.B <- rbind(offon.B,
                 cbind(tx.init.time.B + (1 - prop.time.on.tx.B) * 1:numsteps.B,
                       prop.time.on.tx.B * 1:numsteps.B))
offon.B <- round(offon.B)
exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
exp.onset.aids.B <- nrow(offon.B)
offon.last.B <- offon.B[nrow(offon.B), ]
offon.B <- rbind(offon.B,
                 matrix(c(offon.last.B[1] + (1:vl.aids.int),
                          rep(offon.last.B[2], vl.aids.int)),
                        ncol = 2))
max.possible.inf.time.B <- nrow(offon.B)
offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

# Stage for Whites
prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
  (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                  nrow = tx.init.time.W)
numsteps.W <- (dat$param$max.time.off.tx.full.int - tx.init.time.W) /
  (1 - prop.time.on.tx.W)
offon.W <- rbind(offon.W,
                 cbind(tx.init.time.W + (1 - prop.time.on.tx.W) * 1:numsteps.W,
                       prop.time.on.tx.W * 1:numsteps.W))
offon.W <- round(offon.W)
exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
exp.onset.aids.W <- nrow(offon.W)
offon.last.W <- offon.W[nrow(offon.W), ]
offon.W <- rbind(offon.W,
                 matrix(c(offon.last.W[1] + (1:vl.aids.int),
                          rep(offon.last.W[2], vl.aids.int)),
                        ncol = 2))
max.possible.inf.time.W <- nrow(offon.W)
offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

# Vl for Blacks
selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 4 & race == "B")
max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
dat$attr$inf.time[selected] <- 1 - time.since.inf
dat$attr$cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
dat$attr$cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
dat$attr$stage[selected] <- stage.B[time.since.inf]
dat$attr$stage.time[selected] <- stage.time.B[time.since.inf]
dat$attr$tx.status[selected] <- 0
dat$attr$tx.status[selected][dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0] <-
  rbinom(sum(dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0),
         1, prop.time.on.tx.B)
dat$attr$vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
  ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
  (time.since.inf > vlar.int + vlaf.int) *
  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
  (time.since.inf > exp.onset.aids.B) *
  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
dat$attr$vl[selected][dat$attr$tx.status[selected] == 1] <- dat$param$vl.full.supp

# VL for Whites
selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 4 & race == "W")
max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
dat$attr$inf.time[selected] <- 1 - time.since.inf
dat$attr$cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
dat$attr$cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
dat$attr$stage[selected] <- stage.W[time.since.inf]
dat$attr$stage.time[selected] <- stage.time.W[time.since.inf]
dat$attr$tx.status[selected] <- 0
dat$attr$tx.status[selected][dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0] <-
  rbinom(sum(dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0),
         1, prop.time.on.tx.W)
dat$attr$vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
  ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
  (time.since.inf > vlar.int + vlaf.int) *
  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
  (time.since.inf > exp.onset.aids.W) *
  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
dat$attr$vl[selected][dat$attr$tx.status[selected] == 1] <- dat$param$vl.full.supp

# Diagnosis
selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 4)
if (dat$param$testing.pattern == "interval") {
  ttntest <- ceiling(runif(length(selected),
                           min = 0,
                           max = dat$param$mean.test.B.int * (race[selected] == "B") +
                             dat$param$mean.test.W.int * (race[selected] == "W")))
}
if (dat$param$testing.pattern == "memoryless") {
  ttntest <- rgeom(length(selected),
                   1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
}

dat$attr$diag.status[selected][ttntest > dat$attr$cum.time.off.tx[selected] - twind.int] <- 0
dat$attr$last.neg.test[selected][ttntest > dat$attr$cum.time.off.tx[selected] - twind.int] <-
  -ttntest[ttntest > dat$attr$cum.time.off.tx[selected] - twind.int]
dat$attr$diag.status[selected][ttntest <= dat$attr$cum.time.off.tx[selected] - twind.int] <- 1
dat$attr$diag.status[selected][dat$attr$cum.time.on.tx[selected] > 0] <- 1
dat$attr$last.neg.test[selected][dat$attr$cum.time.on.tx[selected] > 0] <- NA


### Part adherent type

# Create set of expected values for (cum.time.off.tx,cum.time.on.tx)

prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
  (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                  nrow = tx.init.time.B)
while (offon.B[nrow(offon.B), 1] / dat$param$max.time.off.tx.part.int +
       offon.B[nrow(offon.B), 2] / dat$param$max.time.on.tx.part.int < 1) {
  offon.B <- rbind(offon.B,
                   offon.B[nrow(offon.B), ] + c(1 - prop.time.on.tx.B,
                                                prop.time.on.tx.B))
}
offon.B <- round(offon.B)
exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
exp.onset.aids.B <- nrow(offon.B)
offon.last.B <- offon.B[nrow(offon.B), ]
offon.B <- rbind(offon.B,
                 matrix(c(offon.last.B[1] + (1:vl.aids.int),
                          rep(offon.last.B[2], vl.aids.int)),
                        ncol = 2))
max.possible.inf.time.B <- nrow(offon.B)
offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
  (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                  nrow = tx.init.time.W)

while (offon.W[nrow(offon.W), 1] / dat$param$max.time.off.tx.part.int +
       offon.W[nrow(offon.W), 2] / dat$param$max.time.on.tx.part.int < 1) {
  offon.W <- rbind(offon.W,
                   offon.W[nrow(offon.W), ] + c(1 - prop.time.on.tx.W,
                                                prop.time.on.tx.W))
}
offon.W <- round(offon.W)
exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
exp.onset.aids.W <- nrow(offon.W)
offon.last.W <- offon.W[nrow(offon.W), ]
offon.W <- rbind(offon.W,
                 matrix(c(offon.last.W[1] + (1:vl.aids.int),
                          rep(offon.last.W[2], vl.aids.int)),
                        ncol = 2))
max.possible.inf.time.W <- nrow(offon.W)
offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

# VL for Blacks
selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 3 & race == "B")
max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
dat$attr$inf.time[selected] <- 1 - time.since.inf
dat$attr$cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
dat$attr$cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
dat$attr$stage[selected] <- stage.B[time.since.inf]
dat$attr$stage.time[selected] <- stage.time.B[time.since.inf]
dat$attr$tx.status[selected] <- 0
dat$attr$tx.status[selected][dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0] <-
  rbinom(sum(dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0),
         1, prop.time.on.tx.B)
dat$attr$vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
  ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
  (time.since.inf > vlar.int + vlaf.int) *
  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
  (time.since.inf > exp.onset.aids.B) *
  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
dat$attr$vl[selected][dat$attr$tx.status[selected] == 1] <- dat$param$vl.part.supp

# VL for Whites
selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 3 & race == "W")
max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
dat$attr$inf.time[selected] <- 1 - time.since.inf
dat$attr$cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
dat$attr$cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
dat$attr$stage[selected] <- stage.W[time.since.inf]
dat$attr$stage.time[selected] <- stage.time.W[time.since.inf]
dat$attr$tx.status[selected] <- 0
dat$attr$tx.status[selected][dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0] <-
  rbinom(sum(dat$attr$stage[selected] == 3 & dat$attr$cum.time.on.tx[selected] >= 0),
         1, prop.time.on.tx.W)
dat$attr$vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
  ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
  (time.since.inf > vlar.int + vlaf.int) *
  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
  (time.since.inf > exp.onset.aids.W) *
  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
dat$attr$vl[selected][dat$attr$tx.status[selected] == 1] <- dat$param$vl.part.supp

# Implement diagnosis for both
selected <- which(dat$attr$status == 1 & dat$attr$tt.traj == 3)
if (dat$param$testing.pattern == "interval") {
  ttntest <- ceiling(runif(length(selected),
                           min = 0,
                           max = dat$param$mean.test.B.int * (race[selected] == "B") +
                             dat$param$mean.test.W.int * (race[selected] == "W")))
}

if (dat$param$testing.pattern == "memoryless") {
  ttntest <- rgeom(length(selected),
                   1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
}


dat$attr$diag.status[selected][ttntest > dat$attr$cum.time.off.tx[selected] - twind.int] <- 0
dat$attr$last.neg.test[selected][ttntest > dat$attr$cum.time.off.tx[selected] - twind.int] <-
  -ttntest[ttntest > dat$attr$cum.time.off.tx[selected] - twind.int]

dat$attr$diag.status[selected][ttntest <= dat$attr$cum.time.off.tx[selected] - twind.int] <- 1
dat$attr$diag.status[selected][dat$attr$cum.time.on.tx[selected] > 0] <- 1
dat$attr$last.neg.test[selected][dat$attr$cum.time.on.tx[selected] > 0] <- NA


# Last neg test before present for negatives
selected <- which(dat$attr$status == 0 & dat$attr$tt.traj %in% c(2, 3, 4))

if (dat$param$testing.pattern == "interval") {
  tslt <- ceiling(runif(length(selected),
                        min = 0,
                        max = dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
}
if (dat$param$testing.pattern == "memoryless") {
  tslt <- rgeom(length(selected),
                1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                       dat$param$mean.test.W.int * (race[selected] == "W")))
}
dat$attr$last.neg.test[selected] <- -tslt


# ## Set all onto dat$attr -- already did this along the way because was not recognizing object without dat$attr$
# dat$attr$stage <- stage
# dat$attr$stage.time <- stage.time
# dat$attr$inf.time <- inf.time
# dat$attr$vl <- vl
# dat$attr$diag.status <- diag.status
# dat$attr$diag.time <- diag.time         #NA
# dat$attr$last.neg.test <- last.neg.test
# dat$attr$tx.status <- tx.status
# dat$attr$tx.init.time <- tx.init.time   #NA, not race specific
# dat$attr$cum.time.on.tx <- cum.time.on.tx
# dat$attr$cum.time.off.tx <- cum.time.off.tx
# dat$attr$infector <- infector           #NA
# dat$attr$inf.role <- inf.role           #NA
# dat$attr$inf.type <- inf.type           #NA
# dat$attr$inf.diag <- inf.diag           #NA
# dat$attr$inf.tx <- inf.tx               #NA
# dat$attr$inf.stage <- inf.stage         #NA
# dat$attr$inf.agecat2 <- inf.agecat2 ###added, NA
# dat$attr$infd.agecat2 <- infd.agecat2 ###added, NA
