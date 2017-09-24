
#' @title Births Module
#'
#' @description Module function for births or entries into the sexually active
#'              population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries among
#' black and white MSM, stochastically determined with draws from Poisson
#' distributions. For each new entry, a set of attributes is added for that node,
#' and the nodes are added onto the network objects. Only attributes that are
#' a part of the network model formulae are updated as vertex attributes on the
#' network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
births_msm <- function(dat, at){

  ## Variables

  # Parameters
  b.B.rate <- dat$param$b.B.rate
  b.W.rate <- dat$param$b.W.rate
  b.method <- dat$param$b.method


  ## Process
  if (b.method == "fixed") {
    numB <- dat$epi$num.B[1]
    numW <- dat$epi$num.W[1]
  }
  if (b.method == "varying") {
    numB <- dat$epi$num.B[at - 1]
    numW <- dat$epi$num.W[at - 1]
  }

  nBirths.B <- rpois(1, b.B.rate * numB) # tho these are separate, same rate of 0.001/7
  nBirths.W <- rpois(1, b.W.rate * numW)
  nBirths <- nBirths.B + nBirths.W

  ## Update Attr
  if (nBirths > 0) {
    dat <- setBirthAttr_msm(dat, at, nBirths.B, nBirths.W)  #runs function below, gives attr to new births;
  }


  # Update Networks
  if (nBirths > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- add_vertices(dat$el[[i]], nBirths)  #addes the new vertices into dataset
    }
  }


  ## Output
  dat$epi$nBirths[at] <- nBirths

  return(dat)
}


setBirthAttr_msm <- function(dat, at, nBirths.B, nBirths.W) {

  nBirths <- nBirths.B + nBirths.W

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nBirths) #uid doesn't get repeated then
  dat$temp$max.uid <- dat$temp$max.uid + nBirths         #uid always add on to n + births
                                                         #previous max (10000) + 8 new births last time + 14 births this time = 10022
                                                         #but newIDs used below -- just for setting initial attributes
                                                         #uid should be used for edgelists and referring to correct/same people
                                                         #no matter the timestep (may infect and then die, don't refer to new birth guy who took over that newID)

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)

  race <- sample(rep(c("B", "W"), c(nBirths.B, nBirths.W)))
  newB <- which(race == "B")
  newW <- which(race == "W")
  dat$attr$race[newIds] <- race

  dat$attr$age[newIds] <- rep(dat$param$birth.age, nBirths)
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds])


  # Disease status
  #dat$attr$status[newIds] <- rep(0, nBirths)

  dat$attr$status[newIds]<-rbinom(nBirths, 1, 0.07)

  dat$attr$inst.ai.class[newIds] <- sample(1:dat$param$num.inst.ai.classes,
                                           nBirths, replace = TRUE)

  dat$attr$tt.traj[newIds[newB]] <- sample(c(1, 2, 3, 4),
                                           nBirths.B, replace = TRUE,
                                           prob = dat$param$tt.traj.B.prob)
  dat$attr$tt.traj[newIds[newW]] <- sample(c(1, 2, 3, 4),
                                           nBirths.W, replace = TRUE,
                                           prob = dat$param$tt.traj.W.prob)

  ###############################################################################
  ## Infection-related attributes

  # #initialize all vectors present in "initialize module;
  # dat$attr$stage[newIds] <- rep(NA, nBirths)
  # dat$attr$stage.time[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.time[newIds] <- rep(NA, nBirths)
  # dat$attr$vl[newIds] <- rep(NA, nBirths)
  # dat$attr$diag.status[newIds] <- rep(NA, nBirths)
  # dat$attr$diag.time[newIds] <- rep(NA, nBirths)
  # dat$attr$tx.status[newIds] <- rep(NA, nBirths)
  # dat$attr$cum.time.on.tx[newIds] <- rep(NA, nBirths)
  # dat$attr$cum.time.off.tx[newIds] <- rep(NA, nBirths)
  # dat$attr$time.since.inf[newIds] <- rep(NA, nBirths)
  #
  # dat$attr$last.neg.test[newIds] <- rep(NA, nBirths)
  # dat$attr$tx.init.time[newIds] <- rep(NA, nBirths)
  # dat$attr$infector[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.role[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.type[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.diag[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.tx[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.stage[newIds] <- rep(NA, nBirths)
  # dat$attr$inf.agecat2[newIds] <- rep(NA, nBirths) ###added
  # dat$attr$infd.agecat2[newIds] <- rep(NA, nBirths) ###added -- should i add "Y"?


  # selected <- newIds[which(dat$attr$status[newIds] == 1)] #6, 6th position of newIds is TRUE
  selected <- intersect(which(dat$attr$status == 1), newIds)


  # attempt 1 -- newly infected but chronic stage (don't want all coming in as acute)
  # dat$attr$stage[newIdS[selected]] <- 3
  # dat$attr$stage.time[selected] <- 1
  # dat$attr$inf.time[selected] <- 1
  # dat$attr$vl[selected] <- dat$param$vl.set.point
  # dat$attr$diag.status[selected] <- 0 #rbinom(length(selected), 1, 0.60)
  # #dat$attr$diag.time[selected] <- rep(NA, nBirths)
  # dat$attr$tx.status[selected] <- 0
  # dat$attr$cum.time.on.tx[selected] <- rep(NA, nBirths)
  # dat$attr$cum.time.off.tx[selected] <- 1

  # attempt #2 -- coming in as chronic
  # dat$attr$stage[selected] <- 3 #chronic
  # dat$attr$stage.time[selected] <- 10 #person-time in  current HIV stage -- in days
  # dat$attr$inf.time[selected] <- at-100 ## if 1, is this going to give trouble with incoming at stage 3 but only given timestep infected?
  # dat$attr$vl[selected] <- dat$param$vl.set.point
  # dat$attr$diag.status[selected] <- rbinom(length(selected), 1, 0.60)
  # diagnosed <- which(dat$attr$status[selected] == 1 & dat$attr$diag.status[selected]==1)
  # dat$attr$diag.time[diagnosed] <- at
  #
  # dat$attr$tx.status[selected] <- 0
  # dat$attr$cum.time.on.tx[selected] <- 0  #treatment time to current time (0/NA bc not on tx) -- in trans newly infected are set to 0
  # dat$attr$cum.time.off.tx[selected] <- 100 #infect time to current time (-at)
  # dat$attr$time.since.inf[selected] <- 100

  #attempt 3 -- newly infected person;
  dat$attr$stage[selected]       <- 1
  dat$attr$stage.time[selected]  <- 0
  dat$attr$inf.time[selected]    <- at
  dat$attr$vl[selected]          <- 0
  dat$attr$diag.status[selected] <- 0
  dat$attr$tx.status[selected]   <- 0
  dat$attr$cum.time.on.tx[selected] <- 0
  dat$attr$cum.time.off.tx[selected] <- 0
  #no time since infection;

  ################################################################################
  # Circumcision
  dat$attr$circ[newIds[newB]] <- rbinom(nBirths.B, 1, dat$param$circ.B.prob)
  dat$attr$circ[newIds[newW]] <- rbinom(nBirths.W, 1, dat$param$circ.W.prob)



  # Role
  # dat$attr$role.class[newIds[newB]] <- sample(c("I", "R", "V"),
  #                                             nBirths.B, replace = TRUE,
  #                                             prob = dat$param$role.B.prob)
  # dat$attr$role.class[newIds[newW]] <- sample(c("I", "R", "V"),
  #                                             nBirths.W, replace = TRUE,
  #                                             prob = dat$param$role.W.prob)

  dat$attr$role.class[newIds] <- sample(c("I", "R", "V"),
                                        nBirths, replace = TRUE,
                                        prob = dat$param$role.Y.prob)  #assigning role to new Births in line w/ young role distribution;


  ins.quot <- rep(NA, nBirths)
  ins.quot[dat$attr$role.class[newIds] == "I"]  <- 1
  ins.quot[dat$attr$role.class[newIds] == "R"]  <- 0
  ins.quot[dat$attr$role.class[newIds] == "V"]  <-
                                  runif(sum(dat$attr$role.class[newIds] == "V"))
  dat$attr$ins.quot[newIds] <- ins.quot

  # CCR5 -- need to change to status=0 new IDS
  ccr5.B.prob <- dat$param$ccr5.B.prob
  ccr5.W.prob <- dat$param$ccr5.W.prob
  dat$attr$ccr5[newIds[newB]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B, replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.prob),
                                                 ccr5.B.prob[2], ccr5.B.prob[1]))
  dat$attr$ccr5[newIds[newW]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.W, replace = TRUE,
                                        prob = c(1 - sum(ccr5.W.prob),
                                                 ccr5.W.prob[2], ccr5.W.prob[1]))

  # Degree
  dat$attr$deg.main[newIds] <- 0
  dat$attr$deg.pers[newIds] <- 0

  # One-off risk group
  dat$attr$riskg[newIds] <- sample(1:5, nBirths, TRUE)

  # UAI group
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(nBirths, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers[newIds] <- uai.always[, 1]
  dat$attr$cond.always.inst[newIds] <- uai.always[, 2]

  # PrEP
  dat$attr$prepStat[newIds] <- 0

  # Agecat2
  dat$attr$agecat2[newIds] <- "Y"

  if (length(which(dat$attr$status == 1 & is.na(dat$attr$vl))) > 0) browser()

  return(dat)
}



#' @title Births Module
#'
#' @description Module for simulating births/entries into the population, including
#'              initialization of attributes for incoming nodes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
births_het <- function(dat, at) {

  # Variables
  b.rate.method <- dat$param$b.rate.method
  b.rate <- dat$param$b.rate
  active <- dat$attr$active


  # Process
  nBirths <- 0
  if (b.rate.method == "stgrowth") {
    exptPopSize <- dat$epi$num[1] * (1 + b.rate*at)
    numNeeded <- exptPopSize - sum(active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    }
  }
  if (b.rate.method == "totpop") {
    nElig <- dat$epi$num[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }
  if (b.rate.method == "fpop") {
    nElig <- dat$epi$num.feml[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }


  # Update Population Structure
  if (nBirths > 0) {
    dat <- setBirthAttr_het(dat, at, nBirths)
    dat$el <- add_vertices(dat$el, nBirths)
  }

  if (unique(sapply(dat$attr, length)) != attributes(dat$el)$n) {
    stop("mismatch between el and attr length in births mod")
  }

  # Output
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}


#' @title Assign Vertex Attributes at Network Entry
#'
#' @description Assigns vertex attributes to incoming nodes at birth/entry into
#'              the network.
#'
#' @inheritParams births_het
#' @param nBirths Number of new births as determined by \code{\link{births_het}}.
#'
#' @keywords het
#'
#' @export
#'
#'
setBirthAttr_het <- function(dat, at, nBirths) {

  # Set attributes for new births to NA
  dat$attr <- lapply(dat$attr, function(x) c(x, rep(NA, nBirths)))
  newIds <- which(is.na(dat$attr$active))


  # Network Status
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$entTime[newIds] <- rep(at, nBirths)


  # Demography
  prop.male <- ifelse(is.null(dat$param$b.propmale),
                      dat$epi$propMale[1],
                      dat$param$b.propmale)
  dat$attr$male[newIds] <- rbinom(nBirths, 1, prop.male)

  dat$attr$age[newIds] <- rep(18, nBirths)

  # Circumcision
  entTime <- dat$attr$entTime

  idsNewMale <- which(dat$attr$male == 1 & entTime == at)

  if (length(idsNewMale) > 0) {
    age <- dat$attr$age[idsNewMale]
    newCirc <- rbinom(length(idsNewMale), 1, dat$param$circ.prob.birth)
    isCirc <- which(newCirc == 1)

    newCircTime <- rep(NA, length(idsNewMale))
    newCircTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

    dat$attr$circStat[idsNewMale] <- newCirc
    dat$attr$circTime[idsNewMale] <- newCircTime
  }


  # Epi/Clinical
  dat$attr$status[newIds] <- rep(0, nBirths)

  if (length(unique(sapply(dat$attr, length))) != 1) {
    sapply(dat$attr, length)
    stop("Attribute dimensions not unique")
  }

  return(dat)
}
