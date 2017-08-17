
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and age combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, disclosure of status, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module msm
#' @export
#'
condoms_msm <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    race <- dat$attr$race
    agecat2 <- dat$attr$agecat2

    # Parameters
    cond.rr.YY <- dat$param$cond.rr.YY
    cond.rr.OY <- dat$param$cond.rr.OY
    cond.rr.OO <- dat$param$cond.rr.OO
    # cond.rr.BB <- dat$param$cond.rr.BB
    # cond.rr.BW <- dat$param$cond.rr.BW
    # cond.rr.WW <- dat$param$cond.rr.WW

    if (type == "main") {
      cond.YY.prob <- dat$param$cond.main.YY.prob #0.5194
      cond.OY.prob <- dat$param$cond.main.OY.prob #0.3318
      cond.OO.prob <- dat$param$cond.main.OO.prob #0.3037
      #   cond.BB.prob <- dat$param$cond.main.BB.prob
      #   cond.BW.prob <- dat$param$cond.main.BW.prob
      #   cond.WW.prob <- dat$param$cond.main.WW.prob
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      cond.always <- NULL
      ptype <- 1
    }
    if (type == "pers") {
      cond.YY.prob <- dat$param$cond.pers.YY.prob #0.2923
      cond.OY.prob <- dat$param$cond.pers.OY.prob #0.3286
      cond.OO.prob <- dat$param$cond.pers.OO.prob #0.2363
      #   cond.BB.prob <- dat$param$cond.pers.BB.prob
      #   cond.BW.prob <- dat$param$cond.pers.BW.prob
      #   cond.WW.prob <- dat$param$cond.pers.WW.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      cond.always <- dat$attr$cond.always.pers
      ptype <- 2
    }
    if (type == "inst") {
      cond.YY.prob <- dat$param$cond.inst.YY.prob #0.2828
      cond.OY.prob <- dat$param$cond.inst.OY.prob #0.3381
      cond.OO.prob <- dat$param$cond.inst.OO.prob #0.2201
      #   cond.BB.prob <- dat$param$cond.inst.BB.prob
      #   cond.BW.prob <- dat$param$cond.inst.BW.prob
      #   cond.WW.prob <- dat$param$cond.inst.WW.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      cond.always <- dat$attr$cond.always.inst
      ptype <- 3
    }

    el <- dat$temp$el
    elt <- el[el[, "ptype"] == ptype, ]  #elt=edgelist by type, where t1=main, etc;

    ## Process ##

    # # Base condom probs
    # race.p1 <- race[elt[, 1]]  #condom prob in main partnerships;  #race?
    # race.p2 <- race[elt[, 2]]  #condom prob in casual partnerships;
    # num.B <- (race.p1 == "B") + (race.p2 == "B")
    # cond.prob <- (num.B == 2) * (cond.BB.prob * cond.rr.BB) +
    #              (num.B == 1) * (cond.BW.prob * cond.rr.BW) +
    #              (num.B == 0) * (cond.WW.prob * cond.rr.WW)

    # Base condom probs
    agecat2.p1 <- agecat2[elt[, 1]]  #taking partner 1
    agecat2.p2 <- agecat2[elt[, 2]]  #taking partner 2
    num.Y <- (agecat2.p1 == "Y") + (agecat2.p2 == "Y")  #concordant age edges
    cond.prob <- (num.Y == 2) * (cond.YY.prob * cond.rr.YY) + #if num.Y=2 then both YY -- all of these some value in params #0.5194 0.2923
                 (num.Y == 1) * (cond.OY.prob * cond.rr.OY) + #0.3318 0.3286
                 (num.Y == 0) * (cond.OO.prob * cond.rr.OO) #0.3037 0.2363


    # Transform to UAI logit
    uai.prob <- 1 - cond.prob
    uai.logodds <- log(uai.prob / (1 - uai.prob))

    # Diagnosis modifier
    pos.diag <- diag.status[elt[, 1]]   #el=p1, p2, st1, st2, ptype, ai
    isDx <- which(pos.diag == 1)        #so take those whose p1 is positive diag -- we did switch those pos to p1, but could be that p2 is also pos
                                        #and we don't know whether they are diagnosed...?
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta #uai prob for those who are diagnosed positive;

    # Disclosure modifier
    isDiscord <- which((elt[, "st1"] - elt[, "st2"]) == 1)
    delt <- elt[isDiscord, ]          #creating discordant el;
    discl.list <- dat$temp$discl.list  #list of those disclosed and at what time (first time is only 2)
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]  #creating some unique uid for dyad (XXXX000XXXX)
    discl.disc <- (delt.cdl %in% disclose.cdl)

    discl <- rep(NA, nrow(elt))
    discl[isDiscord] <- discl.disc

    isDisc <- which(discl == 1)
    uai.logodds[isDisc] <- uai.logodds[isDisc] + discl.beta

    # Back transform to prob
    old.uai.prob <- uai.prob
    uai.prob <- exp(uai.logodds) / (1 + exp(uai.logodds))

    uai.prob[is.na(uai.prob) & old.uai.prob == 0] <- 0
    uai.prob[is.na(uai.prob) & old.uai.prob == 1] <- 1

    # UAI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always[elt[, 1]]
      ca2 <- cond.always[elt[, 2]]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob)
      if (type == "pers") {
        dat$epi$cprob.always.pers[at] <- mean(uai.prob == 0) #% of dyads that had a partner with always cond use
      } else {
        dat$epi$cprob.always.inst[at] <- mean(uai.prob == 0)
      }
    }

    ai.vec <- elt[, "ai"]
    pos <- rep(elt[, "p1"], ai.vec)
    neg <- rep(elt[, "p2"], ai.vec)  #why does this assign neg to partners who may be positive but in st2 bc other partner st1 also positive?
    ptype <- rep(elt[, "ptype"], ai.vec)  #should be p1 and p2 not pos/neg since p2 could be pos if seroconcordant pos but still want uai prob for these;

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(pos), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)  #numbering out ids for all acts, repeat id for same pair multiple acts;
      al <- cbind(pos, neg, ptype, uai, pid) #[1620,] 5885 5874     1   1  902
                                             #[1621,] 5885 5874     1   0  902
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(pos, neg, ptype, uai, pid)
      al <- rbind(al, tmp.al)   #concatenate casuals to main for act list -- this is not necessarily discordant el -- just full right?;
    }

  } # end ptype loop

  dat$temp$al <- al

  return(dat)
}
