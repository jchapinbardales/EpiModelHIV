
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the age of
#' both members in a pair and the expected values within young-young, old-young,
#' and old-old combinations. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @return
#' This function returns the \code{dat} object with the updated discordant act
#' list (\code{dal}). Each element of \code{dal} is a data frame with the ids of the
#' discordant pair repeated the number of times they have AI.
#'
#' @keywords module msm
#' @export
#'
acts_msm <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    status <- dat$attr$status
    race <- dat$attr$race
    agecat2 <- dat$attr$agecat2

    # Parameters
    ai.scale <- dat$param$ai.scale
    if (type == "main") {
      base.ai.YY.rate <- dat$param$base.ai.main.YY.rate
      base.ai.OY.rate <- dat$param$base.ai.main.OY.rate
      base.ai.OO.rate <- dat$param$base.ai.main.OO.rate
      # base.ai.BB.rate <- dat$param$base.ai.main.BB.rate
      # base.ai.BW.rate <- dat$param$base.ai.main.BW.rate
      # base.ai.WW.rate <- dat$param$base.ai.main.WW.rate
      fixed <- FALSE
      ptype <- 1
      el <- dat$el[[1]]
    }
    if (type == "pers") {
      base.ai.YY.rate <- dat$param$base.ai.pers.YY.rate
      base.ai.OY.rate <- dat$param$base.ai.pers.OY.rate
      base.ai.OO.rate <- dat$param$base.ai.pers.OO.rate
      # base.ai.BB.rate <- dat$param$base.ai.pers.BB.rate
      # base.ai.BW.rate <- dat$param$base.ai.pers.BW.rate
      # base.ai.WW.rate <- dat$param$base.ai.pers.WW.rate
      fixed <- FALSE
      ptype <- 2
      el <- dat$el[[2]]
    }
    if (type == "inst") {
      base.ai.YY.rate <- 1
      base.ai.OY.rate <- 1
      base.ai.OO.rate <- 1
      # base.ai.BB.rate <- 1
      # base.ai.BW.rate <- 1
      # base.ai.WW.rate <- 1
      fixed <- ifelse(ai.scale != 1, FALSE, TRUE)
      ptype <- 3
      el <- dat$el[[3]]
    }

    ## Processes ##

    # Construct edgelist

    st1 <- status[el[, 1]]
    st2 <- status[el[, 2]]
    disc <- abs(st1 - st2) == 1 #discordant ;
    el[which(disc == 1 & st2 == 1), ] <- el[which(disc == 1 & st2 == 1), 2:1] #switches second column to first for those disc pairs that the st2 position is the positive one
    el <- cbind(el, status[el[, 1]], status[el[, 2]])                         #this makes is so either: 0 0, 1 0, 1 1 -- then subtract st1-st2=1 to get disc edgelist later (condoms)
    colnames(el) <- c("p1", "p2", "st1", "st2")

    if (nrow(el) > 0) {

      # # Base AI rates
      # ai.rate <- rep(NA, nrow(el))
      # race.p1 <- race[el[, 1]]
      # race.p2 <- race[el[, 2]]
      # num.B <- (race.p1 == "B") + (race.p2 == "B")
      # ai.rate <- (num.B == 2) * base.ai.BB.rate +
      #            (num.B == 1) * base.ai.BW.rate +
      #            (num.B == 0) * base.ai.WW.rate
      # ai.rate <- ai.rate * ai.scale

      # Base AI rates
      ai.rate <- rep(NA, nrow(el))
      agecat2.p1 <- agecat2[el[, 1]]
      agecat2.p2 <- agecat2[el[, 2]]
      num.Y <- (agecat2.p1 == "Y") + (agecat2.p2 == "Y")  #does this create prob w/ num.Y= num of y nodes? no bc we dont set it to dat
      ai.rate <- (num.Y == 2) * base.ai.YY.rate +    #YY 1.0955 0.7273
        (num.Y == 1) * base.ai.OY.rate +    #OY 1.1466 0.7609
        (num.Y == 0) * base.ai.OO.rate      #OO 1.4280 0.9478
      ai.rate <- ai.rate * ai.scale

      # Final act number
      if (fixed == FALSE) {
        ai <- rpois(length(ai.rate), ai.rate)
      } else {
        ai <- round(ai.rate)
      }

      # Full edge list
      el <- cbind(el, ptype, ai)
      colnames(el)[5:6] <- c("ptype", "ai")

      if (type == "main") {
        dat$temp$el <- el
      } else {
        dat$temp$el <- rbind(dat$temp$el, el) #concatenates other el from cas, inst
      }
    }

  } # loop over type end

  # Remove inactive edges (acts = 0 this timestep) from el
  dat$temp$el <- dat$temp$el[-which(dat$temp$el[, "ai"] == 0), ]

  return(dat)
}

