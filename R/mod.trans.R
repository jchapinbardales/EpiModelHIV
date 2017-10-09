
# MSM -----------------------------------------------------------------

#' @title Transmission Module
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This is the final substantive function that occurs within the time loop at
#' each time step. This function takes the discordant edgelist and calculates a
#' transmission probability for each row (one sexual act) between dyads on the
#' network. After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' The per-act transmission probability depends on the following elements:
#' insertive versus receptive role, viral load of the infected partner, an
#' acute stage infection excess risk, condom use, and the CCR5 genetic allele.
#' Given these transmission probabilities, transmission is stochastically
#' simulating by drawing from a binomial distribution for each act conditional
#' on the per-act probability.
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race and age groups are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
trans_msm <- function(dat, at){

    # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  agecat2 <- dat$attr$agecat2 ###added,
  age <- dat$attr$age


  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
  acute.rr <- dat$param$acute.rr
  condom.rr <- dat$param$condom.rr
  circ.rr <- dat$param$circ.rr
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  prep.hr <- dat$param$prep.class.hr

  # Data
  dal <- dat$temp$dal
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]

  if (nrow(dal) == 0) {
    return(dat)
  }

  ## Reorder by role: ins on the left, rec on the right,
  ##                  with flippers represented twice
  disc.ip <- dal[dal[, "ins"] %in% 1:2, ]
  disc.rp <- dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- c("i", "r")
  colnames(disc.rp)[1:2] <- c("i", "r")


  # PATP: Insertive Man Infector (Col 1) --------------------------------

  # Attributes of infected

  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- URAI.prob * 2.45^(ip.vl - 4.5)

  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1-ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr)

  # CCR5
  ip.tlo[ip.ccr5 == "DD"] <- ip.tlo[ip.ccr5 == "DD"] + -Inf
  ip.tlo[ip.ccr5 == "DW"] <- ip.tlo[ip.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    ip.tlo[temp.ids] <- ip.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c(1, 2))
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)

  # Retransformation to probability
  ip.tprob <- plogis(ip.tlo)            # 8.938000e-03 0.000000e+00
  if (any(is.na(ip.tprob))) browser()
  if (any(ip.tprob < 0 | ip.tprob > 1)) browser()
  # stopifnot(ip.tprob >= 0, ip.tprob <= 1)


  # PATP: Receptive Man Infector (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible - no longer on right, reorder for ins, rec
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- UIAI.prob * 2.45^(rp.vl - 4.5)

  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1-rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr)

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr)

  # CCR5
  rp.tlo[rp.ccr5 == "DD"] <- rp.tlo[rp.ccr5 == "DD"] + -Inf
  rp.tlo[rp.ccr5 == "DW"] <- rp.tlo[rp.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    rp.tlo[temp.ids] <- rp.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c(1, 2))
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  # Retransformation to probability
  rp.tprob <- plogis(rp.tlo)
  if (any(is.na(rp.tprob))) browser()
  if (any(rp.tprob < 0 | rp.tprob > 1)) browser()
  # stopifnot(rp.tprob >= 0, rp.tprob <= 1)

  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob) #vector of 0s and 1's, binomial draw based on tp for insertive partner;
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob) #vector of 0s and 1's, binomial draw based on tp for receptive partner;


  # Output --------------------------------------------------------------

  # Update attributes

    infected <- infector <- inf.type <- inf.role <- inf.stage <- inf.condoms <-
      inf.diag <- inf.cum.time.on.tx <- inf.vl <- inf.agecat2 <- infd.agecat2 <- NULL

  if (sum(trans.ip, trans.rp) > 0) {

    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])
    infector <- c(disc.ip[trans.ip == 1, 1], #vector of uid's length of # who transmitted as insertive partner
                  disc.rp[trans.rp == 1, 2]) #vector of uid's length of # who transmitted as receptive partner

    inf.role <- c(rep(0, sum(trans.ip)), rep(1, sum(trans.rp))) #role of person who acquired (susceptible)
                #for as many who transmitted as insertive partner, put 0 for receptive role of person who acquired
                #for as many who transmitted as receptive partner, put 1 for insertive role of person who acquired
    inf.type <- c(disc.ip[trans.ip == 1, "ptype"],
                  disc.rp[trans.rp == 1, "ptype"])

    inf.stage <- stage[infector]
    inf.diag  <- diag.status[infector]
    inf.tx    <- tx.status[infector]
    inf.cum.time.on.tx <- dat$attr$cum.time.on.tx[infector]
    inf.vl <- vl[infector]

    browser()

    #age
    inf.agecat2 <- agecat2[infector]  #infectors age
    infd.agecat2 <- agecat2[infected]  #infected persons age
    inf.age <- age[infector]
    infd.age <- age[infected]



    dat$attr$status[infected]      <- 1
    dat$attr$inf.time[infected]    <- at
    dat$attr$vl[infected]          <- 0
    dat$attr$stage[infected]       <- 1
    dat$attr$stage.time[infected]  <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected]   <- 0

    dat$attr$infector[infected]  <- infector
    dat$attr$inf.role[infected]  <- inf.role #0=sus rec, 1=sus ins, so yes, role of infected/one who acquired
    dat$attr$inf.type[infected]  <- inf.type
    dat$attr$inf.diag[infected]  <- inf.diag
    dat$attr$inf.tx[infected]    <- inf.tx     #even though tx of infector, doesn't matter, define as tx of infector
                                            #for that infected/acquiring person - dont take literal as their tx status
    dat$attr$inf.stage[infected] <- inf.stage

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0

    #age
    dat$attr$inf.age[infected] <- inf.age  #give infected person value for age of infector
    dat$attr$infd.age[infected] <- infd.age  #give infected person value for age of infected
    dat$attr$inf.agecat2[infected] <- inf.agecat2  #give infected person value for agecat of infector
    dat$attr$infd.agecat2[infected] <- infd.agecat2  #give infected person value for agecat of infected
    #dont need to carry over in dat file to next simulation? but want to set this in epi data
    #to look at mean age of infection, mean age of transmitting overall at end, but don't need mean tracked [at] each timestep;


    # Summary Output
    dat$epi$incid[at] <- length(infected)

    dat$epi$trans.main[at] <- sum(inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.casl[at] <- sum(inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst[at] <- sum(inf.type == 3, na.rm = TRUE) / length(infected)

    dat$epi$trans.recpt.sus[at] <- sum(inf.role == 0, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus[at]  <- sum(inf.role == 1, na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act[at]  <- sum(inf.stage %in% 1:2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr[at]  <- sum(inf.stage == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids[at] <- sum(inf.stage == 4, na.rm = TRUE) / length(infected)

    #dat$epi$trans.condoms[at] <- sum(inf.condoms == 1, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx[at]         <- sum(inf.diag == 0, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                          inf.vl >= 4.5, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                          inf.vl < 4.5 & inf.vl > 1.5, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                          inf.vl <= 1.5, na.rm = TRUE) / length(infected)

      #among young infections -- total row

      dat$epi$trans.recpt.sus.amongY[at] <- sum(inf.role == 0 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.inst.sus.amongY[at]  <- sum(inf.role == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.stage.act.amongY[at]  <- sum(inf.stage %in% 1:2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.stage.chr.amongY[at]  <- sum(inf.stage == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.stage.aids.amongY[at] <- sum(inf.stage == 4 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.undx.amongY[at]         <- sum(inf.diag == 0 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.amongY[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notretained.amongY[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl >= 4.5 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.partsup.amongY[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl < 4.5 & inf.vl > 1.5 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.amongY[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl <= 1.5 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)


######################
#AGE SPECIFIC OUTPUTS;
######################

  #means/medians
    dat$epi$meanage.inf[at]  <- mean(inf.age, na.rm = TRUE)
    dat$epi$meanage.infd[at] <- mean(infd.age, na.rm = TRUE)
    dat$epi$medianage.inf[at]  <- median(inf.age, na.rm = TRUE)
    dat$epi$medianage.infd[at] <- median(infd.age, na.rm = TRUE)

  #incidence by age;
    dat$epi$incid.inf.Y[at]   <- sum(inf.agecat2 == "Y", na.rm = TRUE)
    dat$epi$incid.inf.O[at]   <- sum(inf.agecat2 == "O", na.rm = TRUE)
    dat$epi$incid.infd.Y[at]   <- sum(infd.agecat2 == "Y", na.rm = TRUE)
    dat$epi$incid.infd.O[at]   <- sum(infd.agecat2 == "O", na.rm = TRUE)

    dat$epi$incid.YY[at]  <- sum(inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE)
    dat$epi$incid.OY[at]  <- sum(sum(inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                 sum(inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE)
    dat$epi$incid.OO[at]  <- sum(inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE)
      #directional
      dat$epi$incid.OYd[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$incid.YOd[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "O", na.rm = TRUE)

  #age of infector
    dat$epi$trans.Y[at] <- sum(inf.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.O[at] <- sum(inf.agecat2 == "O", na.rm = TRUE) / length(infected)

  #age combo
    dat$epi$trans.YY[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.OY[at] <- sum(sum(inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                sum(inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.OO[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
        #directional
        dat$epi$trans.OYd[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
        dat$epi$trans.YOd[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)

        #among infections to young people
        dat$epi$trans.OYd.Y[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.YY.Y[at]  <- sum(inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

  #age of infector + PT
    dat$epi$trans.Ymain[at] <- sum(inf.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.Omain[at] <- sum(inf.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.Ycasl[at] <- sum(inf.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.Ocasl[at] <- sum(inf.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.Yinst[at] <- sum(inf.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.Oinst[at] <- sum(inf.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)

  #agecombo + PT
    dat$epi$trans.YYmain[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.OYmain[at] <- sum(sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1),
                                    sum(inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1), na.rm = TRUE) / length(infected)
    dat$epi$trans.OOmain[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)

    dat$epi$trans.YYcasl[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.OYcasl[at] <- sum(sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2),
                                    sum(inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2), na.rm = TRUE) / length(infected)
    dat$epi$trans.OOcasl[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)

    dat$epi$trans.YYinst[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.OYinst[at] <- sum(sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3),
                                    sum(inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3), na.rm = TRUE) / length(infected)
    dat$epi$trans.OOinst[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)

        #may want directional?;
        dat$epi$trans.OYdmain[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
        dat$epi$trans.YOdmain[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
        dat$epi$trans.OYdcasl[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
        dat$epi$trans.YOdcasl[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
        dat$epi$trans.OYdinst[at] <- sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
        dat$epi$trans.YOdinst[at] <- sum(inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)

          #among infections to young people
          dat$epi$trans.main.Y[at] <- sum(infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.casl.Y[at] <- sum(infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.inst.Y[at] <- sum(infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

          dat$epi$trans.OYdmain.Y[at] <-sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.YYmain.Y[at]  <-sum(inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.OYdcasl.Y[at] <-sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.YYcasl.Y[at]  <-sum(inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.OYdinst.Y[at] <-sum(inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.YYinst.Y[at]  <-sum(inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)



#####################################################################
#CONTINUOUS AGE & PT for GRAPHS - sum of infections, not PAFs;
#####################################################################

#cont age of infector ;
   dat$epi$trans.18[at] <- sum(floor(inf.age)==18, na.rm = TRUE)
   dat$epi$trans.19[at] <- sum(floor(inf.age)==19, na.rm = TRUE)
   dat$epi$trans.20[at] <- sum(floor(inf.age)==20, na.rm = TRUE)
   dat$epi$trans.21[at] <- sum(floor(inf.age)==21, na.rm = TRUE)
   dat$epi$trans.22[at] <- sum(floor(inf.age)==22, na.rm = TRUE)
   dat$epi$trans.23[at] <- sum(floor(inf.age)==23, na.rm = TRUE)
   dat$epi$trans.24[at] <- sum(floor(inf.age)==24, na.rm = TRUE)
   dat$epi$trans.25[at] <- sum(floor(inf.age)==25, na.rm = TRUE)
   dat$epi$trans.26[at] <- sum(floor(inf.age)==26, na.rm = TRUE)
   dat$epi$trans.27[at] <- sum(floor(inf.age)==27, na.rm = TRUE)
   dat$epi$trans.28[at] <- sum(floor(inf.age)==28, na.rm = TRUE)
   dat$epi$trans.29[at] <- sum(floor(inf.age)==29, na.rm = TRUE)
   dat$epi$trans.30[at] <- sum(floor(inf.age)==30, na.rm = TRUE)
   dat$epi$trans.31[at] <- sum(floor(inf.age)==31, na.rm = TRUE)
   dat$epi$trans.32[at] <- sum(floor(inf.age)==32, na.rm = TRUE)
   dat$epi$trans.33[at] <- sum(floor(inf.age)==33, na.rm = TRUE)
   dat$epi$trans.34[at] <- sum(floor(inf.age)==34, na.rm = TRUE)
   dat$epi$trans.35[at] <- sum(floor(inf.age)==35, na.rm = TRUE)
   dat$epi$trans.36[at] <- sum(floor(inf.age)==36, na.rm = TRUE)
   dat$epi$trans.37[at] <- sum(floor(inf.age)==37, na.rm = TRUE)
   dat$epi$trans.38[at] <- sum(floor(inf.age)==38, na.rm = TRUE)
   dat$epi$trans.39[at] <- sum(floor(inf.age)==39, na.rm = TRUE)
   dat$epi$trans.40[at] <- sum(floor(inf.age)==40, na.rm = TRUE)

#cont age of infector & PT ;
   #main;
   dat$epi$trans.18.main[at] <- sum(floor(inf.age)==18 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.19.main[at] <- sum(floor(inf.age)==19 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.20.main[at] <- sum(floor(inf.age)==20 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.21.main[at] <- sum(floor(inf.age)==21 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.22.main[at] <- sum(floor(inf.age)==22 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.23.main[at] <- sum(floor(inf.age)==23 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.24.main[at] <- sum(floor(inf.age)==24 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.25.main[at] <- sum(floor(inf.age)==25 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.26.main[at] <- sum(floor(inf.age)==26 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.27.main[at] <- sum(floor(inf.age)==27 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.28.main[at] <- sum(floor(inf.age)==28 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.29.main[at] <- sum(floor(inf.age)==29 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.30.main[at] <- sum(floor(inf.age)==30 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.31.main[at] <- sum(floor(inf.age)==31 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.32.main[at] <- sum(floor(inf.age)==32 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.33.main[at] <- sum(floor(inf.age)==33 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.34.main[at] <- sum(floor(inf.age)==34 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.35.main[at] <- sum(floor(inf.age)==35 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.36.main[at] <- sum(floor(inf.age)==36 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.37.main[at] <- sum(floor(inf.age)==37 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.38.main[at] <- sum(floor(inf.age)==38 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.39.main[at] <- sum(floor(inf.age)==39 & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.40.main[at] <- sum(floor(inf.age)==40 & inf.type == 1, na.rm = TRUE)

   #casl;
   dat$epi$trans.18.casl[at] <- sum(floor(inf.age)==18 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.19.casl[at] <- sum(floor(inf.age)==19 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.20.casl[at] <- sum(floor(inf.age)==20 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.21.casl[at] <- sum(floor(inf.age)==21 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.22.casl[at] <- sum(floor(inf.age)==22 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.23.casl[at] <- sum(floor(inf.age)==23 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.24.casl[at] <- sum(floor(inf.age)==24 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.25.casl[at] <- sum(floor(inf.age)==25 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.26.casl[at] <- sum(floor(inf.age)==26 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.27.casl[at] <- sum(floor(inf.age)==27 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.28.casl[at] <- sum(floor(inf.age)==28 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.29.casl[at] <- sum(floor(inf.age)==29 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.30.casl[at] <- sum(floor(inf.age)==30 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.31.casl[at] <- sum(floor(inf.age)==31 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.32.casl[at] <- sum(floor(inf.age)==32 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.33.casl[at] <- sum(floor(inf.age)==33 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.34.casl[at] <- sum(floor(inf.age)==34 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.35.casl[at] <- sum(floor(inf.age)==35 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.36.casl[at] <- sum(floor(inf.age)==36 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.37.casl[at] <- sum(floor(inf.age)==37 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.38.casl[at] <- sum(floor(inf.age)==38 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.39.casl[at] <- sum(floor(inf.age)==39 & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.40.casl[at] <- sum(floor(inf.age)==40 & inf.type == 2, na.rm = TRUE)

   #one off;
   dat$epi$trans.18.inst[at] <- sum(floor(inf.age)==18 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.19.inst[at] <- sum(floor(inf.age)==19 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.20.inst[at] <- sum(floor(inf.age)==20 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.21.inst[at] <- sum(floor(inf.age)==21 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.22.inst[at] <- sum(floor(inf.age)==22 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.23.inst[at] <- sum(floor(inf.age)==23 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.24.inst[at] <- sum(floor(inf.age)==24 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.25.inst[at] <- sum(floor(inf.age)==25 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.26.inst[at] <- sum(floor(inf.age)==26 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.27.inst[at] <- sum(floor(inf.age)==27 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.28.inst[at] <- sum(floor(inf.age)==28 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.29.inst[at] <- sum(floor(inf.age)==29 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.30.inst[at] <- sum(floor(inf.age)==30 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.31.inst[at] <- sum(floor(inf.age)==31 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.32.inst[at] <- sum(floor(inf.age)==32 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.33.inst[at] <- sum(floor(inf.age)==33 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.34.inst[at] <- sum(floor(inf.age)==34 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.35.inst[at] <- sum(floor(inf.age)==35 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.36.inst[at] <- sum(floor(inf.age)==36 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.37.inst[at] <- sum(floor(inf.age)==37 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.38.inst[at] <- sum(floor(inf.age)==38 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.39.inst[at] <- sum(floor(inf.age)==39 & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.40.inst[at] <- sum(floor(inf.age)==40 & inf.type == 3, na.rm = TRUE)

# to a young partner in X type;
   #main;
   dat$epi$trans.18.Ymain[at] <- sum(floor(inf.age)==18 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.19.Ymain[at] <- sum(floor(inf.age)==19 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.20.Ymain[at] <- sum(floor(inf.age)==20 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.21.Ymain[at] <- sum(floor(inf.age)==21 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.22.Ymain[at] <- sum(floor(inf.age)==22 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.23.Ymain[at] <- sum(floor(inf.age)==23 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.24.Ymain[at] <- sum(floor(inf.age)==24 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.25.Ymain[at] <- sum(floor(inf.age)==25 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.26.Ymain[at] <- sum(floor(inf.age)==26 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.27.Ymain[at] <- sum(floor(inf.age)==27 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.28.Ymain[at] <- sum(floor(inf.age)==28 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.29.Ymain[at] <- sum(floor(inf.age)==29 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.30.Ymain[at] <- sum(floor(inf.age)==30 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.31.Ymain[at] <- sum(floor(inf.age)==31 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.32.Ymain[at] <- sum(floor(inf.age)==32 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.33.Ymain[at] <- sum(floor(inf.age)==33 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.34.Ymain[at] <- sum(floor(inf.age)==34 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.35.Ymain[at] <- sum(floor(inf.age)==35 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.36.Ymain[at] <- sum(floor(inf.age)==36 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.37.Ymain[at] <- sum(floor(inf.age)==37 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.38.Ymain[at] <- sum(floor(inf.age)==38 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.39.Ymain[at] <- sum(floor(inf.age)==39 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.40.Ymain[at] <- sum(floor(inf.age)==40 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE)

   #casl;
   dat$epi$trans.18.Ycasl[at] <- sum(floor(inf.age)==18 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.19.Ycasl[at] <- sum(floor(inf.age)==19 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.20.Ycasl[at] <- sum(floor(inf.age)==20 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.21.Ycasl[at] <- sum(floor(inf.age)==21 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.22.Ycasl[at] <- sum(floor(inf.age)==22 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.23.Ycasl[at] <- sum(floor(inf.age)==23 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.24.Ycasl[at] <- sum(floor(inf.age)==24 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.25.Ycasl[at] <- sum(floor(inf.age)==25 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.26.Ycasl[at] <- sum(floor(inf.age)==26 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.27.Ycasl[at] <- sum(floor(inf.age)==27 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.28.Ycasl[at] <- sum(floor(inf.age)==28 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.29.Ycasl[at] <- sum(floor(inf.age)==29 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.30.Ycasl[at] <- sum(floor(inf.age)==30 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.31.Ycasl[at] <- sum(floor(inf.age)==31 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.32.Ycasl[at] <- sum(floor(inf.age)==32 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.33.Ycasl[at] <- sum(floor(inf.age)==33 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.34.Ycasl[at] <- sum(floor(inf.age)==34 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.35.Ycasl[at] <- sum(floor(inf.age)==35 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.36.Ycasl[at] <- sum(floor(inf.age)==36 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.37.Ycasl[at] <- sum(floor(inf.age)==37 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.38.Ycasl[at] <- sum(floor(inf.age)==38 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.39.Ycasl[at] <- sum(floor(inf.age)==39 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.40.Ycasl[at] <- sum(floor(inf.age)==40 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE)

   #one off;
   dat$epi$trans.18.Yinst[at] <- sum(floor(inf.age)==18 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.19.Yinst[at] <- sum(floor(inf.age)==19 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.20.Yinst[at] <- sum(floor(inf.age)==20 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.21.Yinst[at] <- sum(floor(inf.age)==21 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.22.Yinst[at] <- sum(floor(inf.age)==22 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.23.Yinst[at] <- sum(floor(inf.age)==23 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.24.Yinst[at] <- sum(floor(inf.age)==24 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.25.Yinst[at] <- sum(floor(inf.age)==25 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.26.Yinst[at] <- sum(floor(inf.age)==26 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.27.Yinst[at] <- sum(floor(inf.age)==27 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.28.Yinst[at] <- sum(floor(inf.age)==28 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.29.Yinst[at] <- sum(floor(inf.age)==29 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.30.Yinst[at] <- sum(floor(inf.age)==30 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.31.Yinst[at] <- sum(floor(inf.age)==31 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.32.Yinst[at] <- sum(floor(inf.age)==32 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.33.Yinst[at] <- sum(floor(inf.age)==33 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.34.Yinst[at] <- sum(floor(inf.age)==34 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.35.Yinst[at] <- sum(floor(inf.age)==35 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.36.Yinst[at] <- sum(floor(inf.age)==36 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.37.Yinst[at] <- sum(floor(inf.age)==37 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.38.Yinst[at] <- sum(floor(inf.age)==38 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.39.Yinst[at] <- sum(floor(inf.age)==39 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.40.Yinst[at] <- sum(floor(inf.age)==40 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE)

# to a older partner in X type;
   #main;
   dat$epi$trans.18.Omain[at] <- sum(floor(inf.age)==18 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.19.Omain[at] <- sum(floor(inf.age)==19 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.20.Omain[at] <- sum(floor(inf.age)==20 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.21.Omain[at] <- sum(floor(inf.age)==21 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.22.Omain[at] <- sum(floor(inf.age)==22 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.23.Omain[at] <- sum(floor(inf.age)==23 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.24.Omain[at] <- sum(floor(inf.age)==24 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.25.Omain[at] <- sum(floor(inf.age)==25 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.26.Omain[at] <- sum(floor(inf.age)==26 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.27.Omain[at] <- sum(floor(inf.age)==27 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.28.Omain[at] <- sum(floor(inf.age)==28 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.29.Omain[at] <- sum(floor(inf.age)==29 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.30.Omain[at] <- sum(floor(inf.age)==30 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.31.Omain[at] <- sum(floor(inf.age)==31 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.32.Omain[at] <- sum(floor(inf.age)==32 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.33.Omain[at] <- sum(floor(inf.age)==33 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.34.Omain[at] <- sum(floor(inf.age)==34 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.35.Omain[at] <- sum(floor(inf.age)==35 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.36.Omain[at] <- sum(floor(inf.age)==36 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.37.Omain[at] <- sum(floor(inf.age)==37 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.38.Omain[at] <- sum(floor(inf.age)==38 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.39.Omain[at] <- sum(floor(inf.age)==39 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)
   dat$epi$trans.40.Omain[at] <- sum(floor(inf.age)==40 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE)

   #casl;
   dat$epi$trans.18.Ocasl[at] <- sum(floor(inf.age)==18 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.19.Ocasl[at] <- sum(floor(inf.age)==19 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.20.Ocasl[at] <- sum(floor(inf.age)==20 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.21.Ocasl[at] <- sum(floor(inf.age)==21 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.22.Ocasl[at] <- sum(floor(inf.age)==22 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.23.Ocasl[at] <- sum(floor(inf.age)==23 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.24.Ocasl[at] <- sum(floor(inf.age)==24 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.25.Ocasl[at] <- sum(floor(inf.age)==25 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.26.Ocasl[at] <- sum(floor(inf.age)==26 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.27.Ocasl[at] <- sum(floor(inf.age)==27 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.28.Ocasl[at] <- sum(floor(inf.age)==28 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.29.Ocasl[at] <- sum(floor(inf.age)==29 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.30.Ocasl[at] <- sum(floor(inf.age)==30 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.31.Ocasl[at] <- sum(floor(inf.age)==31 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.32.Ocasl[at] <- sum(floor(inf.age)==32 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.33.Ocasl[at] <- sum(floor(inf.age)==33 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.34.Ocasl[at] <- sum(floor(inf.age)==34 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.35.Ocasl[at] <- sum(floor(inf.age)==35 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.36.Ocasl[at] <- sum(floor(inf.age)==36 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.37.Ocasl[at] <- sum(floor(inf.age)==37 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.38.Ocasl[at] <- sum(floor(inf.age)==38 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.39.Ocasl[at] <- sum(floor(inf.age)==39 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)
   dat$epi$trans.40.Ocasl[at] <- sum(floor(inf.age)==40 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE)

   #one off;
   dat$epi$trans.18.Oinst[at] <- sum(floor(inf.age)==18 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.19.Oinst[at] <- sum(floor(inf.age)==19 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.20.Oinst[at] <- sum(floor(inf.age)==20 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.21.Oinst[at] <- sum(floor(inf.age)==21 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.22.Oinst[at] <- sum(floor(inf.age)==22 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.23.Oinst[at] <- sum(floor(inf.age)==23 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.24.Oinst[at] <- sum(floor(inf.age)==24 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.25.Oinst[at] <- sum(floor(inf.age)==25 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.26.Oinst[at] <- sum(floor(inf.age)==26 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.27.Oinst[at] <- sum(floor(inf.age)==27 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.28.Oinst[at] <- sum(floor(inf.age)==28 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.29.Oinst[at] <- sum(floor(inf.age)==29 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.30.Oinst[at] <- sum(floor(inf.age)==30 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.31.Oinst[at] <- sum(floor(inf.age)==31 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.32.Oinst[at] <- sum(floor(inf.age)==32 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.33.Oinst[at] <- sum(floor(inf.age)==33 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.34.Oinst[at] <- sum(floor(inf.age)==34 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.35.Oinst[at] <- sum(floor(inf.age)==35 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.36.Oinst[at] <- sum(floor(inf.age)==36 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.37.Oinst[at] <- sum(floor(inf.age)==37 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.38.Oinst[at] <- sum(floor(inf.age)==38 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.39.Oinst[at] <- sum(floor(inf.age)==39 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)
   dat$epi$trans.40.Oinst[at] <- sum(floor(inf.age)==40 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE)



###############################################
#AGE & PT BY HIV STAGE OF INFECTION OF INFECTOR;
###############################################

  #age;
    dat$epi$trans.stage.act.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.O[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.O[at]  <- sum(inf.stage == 3 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.O[at] <- sum(inf.stage == 4 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
        #among young infections;
        dat$epi$trans.stage.act.Y.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.chr.Y.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.aids.Y.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.act.O.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.chr.O.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.aids.O.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

  #PT;
    dat$epi$trans.stage.act.main[at]  <- sum(inf.stage %in% 1:2 & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.main[at]  <- sum(inf.stage == 3 & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.main[at] <- sum(inf.stage == 4 & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.casl[at]  <- sum(inf.stage %in% 1:2 & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.casl[at]  <- sum(inf.stage == 3 & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.casl[at] <- sum(inf.stage == 4 & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.inst[at]  <- sum(inf.stage %in% 1:2 & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.inst[at]  <- sum(inf.stage == 3 & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.inst[at] <- sum(inf.stage == 4 & inf.type == 3, na.rm = TRUE) / length(infected)
        #among young infections;
        dat$epi$trans.stage.act.main.Y[at]  <- sum(inf.stage %in% 1:2 & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.chr.main.Y[at]  <- sum(inf.stage == 3 & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.aids.main.Y[at] <- sum(inf.stage == 4 & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.act.casl.Y[at]  <- sum(inf.stage %in% 1:2 & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.chr.casl.Y[at]  <- sum(inf.stage == 3 & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.aids.casl.Y[at] <- sum(inf.stage == 4 & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.act.inst.Y[at]  <- sum(inf.stage %in% 1:2 & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.chr.inst.Y[at]  <- sum(inf.stage == 3 & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.stage.aids.inst.Y[at] <- sum(inf.stage == 4 & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)


  #age combo;
    dat$epi$trans.stage.act.YY[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.YY[at]  <- sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.YY[at] <- sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act.OY[at]  <- sum(sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                           sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OY[at]  <- sum(sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                           sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OY[at] <- sum(sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                           sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act.OO[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OO[at]  <- sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OO[at] <- sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    #dont really care about directionality here? -- assume there are going to be more old chronics transmitting
    # to young;

  #age & PT;
    dat$epi$trans.stage.act.Ymain[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Ymain[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Ymain[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.Omain[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Omain[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Omain[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & inf.type == 1, na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act.Ycasl[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Ycasl[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Ycasl[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.Ocasl[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Ocasl[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Ocasl[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & inf.type == 2, na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act.Yinst[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Yinst[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Yinst[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.Oinst[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.Oinst[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.Oinst[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & inf.type == 3, na.rm = TRUE) / length(infected)

          #among young infections;
          dat$epi$trans.stage.act.Ymain.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.chr.Ymain.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.aids.Ymain.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.act.Omain.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.chr.Omain.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.aids.Omain.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & inf.type == 1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

          dat$epi$trans.stage.act.Ycasl.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.chr.Ycasl.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.aids.Ycasl.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.act.Ocasl.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.chr.Ocasl.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.aids.Ocasl.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & inf.type == 2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

          dat$epi$trans.stage.act.Yinst.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="Y" & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.chr.Yinst.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="Y" & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.aids.Yinst.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="Y" & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.act.Oinst.Y[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2=="O" & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.chr.Oinst.Y[at]  <- sum(inf.stage == 3 & inf.agecat2=="O" & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
          dat$epi$trans.stage.aids.Oinst.Y[at] <- sum(inf.stage == 4 & inf.agecat2=="O" & inf.type == 3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

  #age combo & PT;
    dat$epi$trans.stage.act.YYmain[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.YYmain[at]  <- sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.YYmain[at] <- sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.OYmain[at]  <- sum(sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1),
                                               sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OYmain[at]  <- sum(sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1),
                                               sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OYmain[at] <- sum(sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1),
                                               sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.OOmain[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OOmain[at]  <- sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OOmain[at] <- sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act.YYcasl[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.YYcasl[at]  <- sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.YYcasl[at] <- sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.OYcasl[at]  <- sum(sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2),
                                               sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OYcasl[at]  <- sum(sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2),
                                               sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OYcasl[at] <- sum(sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2),
                                               sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.OOcasl[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OOcasl[at]  <- sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OOcasl[at] <- sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)

    dat$epi$trans.stage.act.YYinst[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.YYinst[at]  <- sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.YYinst[at] <- sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.OYinst[at]  <- sum(sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3),
                                               sum(inf.stage %in% 1:2 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OYinst[at]  <- sum(sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3),
                                               sum(inf.stage == 3 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OYinst[at] <- sum(sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3),
                                               sum(inf.stage == 4 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3), na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.act.OOinst[at]  <- sum(inf.stage %in% 1:2 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.chr.OOinst[at]  <- sum(inf.stage == 3 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.stage.aids.OOinst[at] <- sum(inf.stage == 4 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
    #dont really care about directionality here? -- assume there are going to be more old chronics transmitting
    # to young;



############################
#AGE & PT BY HIV CARE STAGE;
############################

  #age
    dat$epi$trans.undx.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                            inf.vl >= 4.5 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                            inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl <= 1.5 & inf.agecat2=="Y", na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.O[at]         <- sum(inf.diag == 0 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.O[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.O[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                            inf.vl >= 4.5 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.O[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                            inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.O[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                            inf.vl <= 1.5 & inf.agecat2=="O", na.rm = TRUE) / length(infected)
        #among young infections;
        dat$epi$trans.undx.Y.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notinitiated.Y.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notretained.Y.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                  inf.vl >= 4.5 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.partsup.Y.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                  inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.fullsup.Y.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                  inf.vl <= 1.5 & inf.agecat2=="Y" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

        dat$epi$trans.undx.O.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notinitiated.O.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notretained.O.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                  inf.vl >= 4.5 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.partsup.O.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                  inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.fullsup.O.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                  inf.vl <= 1.5 & inf.agecat2=="O" & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)




    #PT
    dat$epi$trans.undx.main[at]         <- sum(inf.diag == 0 & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.main[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.main[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl >= 4.5 & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.main[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl < 4.5 & inf.vl > 1.5 & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.main[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl <= 1.5 & inf.type==1, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.casl[at]         <- sum(inf.diag == 0 & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.casl[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.casl[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl >= 4.5 & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.casl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl < 4.5 & inf.vl > 1.5 & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.casl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl <= 1.5 & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.inst[at]         <- sum(inf.diag == 0 & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.inst[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.inst[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl >= 4.5 & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.inst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl < 4.5 & inf.vl > 1.5 & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.inst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl <= 1.5 & inf.type==3, na.rm = TRUE) / length(infected)

        #among young infections
        dat$epi$trans.undx.main.Y[at]         <- sum(inf.diag == 0 & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notinitiated.main.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notretained.main.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl >= 4.5 & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.partsup.main.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl < 4.5 & inf.vl > 1.5 & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.fullsup.main.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl <= 1.5 & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

        dat$epi$trans.undx.casl.Y[at]         <- sum(inf.diag == 0 & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notinitiated.casl.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notretained.casl.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl >= 4.5 & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.partsup.casl.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl < 4.5 & inf.vl > 1.5 & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.fullsup.casl.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl <= 1.5 & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

        dat$epi$trans.undx.inst.Y[at]         <- sum(inf.diag == 0 & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notinitiated.inst.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.notretained.inst.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl >= 4.5 & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.partsup.inst.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl < 4.5 & inf.vl > 1.5 & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.fullsup.inst.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                     inf.vl <= 1.5 & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)




  #agecombo
    dat$epi$trans.undx.YY[at]         <- sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.YY[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.YY[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.YY[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.YY[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OY[at]         <- sum(sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                             sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OY[at] <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                             sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OY[at]  <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                             sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OY[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                             sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OY[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                             sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OO[at]         <- sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OO[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OO[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OO[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OO[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                             inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)

    #dont really care about directionality here -- assume there are going to be more old chronics transmitting
    # to young;

  #age & PT
    dat$epi$trans.undx.Ymain[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Ymain[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.Ymain[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl >= 4.5 & inf.agecat2=="Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.Ymain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Ymain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                               inf.vl <= 1.5 & inf.agecat2=="Y" & inf.type==1, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.Omain[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Omain[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.Omain[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl >= 4.5 & inf.agecat2=="O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.Omain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Omain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl <= 1.5 & inf.agecat2=="O" & inf.type==1, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.Ycasl[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Ycasl[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.Ycasl[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl >= 4.5 & inf.agecat2=="Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.Ycasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Ycasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl <= 1.5 & inf.agecat2=="Y" & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.Ocasl[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Ocasl[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.notretained.Ocasl[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl >= 4.5 & inf.agecat2=="O" & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.partsup.Ocasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Ocasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl <= 1.5 & inf.agecat2=="O" & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.Yinst[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Yinst[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.Yinst[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl >= 4.5 & inf.agecat2=="Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.Yinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Yinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl <= 1.5 & inf.agecat2=="Y" & inf.type==3, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.Oinst[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.Oinst[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.Oinst[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl >= 4.5 & inf.agecat2=="O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.Oinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.Oinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                inf.vl <= 1.5 & inf.agecat2=="O" & inf.type==3, na.rm = TRUE) / length(infected)

     #among young infections
      dat$epi$trans.undx.Ymain.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.Ymain.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notretained.Ymain.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl >= 4.5 & inf.agecat2=="Y" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.partsup.Ymain.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.Ymain.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl <= 1.5 & inf.agecat2=="Y" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.undx.Omain.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.Omain.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notretained.Omain.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl >= 4.5 & inf.agecat2=="O" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.partsup.Omain.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.Omain.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl <= 1.5 & inf.agecat2=="O" & inf.type==1 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.undx.Ycasl.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.Ycasl.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notretained.Ycasl.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl >= 4.5 & inf.agecat2=="Y" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.partsup.Ycasl.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.Ycasl.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl <= 1.5 & inf.agecat2=="Y" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.undx.Ocasl.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.Ocasl.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.notretained.Ocasl.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl >= 4.5 & inf.agecat2=="O" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.partsup.Ocasl.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.Ocasl.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl <= 1.5 & inf.agecat2=="O" & inf.type==2 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.undx.Yinst.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="Y" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.Yinst.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="Y" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notretained.Yinst.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl >= 4.5 & inf.agecat2=="Y" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.partsup.Yinst.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="Y" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.Yinst.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl <= 1.5 & inf.agecat2=="Y" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)

      dat$epi$trans.undx.Oinst.Y[at]         <- sum(inf.diag == 0 & inf.agecat2=="O" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notinitiated.Oinst.Y[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2=="O" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.notretained.Oinst.Y[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl >= 4.5 & inf.agecat2=="O" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.partsup.Oinst.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2=="O" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
      dat$epi$trans.fullsup.Oinst.Y[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                    inf.vl <= 1.5 & inf.agecat2=="O" & inf.type==3 & infd.agecat2=="Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)



  #age combo & PT
    dat$epi$trans.undx.YYmain[at]         <- sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.YYmain[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.YYmain[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.YYmain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.YYmain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==1, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OYmain[at]         <- sum(sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==1),
                                                 sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==1), na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OYmain[at] <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==1),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==1), na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OYmain[at]  <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==1),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==1), na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OYmain[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==1),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==1), na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OYmain[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==1),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==1), na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OOmain[at]         <- sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OOmain[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OOmain[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OOmain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==1, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OOmain[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==1, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.YYcasl[at]         <- sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.YYcasl[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.YYcasl[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.YYcasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.YYcasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OYcasl[at]         <- sum(sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==2),
                                                 sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==2), na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OYcasl[at] <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==2),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==2), na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OYcasl[at]  <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==2),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==2), na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OYcasl[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==2),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==2), na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OYcasl[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==2),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==2), na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OOcasl[at]         <- sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OOcasl[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OOcasl[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OOcasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==2, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OOcasl[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==2, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.YYinst[at]         <- sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.YYinst[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.YYinst[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.YYinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.YYinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type==3, na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OYinst[at]         <- sum(sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==3),
                                                 sum(inf.diag == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==3), na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OYinst[at] <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==3),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==3), na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OYinst[at]  <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==3),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl >= 4.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==3), na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OYinst[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==3),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==3), na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OYinst[at]      <- sum(sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type==3),
                                                 sum(inf.diag == 1 & inf.cum.time.on.tx > 0 & inf.vl <= 1.5 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type==3), na.rm = TRUE) / length(infected)

    dat$epi$trans.undx.OOinst[at]         <- sum(inf.diag == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notinitiated.OOinst[at] <- sum(inf.diag == 1 & inf.cum.time.on.tx == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.notretained.OOinst[at]  <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl >= 4.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.partsup.OOinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl < 4.5 & inf.vl > 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==3, na.rm = TRUE) / length(infected)
    dat$epi$trans.fullsup.OOinst[at]      <- sum(inf.diag == 1 & inf.cum.time.on.tx > 0 &
                                                 inf.vl <= 1.5 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type==3, na.rm = TRUE) / length(infected)

    #dont really care about directionality here -- assume there are going to be more old chronics transmitting
    # to young;


#########################
#Age & PT BY POSITIONING;
#########################

  #receptive by age -- note that NOT from infector perspective - want to highlight transmissions due to receptive sex;
                       #inf.role = role of person who acquired (susceptible)
    dat$epi$trans.recpt.sus.Y[at] <- sum(inf.role == 0 & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Y[at]  <- sum(inf.role == 1 & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.O[at] <- sum(inf.role == 0 & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.O[at]  <- sum(inf.role == 1 & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
        #among young infections;
        dat$epi$trans.recpt.sus.Y.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Y.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.O.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        #above = what proportions of transmissions to young pp occur during receptive sex with an older partner
        dat$epi$trans.inst.sus.O.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        #recall for nomenclature, above sus.O means older is susceptible but here sus.O.Y means Y is sus and O is transmitter

  #receptive by PT;
    dat$epi$trans.recpt.sus.main[at] <- sum(inf.role == 0 & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.main[at]  <- sum(inf.role == 1 & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.casl[at] <- sum(inf.role == 0 & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.casl[at]  <- sum(inf.role == 1 & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.inst[at] <- sum(inf.role == 0 & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.inst[at]  <- sum(inf.role == 1 & inf.type == 3, na.rm = TRUE) / length(infected)
        #among young infections;
        dat$epi$trans.recpt.sus.main.Y[at] <- sum(inf.role == 0 & inf.type == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.main.Y[at]  <- sum(inf.role == 1 & inf.type == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.casl.Y[at] <- sum(inf.role == 0 & inf.type == 2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.casl.Y[at]  <- sum(inf.role == 1 & inf.type == 2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.inst.Y[at] <- sum(inf.role == 0 & inf.type == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.inst.Y[at]  <- sum(inf.role == 1 & inf.type == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)


  #receptive by age combo;
    dat$epi$trans.recpt.sus.YY[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OY[at] <- sum(sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                          sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OO[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.YY[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OY[at]  <- sum(sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y"),
                                          sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O"), na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OO[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)

  #receptive by age & PT;
    dat$epi$trans.recpt.sus.Ymain[at] <- sum(inf.role == 0 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Ymain[at]  <- sum(inf.role == 1 & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.Omain[at] <- sum(inf.role == 0 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Omain[at]  <- sum(inf.role == 1 & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.Ycasl[at] <- sum(inf.role == 0 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Ycasl[at]  <- sum(inf.role == 1 & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.Ocasl[at] <- sum(inf.role == 0 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Ocasl[at]  <- sum(inf.role == 1 & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.Yinst[at] <- sum(inf.role == 0 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Yinst[at]  <- sum(inf.role == 1 & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.Oinst[at] <- sum(inf.role == 0 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.Oinst[at]  <- sum(inf.role == 1 & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
        #among young infections; sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.Ymain.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & inf.type == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Ymain.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & inf.type == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.Omain.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & inf.type == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Omain.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & inf.type == 1 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.Ycasl.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & inf.type == 2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Ycasl.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & inf.type == 2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.Ocasl.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & inf.type == 2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Ocasl.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & inf.type == 2 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.Yinst.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & inf.type == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Yinst.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & inf.type == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.recpt.sus.Oinst.Y[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & inf.type == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        dat$epi$trans.inst.sus.Oinst.Y[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & inf.type == 3 & infd.agecat2 == "Y", na.rm = TRUE) / sum(infd.agecat2 == "Y", na.rm = TRUE)
        #recall that before sus.Ymain meant Y was susceptible from main partner but here, we are saying Y susceptible from young main partner
        #or Omain.Y is that a Y susceptible from an Omain partner


  #receptive by age combo & PT;
    dat$epi$trans.recpt.sus.YYmain[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OYmain[at] <- sum(sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1),
                                              sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1), na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OOmain[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.YYmain[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OYmain[at]  <- sum(sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1),
                                              sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1), na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OOmain[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)

    dat$epi$trans.recpt.sus.YYcasl[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OYcasl[at] <- sum(sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2),
                                              sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2), na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OOcasl[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.YYcasl[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OYcasl[at]  <- sum(sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2),
                                             sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2), na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OOcasl[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)

    dat$epi$trans.recpt.sus.YYinst[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OYinst[at] <- sum(sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3),
                                              sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3), na.rm = TRUE) / length(infected)
    dat$epi$trans.recpt.sus.OOinst[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.YYinst[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OYinst[at]  <- sum(sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3),
                                              sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3), na.rm = TRUE) / length(infected)
    dat$epi$trans.inst.sus.OOinst[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)

        #may want directional?;
        dat$epi$trans.recpt.sus.OYd[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
        dat$epi$trans.recpt.sus.YOd[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.OYd[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y", na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.YOd[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O", na.rm = TRUE) / length(infected)

        dat$epi$trans.recpt.sus.OYdmain[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
        dat$epi$trans.recpt.sus.YOdmain[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.OYdmain[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 1, na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.YOdmain[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 1, na.rm = TRUE) / length(infected)
        dat$epi$trans.recpt.sus.OYdcasl[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
        dat$epi$trans.recpt.sus.YOdcasl[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.OYdcasl[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 2, na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.YOdcasl[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 2, na.rm = TRUE) / length(infected)
        dat$epi$trans.recpt.sus.OYdinst[at] <- sum(inf.role == 0 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
        dat$epi$trans.recpt.sus.YOdinst[at] <- sum(inf.role == 0 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.OYdinst[at]  <- sum(inf.role == 1 & inf.agecat2 == "O" & infd.agecat2 == "Y" & inf.type == 3, na.rm = TRUE) / length(infected)
        dat$epi$trans.inst.sus.YOdinst[at]  <- sum(inf.role == 1 & inf.agecat2 == "Y" & infd.agecat2 == "O" & inf.type == 3, na.rm = TRUE) / length(infected)


  }



    # Summary Output -- (if no transmissions)
    dat$epi$incid[at] <- length(infected)


    return(dat)
  }





  # HET -----------------------------------------------------------------

  #' @title Infection Module
  #'
  #' @description Module function to simulate transmission over an active
  #'              discordant edgelist.
  #'
  #' @inheritParams aging_het
  #'
  #' @keywords module het
  #'
  #' @export
  #'
  #'
  trans_het <- function(dat, at) {

    ## Discordant Edgelist
    del <- discord_edgelist_het(dat, at)

    nInf <- 0
    idsInf <- idsTrans <- NULL

    if (!is.null(del)) {

      ## Acts
      nedges <- length(del[[1]])

      act.rate.early <- dat$param$act.rate.early
      act.rate.late <- dat$param$act.rate.late
      act.rate.cd4 <- dat$param$act.rate.cd4

      cd4Count <- dat$attr$cd4Count[del$inf]

      isLate <- which(cd4Count < act.rate.cd4)

      rates <- rep(act.rate.early, nedges)
      rates[isLate] <- act.rate.late


      # Process
      act.rand <- dat$param$acts.rand
      if (act.rand == TRUE) {
        numActs <- rpois(nedges, rates)
      } else {
        numActs <- rates
      }

      cond.prob <- dat$param$cond.prob
      cond.prob <- rep(cond.prob, nedges)

      del$numActs <- numActs

      if (act.rand == TRUE) {
        del$protActs <- rbinom(nedges, rpois(nedges, numActs), cond.prob)
      } else {
        del$protActs <- numActs * cond.prob
      }

      del$protActs <- pmin(numActs, del$protActs)
      del$unprotActs <- numActs - del$protActs

      stopifnot(all(del$unprotActs >= 0))


      ## Transmission

      # Base transmission probability
      vlLevel <- dat$attr$vlLevel[del$inf]
      males <- dat$attr$male[del$sus]
      ages <- dat$attr$age[del$sus]
      circs <- dat$attr$circStat[del$sus]
      prop.male <- dat$epi$propMale[at - 1]
      base.tprob <- hughes_tp(vlLevel, males, ages, circs, prop.male)

      # Acute and aids stage multipliers
      acute.stage.mult <- dat$param$acute.stage.mult
      aids.stage.mult <- dat$param$aids.stage.mult

      isAcute <- which(at - dat$attr$infTime[del$inf] <
                         (dat$param$vl.acute.topeak + dat$param$vl.acute.toset))
      isAIDS <- which(dat$attr$cd4Count[del$inf] < 200)

      base.tprob[isAcute] <- base.tprob[isAcute] * acute.stage.mult
      base.tprob[isAIDS] <- base.tprob[isAIDS] * aids.stage.mult


      # Condoms
      # Probability as a mixture function of protected and unprotected acts
      cond.eff <- dat$param$cond.eff
      prob.stasis.protacts <- (1 - base.tprob*(1 - cond.eff)) ^ del$protActs
      prob.stasis.unptacts <- (1 - base.tprob) ^ del$unprotActs
      prob.stasis <- prob.stasis.protacts * prob.stasis.unptacts
      finl.tprob <- 1 - prob.stasis

      # Transmission
      del$base.tprob <- base.tprob
      del$finl.tprob <- finl.tprob

      stopifnot(length(unique(sapply(del, length))) == 1)

      # Random transmission given final trans prob
      idsTrans <- which(rbinom(nedges, 1, del$finl.tprob) == 1)

      # Subset discord edgelist to transmissions
      del <- keep.attr(del, idsTrans)


      ## Update Nodal Attr
      idsInf <- unique(del$sus)
      idsTrans <- unique(del$inf)
      nInf <- length(idsInf)

      if (nInf > 0) {
        dat$attr$status[idsInf] <- 1
        dat$attr$infTime[idsInf] <- at
        dat$attr$ageInf[idsInf] <- dat$attr$age[idsInf]
        dat$attr$dxStat[idsInf] <- 0
        dat$attr$vlLevel[idsInf] <- 0
        dat$attr$txCD4min[idsInf] <-
          pmin(rnbinom(nInf,
                       size = nbsdtosize(dat$param$tx.init.cd4.mean,
                                         dat$param$tx.init.cd4.sd),
                       mu = dat$param$tx.init.cd4.mean),
               dat$param$tx.elig.cd4)
      }

      ## Transmission data frame
      if (dat$control$save.transmat == TRUE) {
        if (nInf > 0) {
          if (at == 2) {
            dat$stats$transmat <- as.data.frame(del)
          } else {
            dat$stats$transmat <- rbind(dat$stats$transmat, as.data.frame(del))
          }
        }
      }

    }

    ## Incidence vector
    dat$epi$si.flow[at] <- nInf
    dat$epi$si.flow.male[at] <- sum(dat$attr$male[idsInf] == 1, na.rm = TRUE)
    dat$epi$si.flow.feml[at] <- sum(dat$attr$male[idsInf] == 0, na.rm = TRUE)

    return(dat)
  }


  discord_edgelist_het <- function(dat, at) {

    status <- dat$attr$status

    idsInft <- which(status == 1)
    nInft <- length(idsInft)

    del <- NULL

    if (nInft > 0) {

      if (is.null(dat$el)) {
        el <- get.dyads.active(dat$nw, at = at)
      } else {
        el <- dat$el
      }

      if (nrow(el) > 0) {
        el <- el[sample(1:nrow(el)), , drop = FALSE]

        disc <- which(abs(status[el[, 1]] - status[el[, 2]]) == 1)
        if (length(disc) > 0) {
          tmp.del <- el[disc, ]
          tmp.del[status[tmp.del[, 2]] == 1, ] <- tmp.del[status[tmp.del[, 2]] == 1, 2:1]

          del <- list()
          del$sus <- tmp.del[, 2]
          del$inf <- tmp.del[, 1]
        }
      }

    }

    return(del)
  }


  hughes_tp <- function(vls, susmales, susages, suscircs, prop.male, fmat = FALSE) {

    suscircs[is.na(suscircs)] <- 0

    sus.hsv2 <- 0.59*prop.male + 0.86*(1 - prop.male)
    sus.gud <- 0.039*prop.male + 0.053*(1 - prop.male)
    sus.tvagin <- 0.068*prop.male + 0.12*(1 - prop.male)
    sus.cerv <- 0.066*(1 - prop.male)

    interc <- -8.3067
    coef.vl <- 1.062566
    coef.male <- 0.6430989
    coef.age <- -0.0403451
    coef.hsv2 <- 0.7625081
    coef.circ <- -0.6377294
    coef.gud <- 0.9749536
    coef.vagin <- 0.9435334
    coef.cerv <- 1.288279

    tp.full <- exp(interc + coef.vl*(vls - 4) +
                     coef.male*susmales + coef.age*(susages - 35) +
                     coef.hsv2*sus.hsv2 + coef.circ*susmales*suscircs +
                     coef.gud*sus.gud + coef.vagin*sus.tvagin +
                     coef.cerv*sus.cerv)

    if (fmat == TRUE) {
      tp.full <- data.frame(tp.full, vls, susmales, susages, suscircs)
    }

    return(tp.full)
  }

