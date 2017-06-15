
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  race <- dat$attr$race
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  agecat2 <- dat$attr$agecat2

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA      #dat$epi for at=1, nsteps=52 (52 weeks in 1 year) so have 52 values for variables in dat$epi
    dat$epi$num.B <- rNA    #the first value is at initiation (initial num=10000, etc)
    dat$epi$num.W <- rNA    #had to leave num.B/W in because called into births_msm module
    dat$epi$num.Y <- rNA
    dat$epi$num.O <- rNA
    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    # dat$epi$i.num.B <- rNA
    # dat$epi$i.num.W <- rNA
    dat$epi$i.num.Y <- rNA    #i.num.Y = at step 1, at initiation, number of Y nodes that are infected;
    dat$epi$i.num.O <- rNA
    dat$epi$i.prev <- rNA
    # dat$epi$i.prev.B <- rNA
    # dat$epi$i.prev.W <- rNA
    dat$epi$i.prev.Y <- rNA
    dat$epi$i.prev.O <- rNA
    dat$epi$nBirths <- rNA
    dat$epi$dth.gen <- rNA
    dat$epi$dth.dis <- rNA
    dat$epi$incid <- rNA

    dat$epi$prepCurr <- rNA
    dat$epi$prepCov <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$incid.prep0 <- rNA
    dat$epi$incid.prep1 <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA

    dat$epi$cprob.always.pers <- rNA
    dat$epi$cprob.always.inst <- rNA

    #additional variables to initialize here for age*PT:
    dat$epi$trans.main <- rNA
    dat$epi$trans.casl <- rNA
    dat$epi$trans.inst <- rNA

    dat$epi$trans.recpt.sus <- rNA
    dat$epi$trans.inst.sus <- rNA

    dat$epi$trans.stage.act <- rNA
    dat$epi$trans.stage.chr <- rNA
    dat$epi$trans.stage.aids <- rNA

    #dat$epi$trans.condoms[at] <- rNA

    dat$epi$trans.undx <- rNA
    dat$epi$trans.notinitiated <- rNA
    dat$epi$trans.notretained <- rNA
    dat$epi$trans.partsup <- rNA
    dat$epi$trans.fullsup <- rNA


    ############################################################

    ######################
    #AGE SPECIFIC OUTPUTS;
    ######################

    #incidence by age;
      dat$epi$incid.inf.Y   <- rNA
      dat$epi$incid.inf.O   <- rNA
      dat$epi$incid.infd.Y  <- rNA
      dat$epi$incid.infd.O  <- rNA

      dat$epi$incid.YY  <- rNA
      dat$epi$incid.OY  <- rNA
      dat$epi$incid.OO  <- rNA
        #directional
        dat$epi$incid.OYd <- rNA
        dat$epi$incid.YOd <- rNA


    #age of infector
      dat$epi$trans.Y <- rNA
      dat$epi$trans.O <- rNA

    #age combo
      dat$epi$trans.YY <- rNA
      dat$epi$trans.OY <- rNA
      dat$epi$trans.OO <- rNA

    #age of infector + PT
      dat$epi$trans.Ymain <- rNA
      dat$epi$trans.Omain <- rNA
      dat$epi$trans.Ycasl <- rNA
      dat$epi$trans.Ocasl <- rNA
      dat$epi$trans.Yinst <- rNA
      dat$epi$trans.Oinst <- rNA


    #agecombo + PT
      dat$epi$trans.YYmain <- rNA
      dat$epi$trans.OYmain <- rNA
      dat$epi$trans.OOmain <- rNA
      dat$epi$trans.YYcasl <- rNA
      dat$epi$trans.OYcasl <- rNA
      dat$epi$trans.OOcasl <- rNA
      dat$epi$trans.YYinst <- rNA
      dat$epi$trans.OYinst <- rNA
      dat$epi$trans.OOinst <- rNA

        #may want directional?;
        dat$epi$trans.OYd <- rNA
        dat$epi$trans.YOd <- rNA
        dat$epi$trans.OYdmain <- rNA
        dat$epi$trans.YOdmain <- rNA
        dat$epi$trans.OYdcasl <- rNA
        dat$epi$trans.YOdcasl <- rNA
        dat$epi$trans.OYdinst <- rNA
        dat$epi$trans.YOdinst <- rNA



    ###############################################
    #AGE & PT BY HIV STAGE OF INFECTION OF INFECTOR;
    ###############################################

    #age;
      dat$epi$trans.stage.act.Y  <- rNA
      dat$epi$trans.stage.chr.Y  <- rNA
      dat$epi$trans.stage.aids.Y <- rNA
      dat$epi$trans.stage.act.O  <- rNA
      dat$epi$trans.stage.chr.O  <- rNA
      dat$epi$trans.stage.aids.O <- rNA

    #PT;
      dat$epi$trans.stage.act.main  <- rNA
      dat$epi$trans.stage.chr.main  <- rNA
      dat$epi$trans.stage.aids.main <- rNA
      dat$epi$trans.stage.act.casl  <- rNA
      dat$epi$trans.stage.chr.casl  <- rNA
      dat$epi$trans.stage.aids.casl <- rNA
      dat$epi$trans.stage.act.inst  <- rNA
      dat$epi$trans.stage.chr.inst  <- rNA
      dat$epi$trans.stage.aids.inst <- rNA

    #age combo;
      dat$epi$trans.stage.act.YY  <- rNA
      dat$epi$trans.stage.chr.YY  <- rNA
      dat$epi$trans.stage.aids.YY <- rNA

      dat$epi$trans.stage.act.OY  <- rNA
      dat$epi$trans.stage.chr.OY  <- rNA
      dat$epi$trans.stage.aids.OY <- rNA

      dat$epi$trans.stage.act.OO  <- rNA
      dat$epi$trans.stage.chr.OO  <- rNA
      dat$epi$trans.stage.aids.OO <- rNA

    #age & PT;
      dat$epi$trans.stage.act.Ymain  <- rNA
      dat$epi$trans.stage.chr.Ymain  <- rNA
      dat$epi$trans.stage.aids.Ymain <- rNA
      dat$epi$trans.stage.act.Omain  <- rNA
      dat$epi$trans.stage.chr.Omain  <- rNA
      dat$epi$trans.stage.aids.Omain <- rNA

      dat$epi$trans.stage.act.Ycasl  <- rNA
      dat$epi$trans.stage.chr.Ycasl  <- rNA
      dat$epi$trans.stage.aids.Ycasl <- rNA
      dat$epi$trans.stage.act.Ocasl  <- rNA
      dat$epi$trans.stage.chr.Ocasl  <- rNA
      dat$epi$trans.stage.aids.Ocasl <- rNA

      dat$epi$trans.stage.act.Yinst  <- rNA
      dat$epi$trans.stage.chr.Yinst  <- rNA
      dat$epi$trans.stage.aids.Yinst <- rNA
      dat$epi$trans.stage.act.Oinst  <- rNA
      dat$epi$trans.stage.chr.Oinst  <- rNA
      dat$epi$trans.stage.aids.Oinst <- rNA

    #age combo & PT;
      dat$epi$trans.stage.act.YYmain  <- rNA
      dat$epi$trans.stage.chr.YYmain  <- rNA
      dat$epi$trans.stage.aids.YYmain <- rNA
      dat$epi$trans.stage.act.OYmain  <- rNA
      dat$epi$trans.stage.chr.OYmain  <- rNA
      dat$epi$trans.stage.aids.OYmain <- rNA
      dat$epi$trans.stage.act.OOmain  <- rNA
      dat$epi$trans.stage.chr.OOmain  <- rNA
      dat$epi$trans.stage.aids.OOmain <- rNA

      dat$epi$trans.stage.act.YYcasl  <- rNA
      dat$epi$trans.stage.chr.YYcasl  <- rNA
      dat$epi$trans.stage.aids.YYcasl <- rNA
      dat$epi$trans.stage.act.OYcasl  <- rNA
      dat$epi$trans.stage.chr.OYcasl  <- rNA
      dat$epi$trans.stage.aids.OYcasl <- rNA
      dat$epi$trans.stage.act.OOcasl  <- rNA
      dat$epi$trans.stage.chr.OOcasl  <- rNA
      dat$epi$trans.stage.aids.OOcasl <- rNA

      dat$epi$trans.stage.act.YYinst  <- rNA
      dat$epi$trans.stage.chr.YYinst  <- rNA
      dat$epi$trans.stage.aids.YYinst <- rNA
      dat$epi$trans.stage.act.OYinst  <- rNA
      dat$epi$trans.stage.chr.OYinst  <- rNA
      dat$epi$trans.stage.aids.OYinst <- rNA
      dat$epi$trans.stage.act.OOinst  <- rNA
      dat$epi$trans.stage.chr.OOinst  <- rNA
      dat$epi$trans.stage.aids.OOinst <- rNA



    ############################
    #AGE & PT BY HIV CARE STAGE;
    ############################

    #age
      dat$epi$trans.undx.Y         <- rNA
      dat$epi$trans.notinitiated.Y <- rNA
      dat$epi$trans.notretained.Y  <- rNA
      dat$epi$trans.partsup.Y      <- rNA
      dat$epi$trans.fullsup.Y      <- rNA

      dat$epi$trans.undx.O         <- rNA
      dat$epi$trans.notinitiated.O <- rNA
      dat$epi$trans.notretained.O  <- rNA
      dat$epi$trans.partsup.O      <- rNA
      dat$epi$trans.fullsup.O      <- rNA

    #PT
      dat$epi$trans.undx.main         <- rNA
      dat$epi$trans.notinitiated.main <- rNA
      dat$epi$trans.notretained.main  <- rNA
      dat$epi$trans.partsup.main      <- rNA
      dat$epi$trans.fullsup.main      <- rNA

      dat$epi$trans.undx.casl        <- rNA
      dat$epi$trans.notinitiated.casl <- rNA
      dat$epi$trans.notretained.casl  <- rNA
      dat$epi$trans.partsup.casl      <- rNA
      dat$epi$trans.fullsup.casl      <- rNA

      dat$epi$trans.undx.inst         <- rNA
      dat$epi$trans.notinitiated.inst <- rNA
      dat$epi$trans.notretained.inst  <- rNA
      dat$epi$trans.partsup.inst      <- rNA
      dat$epi$trans.fullsup.inst      <- rNA

    #agecombo
      dat$epi$trans.undx.YY         <- rNA
      dat$epi$trans.notinitiated.YY <- rNA
      dat$epi$trans.notretained.YY  <- rNA
      dat$epi$trans.partsup.YY      <- rNA
      dat$epi$trans.fullsup.YY      <- rNA

      dat$epi$trans.undx.OY         <- rNA
      dat$epi$trans.notinitiated.OY <- rNA
      dat$epi$trans.notretained.OY  <- rNA
      dat$epi$trans.partsup.OY      <- rNA
      dat$epi$trans.fullsup.OY      <- rNA

      dat$epi$trans.undx.OO         <- rNA
      dat$epi$trans.notinitiated.OO <- rNA
      dat$epi$trans.notretained.OO  <- rNA
      dat$epi$trans.partsup.OO      <- rNA
      dat$epi$trans.fullsup.OO      <- rNA

    #age & PT
      dat$epi$trans.undx.Ymain         <- rNA
      dat$epi$trans.notinitiated.Ymain <- rNA
      dat$epi$trans.notretained.Ymain  <- rNA
      dat$epi$trans.partsup.Ymain      <- rNA
      dat$epi$trans.fullsup.Ymain      <- rNA

      dat$epi$trans.undx.Omain         <- rNA
      dat$epi$trans.notinitiated.Omain <- rNA
      dat$epi$trans.notretained.Omain  <- rNA
      dat$epi$trans.partsup.Omain      <- rNA
      dat$epi$trans.fullsup.Omain      <- rNA

      dat$epi$trans.undx.Ycasl         <- rNA
      dat$epi$trans.notinitiated.Ycasl <- rNA
      dat$epi$trans.notretained.Ycasl  <- rNA
      dat$epi$trans.partsup.Ycasl      <- rNA
      dat$epi$trans.fullsup.Ycasl      <- rNA

      dat$epi$trans.undx.Ocasl         <- rNA
      dat$epi$trans.notinitiated.Ocasl <- rNA
      dat$epi$trans.notretained.Ocasl  <- rNA
      dat$epi$trans.partsup.Ocasl      <- rNA
      dat$epi$trans.fullsup.Ocasl      <- rNA

      dat$epi$trans.undx.Yinst         <- rNA
      dat$epi$trans.notinitiated.Yinst <- rNA
      dat$epi$trans.notretained.Yinst  <- rNA
      dat$epi$trans.partsup.Yinst      <- rNA
      dat$epi$trans.fullsup.Yinst      <- rNA

      dat$epi$trans.undx.Oinst         <- rNA
      dat$epi$trans.notinitiated.Oinst <- rNA
      dat$epi$trans.notretained.Oinst  <- rNA
      dat$epi$trans.partsup.Oinst      <- rNA
      dat$epi$trans.fullsup.Oinst      <- rNA

    #age combo & PT
      dat$epi$trans.undx.YYmain         <- rNA
      dat$epi$trans.notinitiated.YYmain <- rNA
      dat$epi$trans.notretained.YYmain  <- rNA
      dat$epi$trans.partsup.YYmain      <- rNA
      dat$epi$trans.fullsup.YYmain      <- rNA

      dat$epi$trans.undx.OYmain         <- rNA
      dat$epi$trans.notinitiated.OYmain <- rNA
      dat$epi$trans.notretained.OYmain  <- rNA
      dat$epi$trans.partsup.OYmain      <- rNA
      dat$epi$trans.fullsup.OYmain      <- rNA

      dat$epi$trans.undx.OOmain         <- rNA
      dat$epi$trans.notinitiated.OOmain <- rNA
      dat$epi$trans.notretained.OOmain  <- rNA
      dat$epi$trans.partsup.OOmain      <- rNA
      dat$epi$trans.fullsup.OOmain      <- rNA

      dat$epi$trans.undx.YYcasl         <- rNA
      dat$epi$trans.notinitiated.YYcasl <- rNA
      dat$epi$trans.notretained.YYcasl  <- rNA
      dat$epi$trans.partsup.YYcasl      <- rNA
      dat$epi$trans.fullsup.YYcasl      <- rNA

      dat$epi$trans.undx.OYcasl         <- rNA
      dat$epi$trans.notinitiated.OYcasl <- rNA
      dat$epi$trans.notretained.OYcasl  <- rNA
      dat$epi$trans.partsup.OYcasl      <- rNA
      dat$epi$trans.fullsup.OYcasl      <- rNA

      dat$epi$trans.undx.OOcasl         <- rNA
      dat$epi$trans.notinitiated.OOcasl <- rNA
      dat$epi$trans.notretained.OOcasl  <- rNA
      dat$epi$trans.partsup.OOcasl      <- rNA
      dat$epi$trans.fullsup.OOcasl      <- rNA

      dat$epi$trans.undx.YYinst         <- rNA
      dat$epi$trans.notinitiated.YYinst <- rNA
      dat$epi$trans.notretained.YYinst  <- rNA
      dat$epi$trans.partsup.YYinst      <- rNA
      dat$epi$trans.fullsup.YYinst      <- rNA

      dat$epi$trans.undx.OYinst         <- rNA
      dat$epi$trans.notinitiated.OYinst <- rNA
      dat$epi$trans.notretained.OYinst  <- rNA
      dat$epi$trans.partsup.OYinst      <- rNA
      dat$epi$trans.fullsup.OYinst      <- rNA

      dat$epi$trans.undx.OOinst         <- rNA
      dat$epi$trans.notinitiated.OOinst <- rNA
      dat$epi$trans.notretained.OOinst  <- rNA
      dat$epi$trans.partsup.OOinst      <- rNA
      dat$epi$trans.fullsup.OOinst      <- rNA



    #########################
    #Age & PT BY POSITIONING;
    #########################

    #receptive by age;
      dat$epi$trans.recpt.sus.Y <- rNA
      dat$epi$trans.inst.sus.Y  <- rNA
      dat$epi$trans.recpt.sus.O <- rNA
      dat$epi$trans.inst.sus.O  <- rNA

    #receptive by PT;
      dat$epi$trans.recpt.sus.main <- rNA
      dat$epi$trans.inst.sus.main  <- rNA
      dat$epi$trans.recpt.sus.casl <- rNA
      dat$epi$trans.inst.sus.casl  <- rNA
      dat$epi$trans.recpt.sus.inst <- rNA
      dat$epi$trans.inst.sus.inst  <- rNA

    #receptive by age combo;
      dat$epi$trans.recpt.sus.YY <- rNA
      dat$epi$trans.recpt.sus.OY <- rNA
      dat$epi$trans.recpt.sus.OO <- rNA
      dat$epi$trans.inst.sus.YY  <- rNA
      dat$epi$trans.inst.sus.OY  <- rNA
      dat$epi$trans.inst.sus.OO  <- rNA
        #may want directional?;
        dat$epi$trans.recpt.sus.OYd <- rNA
        dat$epi$trans.recpt.sus.YOd <- rNA
        dat$epi$trans.inst.sus.OYd  <- rNA
        dat$epi$trans.inst.sus.YOd  <- rNA

    #receptive by age & PT;
      dat$epi$trans.recpt.sus.Ymain <- rNA
      dat$epi$trans.inst.sus.Ymain  <- rNA
      dat$epi$trans.recpt.sus.Omain <- rNA
      dat$epi$trans.inst.sus.Omain  <- rNA
      dat$epi$trans.recpt.sus.Ycasl <- rNA
      dat$epi$trans.inst.sus.Ycasl  <- rNA
      dat$epi$trans.recpt.sus.Ocasl <- rNA
      dat$epi$trans.inst.sus.Ocasl  <- rNA
      dat$epi$trans.recpt.sus.Yinst <- rNA
      dat$epi$trans.inst.sus.Yinst  <- rNA
      dat$epi$trans.recpt.sus.Oinst <- rNA
      dat$epi$trans.inst.sus.Oinst  <- rNA

    #receptive by age combo & PT;
      dat$epi$trans.recpt.sus.YYmain <- rNA
      dat$epi$trans.recpt.sus.OYmain <- rNA
      dat$epi$trans.recpt.sus.OOmain <- rNA
      dat$epi$trans.inst.sus.YYmain  <- rNA
      dat$epi$trans.inst.sus.OYmain  <- rNA
      dat$epi$trans.inst.sus.OOmain  <- rNA

      dat$epi$trans.recpt.sus.YYcasl <- rNA
      dat$epi$trans.recpt.sus.OYcasl <- rNA
      dat$epi$trans.recpt.sus.OOcasl <- rNA
      dat$epi$trans.inst.sus.YYcasl  <- rNA
      dat$epi$trans.inst.sus.OYcasl  <- rNA
      dat$epi$trans.inst.sus.OOcasl  <- rNA

      dat$epi$trans.recpt.sus.YYinst <- rNA
      dat$epi$trans.recpt.sus.OYinst <- rNA
      dat$epi$trans.recpt.sus.OOinst <- rNA
      dat$epi$trans.inst.sus.YYinst  <- rNA
      dat$epi$trans.inst.sus.OYinst  <- rNA
      dat$epi$trans.inst.sus.OOinst  <- rNA

          #may want directional?;
          dat$epi$trans.recpt.sus.OYdmain <- rNA
          dat$epi$trans.recpt.sus.YOdmain <- rNA
          dat$epi$trans.inst.sus.OYdmain  <- rNA
          dat$epi$trans.inst.sus.YOdmain  <- rNA
          dat$epi$trans.recpt.sus.OYdcasl <- rNA
          dat$epi$trans.recpt.sus.YOdcasl <- rNA
          dat$epi$trans.inst.sus.OYdcasl  <- rNA
          dat$epi$trans.inst.sus.YOdcasl  <- rNA
          dat$epi$trans.recpt.sus.OYdinst <- rNA
          dat$epi$trans.recpt.sus.YOdinst <- rNA
          dat$epi$trans.inst.sus.OYdinst  <- rNA
          dat$epi$trans.inst.sus.YOdinst  <- rNA


    ############################################################



    ##same as before below


  }


  dat$epi$num[at] <- length(status)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)
  dat$epi$num.Y[at] <- sum(agecat2 == "Y", na.rm = TRUE)
  dat$epi$num.O[at] <- sum(agecat2 == "O", na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)   #of infected nodes, at step 1 weighted avg of prev Y/O and pop Y/O = 2987;
  # dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  # dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.num.Y[at] <- sum(status == 1 & agecat2 == "Y", na.rm = TRUE)  #810 --> 810/3693=prev.Y
  dat$epi$i.num.O[at] <- sum(status == 1 & agecat2 == "O", na.rm = TRUE)  #2177 --> 2177/6307=prev.O
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  # dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  # dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]
  dat$epi$i.prev.Y[at] <- dat$epi$i.num.Y[at] / dat$epi$num.Y[at]
  dat$epi$i.prev.O[at] <- dat$epi$i.num.O[at] / dat$epi$num.O[at]

  #should I be taking prep stuff out?
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(dat$attr$prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) &
                                   status == 1, na.rm = TRUE)

  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
    sum((is.na(prepStat) | prepStat == 0), na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] /
      sum(prepStat == 1, na.rm = TRUE)
  }


  return(dat)
}





#' @title Prevalence Module
#'
#' @description Module function to calculate and store summary statistics for
#'              disease prevalence, demographics, and other epidemiological
#'              outcomes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
