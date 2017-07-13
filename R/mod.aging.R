#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{age} and \code{sqrt.age}. The \code{sqrt.age} vertex attribute is also
#' updated on the three networks.
#'
#' @keywords module msm
#' @export
#'
aging_msm <- function(dat, at) {

  time.unit <- dat$param$time.unit

  age <- dat$attr$age
  active <- dat$attr$active
  agecat2 <- dat$attr$agecat2

  age[active == 1] <- age[active == 1] + time.unit / 365 #5292;


  ##########################################################################

  #pull new olds before changing them to “O”
  ids.get.old <- which(floor(age)==25 & agecat2=='Y' & active==1)

  nolds <- length(ids.get.old) #number of those getting old

  #for all new 25 year olds, going to assign them agecat2="O" -- will happen again and again for
  #25 yr olds through their first year, but oh well; could reassign every time step but that's annoying, agecat2[age >= 25] <- "O"
  #agecat2[floor(age)==25] <- "O"

  #or#agecat2[age>25 & age<=25.0192307692307692] <- "O"  #catches people in that first week

  agecat2[ids.get.old] <- "O"


  role.class <-dat$attr$role.class
  old.class <-dat$attr$role.class

  activenum<-sum(dat$attr$active==1)
  newolds <- rep(0, activenum)
  newolds[ids.get.old] <- 1
  #nolds<-sum(newolds)

  #before: role.class[ids.get.old]
  role.class[ids.get.old] <- sample(c("I", "R", "V"),
                                        nolds, replace = TRUE,
                                        prob = dat$param$role.O.prob)  #assigning role to new Olds in line w/ old role distribution;
#old.class[newolds==1]
# before  [1] "I" "V" "R" "V" "V" "I" "V" "V" "V" "R" "R" "R" "V" "I" "R"
# after   [1] "V" "V" "R" "V" "I" "I" "R" "V" "I" "R" "V" "I" "I" "R" "I"

  ins.quot <- dat$attr$ins.quot
  old.ins.quot <- dat$attr$ins.quot

  ins.quot[newolds==1 & role.class == "I"]  <- 1
  ins.quot[newolds==1 & role.class == "R"]  <- 0
  ins.quot[newolds==1 & role.class == "V"]  <- runif(sum(role.class[ids.get.old] == "V"))


  dat$attr$role.class <- role.class
  dat$attr$ins.quot <- ins.quot


  dat$attr$age <- age
  dat$attr$sqrt.age <- sqrt(age)
  dat$attr$agecat2 <- agecat2

  return(dat)
}


#' @title Aging Module
#'
#' @description This module ages all active nodes in the population by one time
#'              unit at each time step.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @keywords module het
#' @export
#'
aging_het <- function(dat, at) {

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Attributes
  age <- dat$attr$age
  active <- dat$attr$active

  ## Updates
  age[active == 1] <- age[active == 1] + time.unit/365

  ## Save out
  dat$attr$age <- age

  return(dat)
}
