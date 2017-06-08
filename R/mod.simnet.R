
# MSM -----------------------------------------------------------------

#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the main, casual, and one-off
#'              networks for one time step.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
simnet_msm <- function(dat, at) {

  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)

  #if (dat$param$method == 1) {
  #  dat$attr$deg.pers <- get_degree(dat$el[[2]])
  #} else {
    dat$attr$deg.pers <- paste0(dat$attr$agecat2, get_degree(dat$el[[2]])) #changed race to agecat2
  #}
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam.m$coef.form,
                                             coef.diss = nwparam.m$coef.diss$coef.adj,
                                             save.changes = TRUE)

  dat$temp$new.edges <- NULL
  if (at == 2) {
    new.edges.m <- matrix(dat$el[[1]], ncol = 2)
  } else {
    new.edges.m <- attributes(dat$el[[1]])$changes
    new.edges.m <- new.edges.m[new.edges.m[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- matrix(dat$attr$uid[new.edges.m], ncol = 2)  #how does this organize the edges? reorders them bc last one not same...


  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)

  # if (dat$param$method == 1) {
  #   dat$attr$deg.main <- get_degree(dat$el[[1]])
  # } else {
    dat$attr$deg.main <- paste0(dat$attr$agecat2, get_degree(dat$el[[1]]))
  #}
  dat <- tergmLite::updateModelTermInputs(dat, network = 2)

  dat$el[[2]] <- tergmLite::simulate_network(p = dat$p[[2]],
                                             el = dat$el[[2]],
                                             coef.form = nwparam.p$coef.form,
                                             coef.diss = nwparam.p$coef.diss$coef.adj,
                                             save.changes = TRUE)

  if (at == 2) {
    new.edges.p <- matrix(dat$el[[2]], ncol = 2)
  } else {
    new.edges.p <- attributes(dat$el[[2]])$changes
    new.edges.p <- new.edges.p[new.edges.p[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- rbind(dat$temp$new.edges,
                              matrix(dat$attr$uid[new.edges.p], ncol = 2))


  ## One-off network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)

  #if (dat$param$method == 1) {
  #  dat$attr$deg.pers <- get_degree(dat$el[[2]])
  #} else {
    dat$attr$deg.pers <- paste0(dat$attr$agecat2, get_degree(dat$el[[2]]))
  #}
  dat <- tergmLite::updateModelTermInputs(dat, network = 3)

  dat$el[[3]] <- tergmLite::simulate_ergm(p = dat$p[[3]],
                                          el = dat$el[[3]],
                                          coef = nwparam.i$coef.form)

  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_resim_nwstats(dat, at)
  }

  return(dat)
}



calc_resim_nwstats <- function(dat, at) {

  agecat2 <- dat$attr$agecat2

  for (nw in 1:3) {

    nwnum <- attr(dat$el[[nw]], "n")

    agecat2[dat$el[[nw]]]
    agecat2.el <- array(agecat2[dat$el[[nw]]], dim = dim(dat$el[[nw]]))

    nm.OO <- sum(ifelse(agecat2.el[, 1] == "O" & agecat2.el[, 2] == "O", 1, 0))
    nm.YY <- sum(ifelse(agecat2.el[, 1] == "Y" & agecat2.el[, 2] == "Y", 1, 0))
    nm.OY <- nrow(dat$el[[nw]]) - nm.YY - nm.OO
    agedisc <- nm.OY/nrow(dat$el[[nw]])

    #  nodematch <- sum(ifelse(age.el[, 1] == age.el[, 2], 1, 0))

    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges / nwnum, 3)
    concurrent <- round(mean(get_degree(dat$el[[nw]]) > 1), 3)

    mat <- matrix(c(edges, meandeg, concurrent, nm.OO, nm.OY, nm.YY, agedisc), ncol = 7, nrow = 1)



    if (at == 2) {
      dat$stats$nwstats[[nw]] <- mat
      colnames(dat$stats$nwstats[[nw]]) <- c("edges", "meand", "conc", "nm.OO", "nm.OY", "nm.YY", "agedisc")


    }
    if (at > 2) {
      dat$stats$nwstats[[nw]] <- rbind(dat$stats$nwstats[[nw]], mat)
    }
  }

  return(dat)
}




#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficients in a dynamic network model
#'              to preserve the mean degree.
#'
#' @inheritParams aging_msm
#'
#' @details
#' In HIV/STI modeling, there is typically an assumption that changes in
#' population size do not affect one's number of partners, specified as the
#' mean degree for network models. A person would not have 10 times the number
#' of partners should he move from a city 10 times as large. This module uses
#' the adjustment of Krivitsky et al. to adjust the edges coefficients on the
#' three network models to account for varying population size in order to
#' preserve that mean degree.
#'
#' @return
#' The network model parameters stored in \code{dat$nwparam} are updated for
#' each of the three network models.
#'
#' @references
#' Krivitsky PN, Handcock MS, and Morris M. "Adjusting for network size and
#' composition effects in exponential-family random graph models." Statistical
#' Methodology. 2011; 8.4: 319-339.
#'
#' @keywords module msm
#'
#' @export
#'
edges_correct_msm <- function(dat, at) {

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  coef.form.m <- get_nwparam(dat, network = 1)$coef.form
  coef.form.m[1] <- coef.form.m[1] + adjust
  dat$nwparam[[1]]$coef.form <- coef.form.m

  coef.form.p <- get_nwparam(dat, network = 2)$coef.form
  coef.form.p[1] <- coef.form.p[1] + adjust
  dat$nwparam[[2]]$coef.form <- coef.form.p

  coef.form.i <- get_nwparam(dat, network = 3)$coef.form
  coef.form.i[1] <- coef.form.i[1] + adjust
  dat$nwparam[[3]]$coef.form <- coef.form.i

  return(dat)
}




# HET -----------------------------------------------------------------


#' @title Network Resimulation Module
#'
#' @description Module function to resimulate the dynamic network forward one
#'              time step conditional on current network structure and vertex
#'              attributes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
simnet_het <- function(dat, at) {

  # Update edges coefficients
  dat <- edges_correct_het(dat, at)

  # Update internal ergm data
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  # Pull network parameters
  nwparam <- get_nwparam(dat, network = 1)

  # Simulate edgelist
  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam$coef.form,
                                             coef.diss = nwparam$coef.diss$coef.adj)

  return(dat)
}



#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficients in a dynamic network model
#'              to preserve the mean degree.
#'
#' @inheritParams aging_het
#'
#' @details
#' In HIV/STI modeling, there is typically an assumption that changes in
#' population size do not affect one's number of partners, specified as the
#' mean degree for network models. A person would not have 10 times the number
#' of partners should he move from a city 10 times as large. This module uses
#' the adjustment of Krivitsky et al. to adjust the edges coefficients on the
#' three network models to account for varying population size in order to
#' preserve that mean degree.
#'
#' @return
#' The network model parameters stored in \code{dat$nwparam} are updated.
#'
#' @references
#' Krivitsky PN, Handcock MS, and Morris M. "Adjusting for network size and
#' composition effects in exponential-family random graph models." Statistical
#' Methodology. 2011; 8.4: 319-339.
#'
#' @keywords module het
#'
#' @export
#'
edges_correct_het <- function(dat, at) {

  # Popsize
  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)

  # New Coefs
  coef.form <- get_nwparam(dat)$coef.form
  coef.form[1] <- coef.form[1] + log(old.num) - log(new.num)
  dat$nwparam[[1]]$coef.form <- coef.form

  return(dat)
}
