
#' @title Disclosure Module
#'
#' @description Module function for disclosure of HIV status to partners given
#'              non-disclosure in the past.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Persons who are infected may disclose their status to partners at three
#' distinct time points: at relationship onset for newly formed discordant
#' pairs; at diagnosis for pairs starting as both negative but with one newly
#' infected; or post diagnosis for one recently infected. The rates of disclosure
#' vary at these three points, and also by the partnership type.
#'
#' @return
#' This function returns the \code{dat} object with the updated disclosure list,
#' on \code{temp$discl.list}.
#'
#' @keywords module msm
#' @export
#'
disclose_msm <- function(dat, at){


  for (type in c("main", "pers", "inst")) {

    # Variables --------------------------------------------------------------

    # Attributes
    status <- dat$attr$status
    uid <- dat$attr$uid                    #uid off a little -- goes up to 10011? that's because of 10000+11 births=10011
    diag.status <- dat$attr$diag.status    #this is going to change each time you pull from births module (poisson draw)
    diag.time <- dat$attr$diag.time
    race <- dat$attr$race
    agecat2 <- dat$attr$agecat2            #agecat2 missing for new births -- is that okay? race gets assigned to these, where?
                                           #nope, not okay, needed to go back to births module and set agecat2 attribute to dat
                                           #for new births;
    # Parameters and network
    if (type == "main") {
      disc.outset.Y.prob    <- dat$param$disc.outset.main.Y.prob
      disc.at.diag.Y.prob   <- dat$param$disc.at.diag.main.Y.prob
      disc.post.diag.Y.prob <- dat$param$disc.post.diag.main.Y.prob
      disc.outset.O.prob    <- dat$param$disc.outset.main.O.prob
      disc.at.diag.O.prob   <- dat$param$disc.at.diag.main.O.prob
      disc.post.diag.O.prob <- dat$param$disc.post.diag.main.O.prob
      # disc.outset.B.prob <- dat$param$disc.outset.main.B.prob
      # disc.at.diag.B.prob <- dat$param$disc.at.diag.main.B.prob
      # disc.post.diag.B.prob <- dat$param$disc.post.diag.main.B.prob
      # disc.outset.W.prob <- dat$param$disc.outset.main.W.prob
      # disc.at.diag.W.prob <- dat$param$disc.at.diag.main.W.prob
      # disc.post.diag.W.prob <- dat$param$disc.post.diag.main.W.prob
      el <- dat$el[[1]]
    }

    if (type == "pers") {
      disc.outset.Y.prob    <- dat$param$disc.outset.pers.Y.prob
      disc.at.diag.Y.prob   <- dat$param$disc.at.diag.pers.Y.prob
      disc.post.diag.Y.prob <- dat$param$disc.post.diag.pers.Y.prob
      disc.outset.O.prob    <- dat$param$disc.outset.pers.O.prob
      disc.at.diag.O.prob   <- dat$param$disc.at.diag.pers.O.prob
      disc.post.diag.O.prob <- dat$param$disc.post.diag.pers.O.prob
      # disc.outset.B.prob <- dat$param$disc.outset.pers.B.prob
      # disc.at.diag.B.prob <- dat$param$disc.at.diag.pers.B.prob
      # disc.post.diag.B.prob <- dat$param$disc.post.diag.pers.B.prob
      # disc.outset.W.prob <- dat$param$disc.outset.pers.W.prob
      # disc.at.diag.W.prob <- dat$param$disc.at.diag.pers.W.prob
      # disc.post.diag.W.prob <- dat$param$disc.post.diag.pers.W.prob
      el <- dat$el[[2]]
    }

    if (type == "inst") {
      disc.inst.Y.prob <- dat$param$disc.inst.Y.prob
      disc.inst.O.prob <- dat$param$disc.inst.O.prob
      # disc.inst.B.prob <- dat$param$disc.inst.B.prob
      # disc.inst.W.prob <- dat$param$disc.inst.W.prob
      el <- dat$el[[3]]
    }


    # Processes --------------------------------------------------------------

    # Check for discordant rels
    posneg <- el[which(status[el[, 1]] - status[el[, 2]] == 1), , drop = FALSE]
    negpos <- el[which(status[el[, 2]] - status[el[, 1]] == 1), , drop = FALSE]
    disc.el <- rbind(posneg, negpos[, 2:1])

    # Check for not already disclosed
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    discord.cdl <- uid[disc.el[, 1]] * 1e7 + uid[disc.el[, 2]]
    notdiscl <- !(discord.cdl %in% disclose.cdl)


    # data frame of non-disclosed pairs
    nd <- disc.el[notdiscl, , drop = FALSE]

    # Check for positive diagnosis
    notdiscl.dx <- which(diag.status[nd[, 1]] == 1)

    # data frame of non-disclosed pairs where infected is dx'ed
    nd.dx <- nd[notdiscl.dx, , drop = FALSE]

    # If there are any eligible pairs
    if (nrow(nd.dx) > 0) {

      # # Split by race of pos node
      # pos.race <- race[nd.dx[, 1]]

      # Split by age of pos node
      pos.agecat2 <- agecat2[nd.dx[, 1]]

      if (type %in% c("main", "pers")) {

        # Check that rel is new
        # new.edges matrix is expressed in uid, so need to transform nd.dx
        new.edges <- dat$temp$new.edges
        new.rel <- ((uid[nd.dx[, 1]] * 1e7 + uid[nd.dx[, 2]]) %in%
                      (new.edges[, 1] * 1e7 + new.edges[, 2])) |
          ((uid[nd.dx[, 2]] * 1e7 + uid[nd.dx[, 1]]) %in%
             (new.edges[, 1] * 1e7 + new.edges[, 2]))

        # Check if diag is new
        new.dx <- diag.time[nd.dx[, 1]] == at

        #   # Assign disclosure probs
        #   dl.prob <- vector("numeric", length = nrow(nd.dx))
        #   dl.prob[pos.race == "B" & new.rel == TRUE] <- disc.outset.B.prob
        #   dl.prob[pos.race == "B" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.B.prob
        #   dl.prob[pos.race == "B" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.B.prob
        #
        #   dl.prob[pos.race == "W" & new.rel == TRUE] <- disc.outset.W.prob
        #   dl.prob[pos.race == "W" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.W.prob
        #   dl.prob[pos.race == "W" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.W.prob
        # }
        #
        # if (type == "inst") {
        #   dl.prob <- vector("numeric", length = nrow(nd.dx))
        #   dl.prob[pos.race == "B"] <- disc.inst.B.prob
        #   dl.prob[pos.race == "W"] <- disc.inst.W.prob
        # }

        # Assign disclosure probs
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        dl.prob[pos.agecat2 == "Y" & new.rel == TRUE] <- disc.outset.Y.prob
        dl.prob[pos.agecat2 == "Y" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.Y.prob
        dl.prob[pos.agecat2 == "Y" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.Y.prob

        dl.prob[pos.agecat2 == "O" & new.rel == TRUE] <- disc.outset.O.prob
        dl.prob[pos.agecat2 == "O" & new.rel == FALSE & new.dx == TRUE] <- disc.at.diag.O.prob
        dl.prob[pos.agecat2 == "O" & new.rel == FALSE & new.dx == FALSE] <- disc.post.diag.O.prob
      }

      if (type == "inst") {
        dl.prob <- vector("numeric", length = nrow(nd.dx))
        dl.prob[pos.agecat2 == "Y"] <- disc.inst.Y.prob
        dl.prob[pos.agecat2 == "O"] <- disc.inst.O.prob
      }

      # Determine disclosers
      discl <- which(rbinom(length(dl.prob), 1, dl.prob) == 1)

      # Write output
      if (length(discl) > 0) {
        discl.mat <- cbind(pos = uid[nd.dx[discl, 1]],
                           neg = uid[nd.dx[discl, 2]],
                           discl.time = at)
        dat$temp$discl.list <- rbind(dat$temp$discl.list, discl.mat)
      }
    }
  }

  if (at > 2) {
    discl.list <- dat$temp$discl.list
    master.el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])
    m <- which(match(discl.list[, 1] * 1e7 + discl.list[, 2],
                     uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]) |
                 match(discl.list[, 2] * 1e7 + discl.list[, 1],
                       uid[master.el[, 1]] * 1e7 + uid[master.el[, 2]]))
    dat$temp$discl.list <- discl.list[m, ]
  }

  return(dat)
}
