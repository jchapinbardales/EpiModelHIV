# MSM -----------------------------------------------------------------

#' @title Calculate Target Statistics for Network Model Estimation
#'
#' @description Calculates the target statistics for the formation and dissolution
#'              components of the network model to be estimated with \code{netest}.
#'
#' @param time.unit Time unit relative to 1 for daily.
#' @param method Method for calculating target statistics by race, with options of
#'        \code{2} for preserving race-specific statistics and \code{1} for
#'        averaging over the statistics and dropping the race-specific terms.
#' @param num.B Population size of black MSM.
#' @param num.W Population size of white MSM.
#' @param num.Y Population size of Young MSM (15-24).
#' @param num.O Population size of Older MSM (25-40).
#' #@param deg.mp.B Degree distribution matrix for main and casual partners for
#' #       black MSM, as a 2 by 3 matrix.
#' #@param deg.mp.W Degree distribution matrix for main and causal partners for
#' #       white MSM, as a 2 by 3 matrix.
#' @param deg.mp.Y Degree distribution matrix for main and casual partners for
#'        young MSM, as a 2 by 3 matrix.
#' @param deg.mp.O Degree distribution matrix for main and causal partners for
#'        older MSM, as a 2 by 3 matrix.
#' #@param mdeg.inst.B Mean degree, or rate, of one-off partnerships per day
#' #       for black MSM.
#' #@param mdeg.inst.W Mean degree, or rate, of one-off partnerships per day
#' #       for white MSM.
#' @param mdeg.inst.Y Mean degree, or rate, of one-off partnerships per day
#'        for young MSM.
#' @param mdeg.inst.O Mean degree, or rate, of one-off partnerships per day
#'        for older MSM.
#' @param qnts.B Means of one-off rates split into quintiles for white MSM. Use
#'        \code{NA} to ignore these quantiles in the target statistics.
#' @param qnts.W Means of one-off rates split into quintiles for black MSM. Use
#'        \code{NA} to ignore these quantiles in the target statistics.
#' @param prop.hom.mpi.B A vector of length 3 for the proportion of main, casual,
#'        and one-off partnerships in same race for black MSM.
#' @param prop.hom.mpi.W A vector of length 3 for the proportion of main, casual,
#'        and one-off partnerships in same race for white MSM.
#' @param balance Method for balancing of edges by race for number of mixed-race
#'        partnerships, with options of \code{"black"} to apply black MSM counts,
#'        \code{"white"} to apply white MSM counts, and \code{"mean"} to take
#'        the average of the two expectations.
#' @param sqrt.adiff.BB Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and one-off black-black
#'        partnerships.
#' @param sqrt.adiff.WW Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and one-off white-white
#'        partnerships.
#' @param sqrt.adiff.BW Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and one-off black-white
#'        partnerships.
#' @param diss.main Dissolution model formula for main partnerships.
#' @param diss.pers Dissolution model formula for casual partnerships.
#' @param durs.main Vector of length 3 with the duration of BB, BW, and WW main
#'        partnerships in days.
#' @param durs.pers Vector of length 3 with the duration of BB, BW, and WW
#'        casual partnerships in days.
#' @param durs.main.age Vector of length 3 with the duration of YY, YO, and OO main
#'        partnerships in days.
#' @param durs.pers.age Vector of length 3 with the duration of YY, YO, and OO
#'        casual partnerships in days.
#' @param ages Integer vector of ages in years that defines range of possible
#'        initial ages in the population.
#' @param asmr.B Vector of length 40 defining the age-specific
#'        mortality rate for persons within that age slot, for black MSM.
#' @param asmr.W Vector of length 40 defining the age-specific
#'        mortality rate for persons within that age slot, for white MSM.
#' #@param role.B.prob Vector of length 3 for the probability of sexual role as
#' #       insertive, receptive, and versatile, for black MSM.
#' #@param role.W.prob Vector of length 3 for the probability of sexual role as
#' #       insertive, receptive, and versatile, for white MSM.
#' @param role.Y.prob Vector of length 3 for the probability of sexual role as
#'        insertive, receptive, and versatile, for young MSM.
#' @param role.O.prob Vector of length 3 for the probability of sexual role as
#'        insertive, receptive, and versatile, for older MSM.
#'
#' @details
#' This function performs basic calculations to determine the components of the
#' formationa and dissolution models for the network model estimation to be
#' conducted with \code{\link{netest}}. The inputs inputs for this function are
#' calculated externally to the package in a setup scenario file.
#'
#' @keywords msm
#'
#' @seealso
#' Network statistics calculated here are entered into \code{\link{base_nw_msm}}
#' to construct the base network, and then into the parameters in
#' \code{\link{param_msm}}.
#'
#' @export
#'
calc_nwstats_msm <- function(time.unit = 7,
                             #method = 2,
                             num.B, 
                             num.W, 
                             num.Y, 
                             num.O, 
                             deg.mp.Y,
                             deg.mp.O,
                             deg.mp.overall,
                             mdeg.inst.Y,
                             mdeg.inst.O,
                             mdeg.inst.overall,
                             qnts.B, 
                             qnts.W, 
                             prop.hom.mpi.B, 
                             prop.hom.mpi.W, 
                             balance = "mean", 
                             sqrt.adiff.BB,
                             sqrt.adiff.WW,
                             sqrt.adiff.BW, 
                             diss.main, 
                             diss.pers, 
                             durs.main, 
                             durs.pers,
                             #durs.main.age, 
                             #durs.pers.age,
                             ages, 
                             asmr.B, 
                             asmr.W, 
                             role.Y.prob,
                             role.O.prob,
                             role.prob.overall) {
  
  if (sum(deg.mp.Y) != 1) {
    stop("deg.mp.Y must sum to 1.")
  }
  if (sum(deg.mp.O) != 1) {
    stop("deg.mp.O must sum to 1.")
  }
  
  
  
  #NEW PT & AGE VARS
  #num.Y = num.Y, 
  #num.O = num.O,
  #deg.mp.Y = deg.mp.Y,
  #deg.mp.O = deg.mp.O,
  #mdeg.inst.Y = mdeg.inst.Y,
  #mdeg.inst.O = mdeg.inst.O,
  #durs.main.age = durs.main.age,
  #durs.pers.age = durs.pers.age,
  #role.Y.prob = role.Y.prob,
  #role.O.prob = role.O.prob,
  
  #browser()
  
  #num <- num.B + num.W
  #10000
  num <- num.Y + num.O
  
  # deg.pers nodal attribute
  #if (method == 2) {
  #  deg.pers.B <- apportion_lr(num.B, c("B0", "B1", "B2"), colSums(deg.mp.B))
  #  deg.pers.W <- apportion_lr(num.W, c("W0", "W1", "W2"), colSums(deg.mp.W))
  #}
  
  #if (method == 1) {
  #  deg.pers <- apportion_lr(num, 0:2, colSums(deg.mp.W)) #same deg matrix for B & W
  #}
  
  #gives all 10000 people an attribute of casual 0,1,or 2;
  deg.pers.overall <- apportion_lr(num, 0:2, colSums(deg.mp.overall)) #mean=0.3986
  
  #gives all 10000 people an attribute of age+casual 0,1,or 2;
  deg.pers.Y <- apportion_lr(num.Y, c("Y0", "Y1", "Y2"), colSums(deg.mp.Y))  #colsum = totals of 0C, 1C, 2C for Young;
  deg.pers.O <- apportion_lr(num.O, c("O0", "O1", "O2"), colSums(deg.mp.O))  #colsum = totals of 0C, 1C, 2C for Old;
  #sum over of main;
  
  
  
  # deg main nodal attribute
  #if (method == 2) {
  #  deg.main.B <- apportion_lr(num.B, c("B0", "B1"), rowSums(deg.mp.B))
  #  deg.main.W <- apportion_lr(num.W, c("W0", "W1"), rowSums(deg.mp.W))
  #}
  #if (method == 1) {
  #  deg.main <- apportion_lr(num, 0:1, rowSums(deg.mp.W))
  #}
  
  #overall deg.main for fitting model;
  deg.main.overall <- apportion_lr(num, 0:1, rowSums(deg.mp.overall)) #mean=0.2882
  
  #gives all 10000 people an attribute of age+main 0, or 1;
  deg.main.Y <- apportion_lr(num.Y, c("Y0", "Y1"), rowSums(deg.mp.Y))  #rowsum = totals of 0M, 1M for Young;
  deg.main.O <- apportion_lr(num.O, c("O0", "O1"), rowSums(deg.mp.O))  #rowsum = totals of 0M, 1M for Old;
  #sum over levels of casual;
  
  
  # Main partnerships -------------------------------------------------------
  
  # Persons in main partnerships by casual degree
  # if (method == 2) {
  #   totdeg.m.by.dp <- c(num.B * deg.mp.B[2, ], num.W * deg.mp.W[2, ])
  # }
  # if (method == 1) {
  #   totdeg.m.by.dp <- c(num * deg.mp.B[2, ]) ##edges=mean deg*n/2, so here n*deg main,0C ; n*mdeg,1C ; n*mdeg,2C...summed and /2 below;
  # }
  
  #overall, for fitting;
  totdeg.m.by.dp.overall <- c(num * deg.mp.overall[2, ])
  #num * % of people with main partner by deg casual
  #number of people with a main partner among those who have 0,1,2 casual partners
  #2183.982  487.842  210.000
  
  totdeg.m.by.dp <- c(num.Y * deg.mp.Y[2, ], num.O * deg.mp.O[2, ])  #num.Y * % of people with main partner (across row2)
  #745.986  166.185   77.553 1437.996  321.657  132.447
  #Y,0C     Y,1C      Y,2C   O,0C      O,1C     O,2C
  #number of people with a main partner among those who are Y/O with X# casual partners
  #makes sense,  if have no casual partners, probs have main edge
  #if 1C,2C, then going to have less main edges 
  
  # Persons in partnerships by race
  # if (method == 2) {
  #totdeg.m.by.race <- c(sum(totdeg.m.by.dp[1:3]), sum(totdeg.m.by.dp[4:6]))
  # }
  
  # dont need unless want homophily more by age beyond sqrtabsdiff
  #totdeg.m.by.age <- c(sum(totdeg.m.by.dp[1:3]), sum(totdeg.m.by.dp[4:6]))
  #989.724                    1892.100
  
  
  # Number of main partnerships -- weighted average across age (bc dist of age not 50/50)
  #overallmeandeg=dY*num.Y/10000+dO*num.O/10000 = # of Ypp with partners + # Opp with partners(an edge)/10000
  #dY*num.Y+dO*num.O = # of Ypp with partners + # Opp with partners (an edge)
  #BUT also want this by casual degree
  #so this is what is doing in totdeg.m.by.dp
  
  #edges.m <- (sum(totdeg.m.by.dp)) / 2   #divide by 2 because from pop to edges, 2 pp per edge
  
  edges.m <- (sum(totdeg.m.by.dp)) / 2   #sum of number of main partners across ages and those with 0/1/2 casuals
  #1440.912
  
  edges.m.overall <- (sum(totdeg.m.by.dp.overall)) / 2  #sum of number of main partners across those with 0/1/2 casuals
  #1440.912--makes sense add up number of main across levels of casuals by age=same sum overall;
  
  
  # Mixing
  # if (method == 2) {
  #   # Number of mixed-race partnerships, with balancing to decide
  #   edges.m.B2W <- totdeg.m.by.race[1] * (1 - prop.hom.mpi.B[1])
  #   edges.m.W2B <- totdeg.m.by.race[2] * (1 - prop.hom.mpi.W[1])
  #  edges.het.m <- switch(balance,
  #                         black = edges.m.B2W,
  #                          white = edges.m.W2B,
  #                         mean = (edges.m.B2W + edges.m.W2B) / 2)
  # 
  #   # Number of same-race partnerships
  #   edges.hom.m <- (totdeg.m.by.race - edges.het.m) / 2
  # 
  #   # Nodemix target stat: number of BB, BW, WW partnerships
  #   edges.nodemix.m <- c(edges.hom.m[1], edges.het.m, edges.hom.m[2])
  # }
  
  # Sqrt absdiff term for age
  # if (method == 2) {
  #   sqrt.adiff.m <- edges.nodemix.m * c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1])
  # }
  # if (method == 1) {
  #sqrt.adiff.m <- edges.m * mean(c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1]))
  # } #668.1029
  
  #Number of YOUNG men with main partner
  
  sqrt.adiff.m <- edges.m * mean(c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1]))
  #668.1029
  
  sqrt.adiff.m.overall <- edges.m.overall * mean(c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1]))
  #668.1029
  
  totdeg.m.Y <- c(num.Y * sum(deg.mp.Y[2, ]))
  #989.724
  
  # Compile target stats
  # if (method == 2) {
  #  stats.m <- c(edges.m, edges.nodemix.m[2:3], totdeg.m.by.dp[c(2:3, 5:6)], sqrt.adiff.m)
  # }
  #if (method == 1) {
  stats.m.overall <- c(edges.m.overall, totdeg.m.by.dp.overall[2:3], totdeg.m.Y, sqrt.adiff.m.overall)
  # 1440.9120  487.8420         210.0000                             #989.724    668.1029
  # edges,    #Mp people w/ 1C  2C                                   #Ypp w/ Mp  sqrtage (no change)
  
  # }
  # if (method == 1) {
  #stats.m <- c(edges.m, totdeg.m.by.dp[2:3], sqrt.adiff.m)
  # }
  
  
  stats.m <- c(edges.m, totdeg.m.by.dp[c(2:3, 5:6)], sqrt.adiff.m)
  #1440.9120  166.1850        77.5530   321.6570           132.4470  668.1029
  #totedges   #Mp among Y,1C and Y,2C   #Mp among O,1C and O,2C      sqrtabsdiff in M
  
  #no age homophily beyond the sq diffs so dont need edges.nodemix.m[2:3]
  # totdeg.m.by.dp[c(2:3, 5:6)] = totdeg.m.by.dp, # pp with Mp's for Y1C,Y2C,O1C,O2C - we have deg.main.Y/O so 
  # don't need to specify # w/ M for Y0C or O0C, have degree of freedom for each
  # 1440.91200   98.12611  896.98695  166.18500   77.55300  321.65700  132.44700  668.10286
  #   Error in vector.namesmatch(target.stats, names(nw.stats)) : 
  #     Length of "target.stats" is 8 but should be 7.
  
  # Dissolution model
  exp.mort <- (mean(asmr.B[ages]) + mean(asmr.W[ages])) / 2 #3.588345e-05
  
  coef.diss.m <- dissolution_coefs(dissolution = diss.main,
                                   duration = durs.main / time.unit, 
                                   d.rate = exp.mort)
  # Dissolution Coefficients
  # =======================
  #   Dissolution Model: ~offset(edges)
  # Target Statistics: 23.65325
  # Crude Coefficient: 3.120303
  # Mortality/Exit Rate: 3.588345e-05
  # Adjusted Coefficient: 3.122002
  
  
  # Casual partnerships -----------------------------------------------------
  
  # Persons in partnerships by main degree
  # if (method == 2) {
  #   totdeg.p.by.dm <- c(num.B * deg.mp.B[, 2] + num.B * deg.mp.B[, 3] * 2,
  #                       num.W * deg.mp.W[, 2] + num.W * deg.mp.W[, 3] * 2)
  # }
  # if (method == 1) {
  #   totdeg.p.by.dm <- c(num * deg.mp.B[, 2] + num * deg.mp.B[, 3] * 2)
  # }
  totdeg.p.by.dm.overall <- c(num * deg.mp.overall[, 2] + num * deg.mp.overall[, 3] * 2)
  #3079.043  907.842
  #0M        1M - number of casual partners among main 0,1
  
  totdeg.p.by.Y<- c(num.Y * sum(deg.mp.Y[, 2]) + num.Y * sum(deg.mp.Y[, 3]) * 2)
  #1344.252, rather than summing here though, could have left as for overall,
  #gotten 2 values (0MY, 1MY), and then just taken totdeg.p.by.Y[2] in stats.p.overall
  
  
  totdeg.p.by.dm <- c(num.Y * deg.mp.Y[, 2] + num.Y * deg.mp.Y[, 3] * 2,
                      num.O * deg.mp.O[, 2] + num.O * deg.mp.O[, 3] * 2)
  #1022.961  321.291 2056.082  586.551
  #Y,0M      Y,1M    O,0M      O,1M
  #number of casual partners among age+main 0,1 -- need to multiple 3rd col by 2 bc 2C cats;
  #again, makes sense...if have no main partner, a lot of casual edges
  #greater number of casual edges in Older regardless of main status
  
  
  # Persons in partnerships by race
  # if (method == 2) {
  #   totdeg.p.by.race <- c(sum(totdeg.p.by.dm[1:2]), sum(totdeg.p.by.dm[3:4]))
  # }
  
  totdeg.p.by.age <- c(sum(totdeg.p.by.dm[1:2]), sum(totdeg.p.by.dm[3:4]))
  #1344.252 2642.633
  
  # Persons concurrent -- concurrent casuals! -- what about 1C, 1M concurrency?;
  # if (method == 2) {
  #   conc.p.by.race <- c(sum(deg.mp.B[, 3]) * num.B, sum(deg.mp.W[, 3]) * num.W)
  # }
  # if (method == 1) {
  #   conc.p <- sum(deg.mp.B[, 3] * num)
  # }
  
  conc.p.overall <- sum(deg.mp.overall[, 3] * num)
  #920.289 - # of pp(regardless of main) who had 2 concurrent Casuals;
  
  conc.p.by.age <- c(sum(deg.mp.Y[, 3]) * num.Y, sum(deg.mp.O[, 3]) * num.O)
  #sum deg of Y-0M,1M among those with 2C
  # # of Ypp(regardless of main) who had 2 concurrent Casuals, # of Opp who had 2 concurrent Casuals
  # conc.p.by.age = 276.975 643.314  --> more older people with 2Casuals
  
  # Number of partnerships
  edges.p <- sum(totdeg.p.by.dm) / 2
  #1993.442
  
  edges.p.overall <- sum(totdeg.p.by.dm.overall) / 2
  #1993.442 -- makes sense, adding up across ages should be same as overall
  
  # Mixing
  # if (method == 2) {
  #   # Number of mixed-race partnerships, with balancing to decide
  #   edges.p.B2W <- totdeg.p.by.race[1] * (1 - prop.hom.mpi.B[2])
  #   edges.p.W2B <- totdeg.p.by.race[2] * (1 - prop.hom.mpi.W[2])
  #   edges.het.p <- switch(balance,
  #                         black = edges.p.B2W, white = edges.p.W2B,
  #                         mean = (edges.p.B2W + edges.p.W2B) / 2)
  # 
  #   # Number of same-race partnerships
  #   edges.hom.p <- (totdeg.p.by.race - edges.het.p) / 2
  # 
  #   # Nodemix target stat: number of BB, BW, WW partnerships
  #   edges.nodemix.p <- c(edges.hom.p[1], edges.het.p, edges.hom.p[2])
  # }
  
  # Sqrt absdiff term for age
  # if (method == 2) {
  #   sqrt.adiff.p <- edges.nodemix.p * c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2])
  # }
  # if (method == 1) {
  #   sqrt.adiff.p <- edges.p * mean(c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2]))
  # }
  
  sqrt.adiff.p <- edges.p * mean(c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2]))
  #1168.822
  
  sqrt.adiff.p.overall <- edges.p.overall * mean(c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2]))
  #1168.822
  
  # Compile target statistics
  # if (method == 2) {
  #   stats.p <- c(edges.p, edges.nodemix.p[2:3], totdeg.p.by.dm[c(2, 4)],
  #                conc.p.by.race, sqrt.adiff.p)
  # }
  # if (method == 1) {
  #   stats.p <- c(edges.p, totdeg.p.by.dm[2], conc.p, sqrt.adiff.p)
  # }
  
  #stats.p.overall <- c(edges.p.overall, totdeg.p.by.dm.overall[2], totdeg.p.by.Y, conc.p.overall, sqrt.adiff.p.overall)
  #1993.442  907.842            1344.252                              920.289 1168.822
  #edges     #Cpartners of 1M   #Cpartners among Y (summed over Mp)   conc    sqrtage
  stats.p.overall <- c(edges.p.overall, totdeg.p.by.dm.overall[2], totdeg.p.by.Y, conc.p.by.age[1:2], sqrt.adiff.p.overall)
  #1993.442   907.842           1344.252            276.975        643.314          1168.822
  #edges     #Cpartners of 1M   #Cpartners among Y  #Y w/ 2ConCasp #O w/ 2ConCasp   sqrtage
                                                    #need to specify both ts when term w/ "by" in it
  
  stats.p <- c(edges.p, totdeg.p.by.dm[c(2, 4)],  #what about totdeg.p.by.age? is this only used if going to account for mixing beyond sq root?
               conc.p.by.age, sqrt.adiff.p)       #i think don't need unless mixing beyond sqrtabsdiff
  #totdeg.p.by.dm[c(2, 4)] = # pp with Cas partners for Y1M, O1M
  # we have deg.pers.Y/O so don't need to specify # w/ Cas 
  # for Y0M or O0M, have degree of freedom for eac
  #1993.442  321.291  586.551  276.975        643.314 1168.822
  #edges     Cdeg,Y1M Cdeg,O1M #Y2concurrCas  #O2Cs   sqrt.adiff.p
  
  # Dissolution model
  coef.diss.p <- dissolution_coefs(dissolution = diss.pers,
                                   duration = durs.pers / time.unit, #durs.pers, when race;
                                   d.rate = exp.mort)
  
  # Dissolution Coefficients
  # =======================
  #   Dissolution Model: ~offset(edges)
  # Target Statistics: 14.20763
  # Crude Coefficient: 2.580795
  # Mortality/Exit Rate: 3.588345e-05
  # Adjusted Coefficient: 2.581815
  
  
  # Instant partnerships ----------------------------------------------------
  
  # Number of instant partnerships per time step, by main and casl degree
  # if (method == 2) {
  #   num.inst.B <- num.B * deg.mp.B * mdeg.inst.B * time.unit
  #   num.inst.W <- num.W * deg.mp.W * mdeg.inst.W * time.unit
  # }
  # if (method == 1) {
  #   num.inst <- num * deg.mp.W * mdeg.inst.W * time.unit
  # }
  
  num.inst.overall <- num * deg.mp.overall * mdeg.inst * time.unit
  #320.1473 146.56291 60.96462
  #127.1913  28.76227 12.38121
  #one off rates by PT combo overall
  
  num.inst.Y.overall <- num.Y * mdeg.inst.Y.overall * time.unit
  #3693 * 0.009749584 * 7
  #one off rates for Young overall (don't need to number of people across
  #the M/C deg matrix since want collapsed over PT)
  
  num.inst.Y <- num.Y * deg.mp.Y * mdeg.inst.Y * time.unit  
  num.inst.O <- num.O * deg.mp.O * mdeg.inst.O * time.unit
  #num.Y * deg.mp.Y = # Ypp in M/C matrix 
  #                     * daily probability of having a one-off partnership given Y and M/C status
  #                     * 7 to make it per week
  # one-off partners across M/C matrix by age (two 2x3 matrices) - one off rates by PT combo and by age
  #Young
  #112.09436 45.34834 18.018974
  #43.30001 13.80831  6.443879
  #Old
  #206.50253 100.92444 42.494428
  #83.90994  14.41924  5.937334
  
  # Risk quantiles
  # if (!is.na(qnts.B[1]) & !is.na(qnts.W[1])) {
  #   if (method == 2) {
  #     num.riskg.B <- (0.2*num.B) * qnts.B * time.unit
  #     num.riskg.W <- (0.2*num.W) * qnts.W * time.unit
  #   }
  #   if (method == 1) {
  #num.riskg <- 0.2 * num * qnts.B * time.unit
  #   }
  # }
  
  num.riskg <- 0.2 * num * qnts.B * time.unit  #really B and W same
  #20% (ie. quintile) * 10000 * avg. rate of AI based on quintiles * 7
  #2000 people * 0 * 7
  #2000 people * 0.00100413270621469 AI acts per person per day * 7 days/wk = 14.05786 AI acts per week 
  # are these one-off AI act rates?? or overall? says overall...
  # num.riskg = 0.00000  14.05786  75.13392 142.73906 441.29450
  
  # Number of instant partnerships per time step, by race
  # if (method == 2) {
  #   totdeg.i <- c(sum(num.inst.B), sum(num.inst.W))
  # }
  # if (method == 1) {
  #   totdeg.i <- sum(num.inst)
  # }
  
  totdeg.i <- c(sum(num.inst.Y), sum(num.inst.O))
  #239.0139          454.1879  -- summed # one-off partners across M/C matrix
  
  totdeg.i.overall <- sum(num.inst.overall)
  #696.0096
  
  # Number of partnerships
  edges.i <- sum(totdeg.i) / 2
  #346.6009
  
  edges.i.overall <- sum(totdeg.i.overall) / 2
  #348.0048 -- hm, why are edges overall and by age for M and C but not I?
  
  # Mixing
  # if (method == 2) {
  #   # Number of mixed-race partnerships, with balancing to decide
  #   edges.i.B2W <- totdeg.i[1] * (1 - prop.hom.mpi.B[3])
  #   edges.i.W2B <- totdeg.i[2] * (1 - prop.hom.mpi.W[3])
  #   edges.het.i <- switch(balance,
  #                         black = edges.i.B2W, white = edges.i.W2B,
  #                         mean = (edges.i.B2W + edges.i.W2B) / 2)
  # 
  #   # Number of same-race partnerships
  #   edges.hom.i <- edges.i - edges.het.i
  # 
  #   # Nodemix target stat: number of BB, BW, WW partnerships
  #   edges.nodemix.i <- c((totdeg.i[1] - edges.het.i) / 2,
  #                        edges.het.i,
  #                        (totdeg.i[1] - edges.het.i) / 2)
  # }
  # 
  # if (method == 2) {
  #   sqrt.adiff.i <- edges.nodemix.i * c(sqrt.adiff.BB[3], sqrt.adiff.BW[3], sqrt.adiff.WW[3])
  # }
  # if (method == 1) {
  #   sqrt.adiff.i <- edges.i * mean(c(sqrt.adiff.BB[3], sqrt.adiff.BW[3], sqrt.adiff.WW[3]))
  # }
  
  sqrt.adiff.i <- edges.i * mean(c(sqrt.adiff.BB[3], sqrt.adiff.BW[3], sqrt.adiff.WW[3]))
  #188.4354
  
  sqrt.adiff.i.overall <- edges.i.overall * mean(c(sqrt.adiff.BB[3], sqrt.adiff.BW[3], sqrt.adiff.WW[3]))
  #189.1986
  
  # if (!is.na(qnts.B[1]) & !is.na(qnts.W[1])) {
  #   if (method == 2) {
  #     stats.i <- c(edges.i, num.inst.B[-1], num.inst.W,
  #                  num.riskg.B[-3], num.riskg.W[-3],
  #                  edges.hom.i, sqrt.adiff.i)
  #   }
  #   if (method == 1) {
  #     stats.i <- c(edges.i, num.inst[-1], num.riskg[-3], sqrt.adiff.i)
  #   }
  # 
  # } else {
  #   if (method == 2) {
  #     stats.i <- c(edges.i, num.inst.B[-1], num.inst.W, edges.hom.i, sqrt.adiff.i)
  #   }
  #   if (method == 1) {
  #     stats.i <- c(edges.i, num.inst[-1], sqrt.adiff.i)
  #   }
  # }
  
  if (!is.na(qnts.B[1]) & !is.na(qnts.W[1])) {
    
    stats.i <- c(edges.i, num.inst.Y[-1], num.inst.O,   # not num.riskg by race, just overall, no homophily;
                 num.riskg[-3], sqrt.adiff.i)           # num.inst.Y[-1] takes out first value 0M0C
    # num.riskg[-3] takes out the middle group (75.13) 
    
    stats.i.overall <- c(edges.i.overall, num.inst.overall[-1], num.inst.Y.overall, num.riskg[-3], sqrt.adiff.i.overall)
  }
  else {
    stats.i <- c(edges.i, num.inst.Y[-1], num.inst.O,   #no riskg, no homophily;
                 sqrt.adiff.i)
    
    stats.i.overall <- c(edges.i.overall, num.inst.overall[-1], num.inst.Y.overall, sqrt.adiff.i.overall)
  }
  
  # 346.600894             43.300011  45.348341  13.808312  18.018974  6.443879 
  # edges      num.inst.Y: Y1M,0C     Y0M,1C     Y1M,1C     Y0M,2C     Y1M,2C  -- doesnt include 1st one Y0M,0C
  
  #                       206.502533  83.909943 100.924437  14.419240  42.494428  5.937334
  #           num.inst.O: O 0M,0C     O1M,0C     O0M,1C     O1M,1C     O0M,2C     O1M,2C
  
  # 0.000000  14.057858 142.739065 441.294497                                 188.435353
  # num.riskg[-3], doesn't include the middle group (75.13) -- not sure why?  sqrt.adiff.i
  
  #stats.i.overall 
  #348.00482  127.19130 146.56291  28.76227  60.96462  12.38121   0.00000  14.05786 142.73906 441.29450 189.19862
  #edges      #1M,0C    #0M,1C     #1M,1C    #0M,2C    #1M,2C     #num.riskg[-3]                         sqrt.adiff.i
              ##one off partner rates by PT combo overall
  
              ##need a num.inst target stat for agecat2 overall (not by PT)
              #num.inst.Y.overall
  
  
  
  # Compile results ---------------------------------------------------------
  out <- list()
  # out$method <- method
  # if (method == 2) {
  #   out$deg.pers <- c(deg.pers.B, deg.pers.W)
  #   out$deg.main <- c(deg.main.B, deg.main.W)
  # }
  # if (method == 1) {
  #   out$deg.pers <- deg.pers
  #   out$deg.main <- deg.main
  # }
  
  #out$deg.pers <- c(deg.pers.Y, deg.pers.O)  #all nodes, attributes of age+#casual
  #out$deg.main <- c(deg.main.Y, deg.main.O)  #all nodes, attributes of age+#main
  
  out$deg.pers <- deg.pers.overall  #all nodes, attributes of #casual
  out$deg.main <- deg.main.overall  #all nodes, attributes of #main
  
  # out$stats.m <- stats.m
  # out$stats.p <- stats.p
  # out$stats.i <- stats.i
  
  out$stats.m <- stats.m.overall
  out$stats.p <- stats.p.overall
  out$stats.i <- stats.i.overall
  
  out$coef.diss.m <- coef.diss.m
  out$coef.diss.p <- coef.diss.p
  
  out$ages <- ages
  out$asmr.B <- asmr.B
  out$asmr.W <- asmr.W
  
  out$time.unit <- time.unit
  out$num.Y <- num.Y
  out$num.O <- num.O
  
  out$deg.mp.overall <- deg.mp.overall
  out$deg.mp.Y <- deg.mp.Y
  out$deg.mp.O <- deg.mp.O
  
  out$role.Y.prob <- role.Y.prob
  out$role.O.prob <- role.O.prob
  
  class(out) <- "nwstats"
  return(out)
}



#' @title Construct Base Network for Model Estimation and Simulation
#'
#' @description Initializes the base network for model estimation within
#'              \code{netest}.
#'
#' @param nwstats An object of class \code{nwstats}, as output from
#'        \code{\link{calc_nwstats_msm}}.
#'
#' @details
#' This function takes the output of \code{\link{calc_nwstats_msm}} and constructs
#' an empty network with the necessary attributes for race, square root of age,
#' and sexual role class. This base network is used for all three network
#' estimations.
#'
#' @seealso
#' The final vertex attributes on the network for cross-network degree are
#' calculated and set on the network with \code{\link{assign_degree}}.
#'
#' @keywords msm
#' @export
#'
base_nw_msm <- function(nwstats) {
  
  num.Y <- nwstats$num.Y
  num.O <- nwstats$num.O
  
  # Initialize network
  n <- num.Y + num.O
  nw <- network::network.initialize(n, directed = FALSE)
  
  # Calculate attributes
  race <- c(rep("B", num.B), rep("W", num.W))  ## how do I assign both race and age if don't specify num.B?
  race <- sample(race)  # this just randomly shuffles BW race around the 10000 people sample
  # in accordance with num.B and num.W distribution
  
  # let's convert these to percentages so that n can be arbitrary
  pct.Y <- nwstats$num.Y / (nwstats$num.Y + nwstats$num.O)
  pct.O <- 1 - pct.Y
  
  # this is the vector of potential ages to sample from
  ager <- nwstats$ages
  ages <- seq(min(ager), max(ager) + 1, 1 / (365 / nwstats$time.unit))
  
  # then sample once for the young ages, once for older
  age <- c(sample(ages[ages < 25], round(pct.Y * n), TRUE),
           sample(ages[ages >= 25], round(pct.O * n), TRUE))
  
  # then just to be sure, reshuffle the whole list
  age <- sample(age)
  
  # check it
  #mean(age < 25)
  #pct.Y
  #mean(age >= 25)
  #pct.O
  
  # finally, create derived variables
  sqrt.age <- sqrt(age)
  agecat2 <- rep(NA, length(age))
  agecat2[age < 25] <- "Y"
  agecat2[age >= 25] <- "O"
  
  # check it again
  #mean(agecat2 == "Y")
  #mean(agecat2 == "O")
  
  
  role.Y <- sample(apportion_lr(num.Y, c("I", "R", "V"), nwstats$role.Y.prob))
  role.O <- sample(apportion_lr(num.O, c("I", "R", "V"), nwstats$role.O.prob))
  role <- rep(NA, n)
  #role[race == "B"] <- role.B
  #role[race == "W"] <- role.W
  role[agecat2 == "Y"] <- role.Y
  role[agecat2 == "O"] <- role.O
  #check -- mean(role=="I" & agecat2=="Y") = 0.0516 = 0.13(%I given Y)*0.37(%Y) = joint prob
  
  # riskg.B <- sample(apportion_lr(num.Y, 1:5, rep(0.2, 5)))
  # riskg.W <- sample(apportion_lr(num.O, 1:5, rep(0.2, 5)))
  # riskg <- rep(NA, n)
  # riskg[race == "B"] <- riskg.B
  # riskg[race == "W"] <- riskg.W
  
  riskg.Y <- sample(apportion_lr(num.Y, 1:5, rep(0.2, 5)))  #may need to use qnts.B?
  riskg.O <- sample(apportion_lr(num.O, 1:5, rep(0.2, 5)))
  riskg <- rep(NA, n)
  riskg[agecat2 == "Y"] <- riskg.Y
  riskg[agecat2 == "O"] <- riskg.O
  
  attr.names <- c("race", "agecat2", "riskg", "sqrt.age", "role.class")
  attr.values <- list(race, agecat2, riskg, sqrt.age, role)
  nw <- network::set.vertex.attribute(nw, attr.names, attr.values)
  
  return(nw)
}


#' @title Assign Degree Vertex Attribute on Network Objects
#'
#' @description Assigns the degree vertex attributes on network objects
#'              conditional on their values from the other networks.
#'
#' @param nw Object of class \code{network} that is the target for the vertex
#'        attribute.
#' @param deg.type Type of degree to assign to \code{nw}, with options of
#'        \code{"pers"} to assign casual degree onto main network and
#'        \code{"main"} to assign main degree to casual network.
#' @param nwstats Object of class \code{nwstats}.
#'
#' @details
#' This function assigns the degree of other networks as a vertex attribute on the
#' target network given a bivariate degree mixing matrix of main, casual, and
#' one-partnerships contained in the \code{nwstats} data.
#'
#' @keywords msm
#' @export
#'
assign_degree <- function(nw, deg.type, nwstats) {
  
  #trying to do diagnostics on the fitted network with no degrees by PT -- may 
  # need to change this below to be deg.mp.overall to get correct assigned degrees.
  # above should all work, no degree variables in the base network.
  
  
  if (!("network" %in% class(nw))) {
    stop("nw must be of class network")
  }
  
  # if (deg.type == "main") {
  #   attr.name <- "deg.main"
  #   dist.B <- rowSums(nwstats$deg.mp.B)
  #   dist.W <- rowSums(nwstats$deg.mp.W)
  # }
  # if (deg.type == "pers") {
  #   attr.name <- "deg.pers"
  #   dist.B <- colSums(nwstats$deg.mp.B)
  #   dist.W <- colSums(nwstats$deg.mp.W)
  # }
  
  
  # if (deg.type == "main") {
  #   attr.name <- "deg.main"
  #   dist.Y <- rowSums(nwstats$deg.mp.overall) #deg.mp.Y     #0.716 0.284, % of pp with 0M, 1M
  #   dist.O <- rowSums(nwstats$deg.mp.overall) #deg.mp.O
  # }
  # if (deg.type == "pers") {
  #   attr.name <- "deg.pers"
  #   dist.Y <- colSums(nwstats$deg.mp.overall) #deg.mp.Y     #0.6970 0.2145 0.0885, % of pp with 0C, 1C, 2C
  #   dist.O <- colSums(nwstats$deg.mp.overall) #deg.mp.O
  # }
  if (deg.type == "main") {
    attr.name <- "deg.main"
    attr.name2 <- "deg.main.age"
    dist <- rowSums(nwstats$deg.mp.overall)     #0.716 0.284, % of pp with 0M, 1M
  }                        #deg.mp.Y - 0.732 0.268
  
  if (deg.type == "pers") {
    attr.name <- "deg.pers"
    attr.name2 <- "deg.pers.age"
    dist <- colSums(nwstats$deg.mp.overall)      #0.6970 0.2145 0.0885, % of pp with 0C, 1C, 2C
  }                        #deg.mp.Y - 0.711 0.214 0.075
  
  
  # if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.B)), 1, tolerance = 5e-6))) {
  #   stop("B degree distributions do not sum to 1")
  # }
  # 
  # if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.W)), 1, tolerance = 5e-6))) {
  #   stop("W degree distributions do not sum to 1")
  # }
  
  # if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.Y)), 1, tolerance = 5e-6))) {
  #   stop("Y degree distributions do not sum to 1")
  # }
  # 
  # if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.O)), 1, tolerance = 5e-6))) {
  #   stop("O degree distributions do not sum to 1")
  # }
  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.overall)), 1, tolerance = 5e-6))) {
    stop("Degree distributions do not sum to 1")
  }
  
  #browser()
  
  # race <- get.vertex.attribute(nw, "race")
  # vB <- which(race == "B")
  # vW <- which(race == "W")
  # nB <- length(vB)
  # nW <- length(vW)
  # 
  # num.degrees.B <- length(dist.B)
  # num.degrees.W <- length(dist.W)
  # 
  # deg.B <- apportion_lr(nB, 0:(num.degrees.B - 1), dist.B, shuffled = TRUE)
  # deg.W <- apportion_lr(nW, 0:(num.degrees.W - 1), dist.W, shuffled = TRUE)
  
  race <- get.vertex.attribute(nw, "race")
  agecat2 <- get.vertex.attribute(nw, "agecat2")
  vY <- which(agecat2 == "Y") #creates vector with all id numbers of nodes with that agecat value;
  vO <- which(agecat2 == "O")
  nY <- length(vY)  #total # of nodes w/ Y
  nO <- length(vO)
  nT<-nY+nO
  
  #num.degrees.Y <- length(dist.Y)
  #num.degrees.O <- length(dist.O)
  num.degrees <- length(dist)
  
  #  deg.Y <- apportion_lr(nY, 0:(num.degrees.Y - 1), dist.Y, shuffled = TRUE) #vector of degrees of nodes;
  #  deg.O <- apportion_lr(nO, 0:(num.degrees.O - 1), dist.O, shuffled = TRUE)
  deg <- apportion_lr(nT, 0:(num.degrees - 1), dist, shuffled = TRUE)
  #deg.age <- c(apportion_lr(nY, 0:(num.degrees.Y - 1), dist.Y, shuffled = TRUE),
  #            apportion_lr(nO, 0:(num.degrees.O - 1), dist.O, shuffled = TRUE))
  
  # if (nwstats$method == 2) {
  #   deg.B <- paste0("B", deg.B)
  #   deg.W <- paste0("W", deg.W)
  # }
  # 
  # nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B, v = vB)
  # nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.W, v = vW)
  # 
  # return(nw)
  
  # if agecat=Y then do;
  # deg.age=paste0("Y", deg);
  # else if agecat=Y then do;
  # deg.age=paste0("O", deg);
  
  #deg[agecat2=="Y"] <- paste0("Y", deg)
  #deg[agecat2=="O"] <- paste0("O", deg)
  #deg.age <- paste0("Y", deg.Y) #labels deg.Y as Y+deg (Y0,Y1,Y2) for casual
  #deg.O <- paste0("O", deg.O)
  
  #deg.Y <- paste0("Y", deg.Y)
  #deg.O <- paste0("O", deg.O)
  
  #  if (agecat2 == "Y") {
  #    deg.agept <- paste0("Y", deg)
  # }
  #  if (agecat2 == "O") {
  #    deg.agept<-paste0("O", deg)
  #  }
  
  deg <- paste0(deg)
  
  deg.Y <- deg[agecat2=="Y"]
  deg.O <- deg[agecat2=="O"]
  
  deg.Y <- paste0("Y",deg[agecat2=="Y"])
  deg.O <- paste0("O",deg[agecat2=="O"])
  
  nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg, v = seq_len(network.size(nw))) #attr.name = deg.pers;
  nw <- set.vertex.attribute(nw, attrname = attr.name2, value = deg.Y, v = vY)
  nw <- set.vertex.attribute(nw, attrname = attr.name2, value = deg.O, v = vO)
  #nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.Y, v = vY)
  #nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O, v = vO)
  
  
  return(nw)
}