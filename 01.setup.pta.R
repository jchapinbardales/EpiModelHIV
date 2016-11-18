## Partner Type By Age SETUP ##

rm(list = ls()) 
suppressMessages(library(EpiModelHIV)) 

# Time unit for simulation, relative to 1 day 
time.unit <- 7 


# Population size by race 
num.B <- 5000 
num.W <- 5000 

num.Y <- 3693 
num.O <- 6307
#must add to 10000, but size in each will differ...
#36.93% 18-24
#63.07% 25-40
# so...10000*.3693=3693
#      10000*.6307=6307

#done 11/1/2016

# 
# # mean/pers degree distributions matrices. 
# deg.mp.B <- deg.mp.W <- (matrix(c(0.506, 0.151, 0.053, 0.207, 0.061, 0.022),  #.15157=.152 and .054 in notes
#                                 byrow = TRUE, nrow = 2) + 
#                          matrix(c(0.435, 0.184, 0.095, 0.233, 0.033, 0.020),  #.021 in notes
#                                 byrow = TRUE, nrow = 2))/2 

################## PT and AGE ############################
#Y=18-24
#O=25-40

deg.mp.Y <- matrix(c(0.509, 0.169, 0.054, 0.202, 0.045, 0.021),  #% 0MAIN0CAS,0MAIN1CAS,0MAIN2CAS,1MAIN0CAS,1MAIN1CAS,1MAIN2CAS --> sum to 1;
                   byrow = TRUE, nrow = 2)

deg.mp.O <- matrix(c(0.455, 0.164, 0.081, 0.228, 0.051, 0.021),  
                   byrow = TRUE, nrow = 2)

#overall deg.mp for fitting model;
deg.mp.overall <-   (matrix(c(0.509, 0.169, 0.054, 0.202, 0.045, 0.021),  
                              byrow = TRUE, nrow = 2) + 
                     matrix(c(0.455, 0.164, 0.081, 0.228, 0.051, 0.021),  
                              byrow = TRUE, nrow = 2))/2

#done 11/1/2016
# mean deg not captured by mean(deg.mp.Y)=0.166667, that just means the mean of 1/6 groups (ie. sum of 1/6 groups)
# done in excel and got: 0.632=mean deg for Y & 0.719=mean degree for O, slightly higher in older, more 2Cs, more 1M
# Mean degree frequencies by age
# Young 15-24 (0/1) 0.732 0.268
# Older 25-40 (0/1): 0.7 0.3
# 
# Casual degree frequencies by age
# Young 15-24 (0/1/2) 0.711 0.214 0.075
# Older 25-40 (0/1/2): 0.683 0.215 0.102


# deg.mp.B
# deg.mp.B = deg.mp.W
#      [,1]   [,2]    [,3]
# [1,] 0.4705 0.1675  0.074
# [2,] 0.2200 0.0470  0.021
# = 
#        0 CAS  1 CAS   2 CAS
# 0 MAIN 0.4705 0.1675  0.074
# 1 MAIN 0.2200 0.0470  0.021
# and there is constraint that no more than 3 partners and
# and no more than 1 main, so either have main or don't
# and if main, then no more than 2 other casual partners

# which comes from: 
#White or Black?
#       [,1]  [,2]  [,3]
# [1,] 0.506 0.151 0.053
# [2,] 0.207 0.061 0.022
# +
#White or Black?
#       [,1]  [,2]  [,3]
# [1,] 0.435 0.184 0.095
# [2,] 0.233 0.033 0.020
# =
# Marginal over B & W:
#       [,1]  [,2]  [,3]
# [1,] 0.941 0.335 0.148
# [2,] 0.440 0.094 0.042
# /2
#TO GET AVERAGE
# =
#       [,1]   [,2]    [,3]
# [1,] 0.4705 0.1675  0.074
# [2,] 0.2200 0.0470  0.021
# then say that this is the same matrix for W and B

# each of these is a proportion of the population that
# falls into one of these six categories (prop of all MSM that had 1M2C, etc)
# to get target statistics for edges for main model, 
# take total marginal main degree = (.22+.047+0.021)*10000 pop size / 2 to get edges
# = 1440 -- this is first target statistic in main network;
# take total marginal casual degree = (.1675+.0470+2(.074)+2(.021))*10000/2
# = 2022.5 -- this is first target statistic in casual network;
## now where does 470, 210, and 667.68 come from?? 
# 470 = .0470*10000 = # nodes with 1 Main, 1 Cas
# 210 = .0210*10000 = # nodes with 1 Main, 2 Cas
# 667.8 ?
#  if (method == 1) {
#     stats.m <- c(edges.m, totdeg.m.by.dp[2:3], sqrt.adiff.m)
#  }                                      #col 2, col 3?
#  if (method == 1) {
#     totdeg.m.by.dp <- c(num * deg.mp.B[2, ])
#  }                      
#  if (method == 1) {
#     stats.m <- c(edges.m, totdeg.m.by.dp[2:3], sqrt.adiff.m)
#  }                        =10000*.047 for row2,col2; 10000*.021 for row2,col3)
#  if (method == 1) {
#     sqrt.adiff.m <- edges.m * mean(c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1]))
#  }

## Now where does 890, 950, 1185.859 come from?? 
# from estimation.R:
#  if (method == 1) {
#     stats.p <- c(edges.p, totdeg.p.by.dm[2], conc.p, sqrt.adiff.p)
#  }                                      #row 2? yes, across those w/ 1 Main
#  if (method == 1) {
#     totdeg.p.by.dm <- c(num * deg.mp.B[, 2] + num * deg.mp.B[, 3] * 2)
#  }                      =10000*.047 + 1000*.021 *2 = 890
#  conc.p = num concurrent casual?=.074+.021=0.095*10000=950  
#  if (method == 1) {
#     conc.p <- sum(deg.mp.B[, 3] * num)
#  }
#  if (method == 1) {
#     sqrt.adiff.p <- edges.p * mean(c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2]))
#  }

# Instant rates 
# mdeg.inst.B <- mdeg.inst.W <- (matrix(c(0.010402, 0.012954, 0.011485, 0.007912, 0.007424, 0.007424),  
#                                       byrow = TRUE, nrow = 2) + 
#                                matrix(c(0.008186, 0.012017, 0.013024, 0.008151, 0.008341, 0.008341),  
#                                       byrow = TRUE, nrow = 2))/2 

################## PT and AGE ############################
#Y=18-24
#O=25-40
#Daily prob of one-off partnership given 0MAIN0CAS,0MAIN1CAS,0MAIN2CAS,1MAIN0CAS,1MAIN1CAS,1MAIN2CAS ;

mdeg.inst.Y <- matrix(c(0.008519, 0.010380, 0.012908, 0.008292, 0.011870, 0.011870),  
                      byrow = TRUE, nrow = 2)
mdeg.inst.O <- matrix(c(0.010280, 0.013939, 0.011883, 0.008336, 0.006404, 0.006404),  
                      byrow = TRUE, nrow = 2)

mdeg.inst <- (matrix(c(0.008519, 0.010380, 0.012908, 0.008292, 0.011870, 0.011870),  
                      byrow = TRUE, nrow = 2) + 
              matrix(c(0.010280, 0.013939, 0.011883, 0.008336, 0.006404, 0.006404),  
                      byrow = TRUE, nrow = 2))/2


#Daily probs of one-off partnership:
#Appendix Note: We can take the total number of AI partners reported for 6 months, 
#multiply this number by the proportion of partners that are one time in the dyadic section, 
#and then divide by 180 to get back to daily.

#Black
#      [,1]   [,2]   [,3]
# [1,] 0.0104 0.0129 0.0115
# [2,] 0.0079 0.0074 00074    ##repeated for M1,C1 and M1, C2+ -- M1, C2+ not in Appendix
# +
#White
#      [,1]   [,2]   [,3]
# [1,] 0.0082 0.0120 0.0130
# [2,] 0.0082 0.0083 0.0083   ##repeated for M1,C1 and M1, C2+ -- M1, C2+ not in Appendix


# Quintile distribution of overall AI rates - heterogeneity in one offs; 
# rank people based on their rates of AI, from lowest AI rate to highest
# divide into quintiles (5 equal groups)
# take average rate of AI in each group = these 5 values below
# accounts for the fact that some people just aren't having sex with one-offs
# they probably have sex with main/casual, but just aren't going to have one-offs
# so their AI rate is effectively 0 -- this is brought in and used in one-off ergm;

qnts.W <- qnts.B <- c(0,  
                      0.00100413270621469,  
                      0.00536670889830508, 
                      0.0101956474689266,  
                      0.0315210354777778) 


# Proportion in same-race partnerships (main, casl, inst) 
prop.hom.mpi.B <- prop.hom.mpi.W <- (c(0.9484, 0.9019, 0.9085) +  
                                       c(0.9154, 0.8509, 0.8944))/2 

# Mean age diffs (main, casl, inst) 
sqrt.adiff.BB <- c(0.417, 0.498, 0.456) 
sqrt.adiff.BW <- c(0.454, 0.629, 0.585) 
sqrt.adiff.WW <- c(0.520, 0.632, 0.590) 


# Mean durations 
# rates.main <- mean(c(0.00287612937991679,  #BBMAIN
#                      0.00269183371091241,  #BWMAIN
#                      0.00180272348650181)) #WWMAIN
# 
# rates.pers <- mean(c(0.00761700198417522,  #BBCAS ## ex: medianD=91 days, ln2/rate=91, rate=ln2/91
#                      0.00350074333616134,  #BWCAS
#                      0.00693147180559945)) #WWCAS
# 
# durs.main <- 1/rates.main 
# durs.pers <- 1/rates.pers 
##see final duration resolution 7_30_2015 email with excel attachment;

################## PT and AGE ############################
#Y=18-24
#O=25-40
rates.main <- mean(c(0.00787667250636301,  #YYMAIN
                     0.00285833888890699,  #YOMAIN
                     0.00182407152778933)) #OOMAIN
rates.pers <- mean(c(0.00682903626167434,  #YYCAS
                     0.00835117085011982,  #YOCAS
                     0.00572848909553674)) #OOCAS

durs.main <- 1/rates.main  ## should be log2/rates.main.age if separated by age?
durs.pers <- 1/rates.pers  ## should be log2/rates.main.age if separated by age?

#=mean of median durations across ages, one duration for main, one duration for casual;
#durs.main=238.8709 (YY 126, YO|OY 349, OO 548)
#durs.pers=143.481 (YY 146, YO|OY 119, OO 174)

## rates are based on median not mean, so shouldn't transformation to duration
## be ln2/rate to hit the median? not 1/rate if going to have separated by age
## and not taking mean of medians?

# durations separated out by age
# rates.main.YY <- 0.00787667250636301  #YYMAIN
# rates.main.YO <- 0.00285833888890699  #YOMAIN
# rates.main.OO <- 0.00182407152778933  #OOMAIN
# rates.pers.YY <- 0.00682903626167434  #YYCAS
# rates.pers.YO <- 0.00835117085011982  #YOCAS
# rates.pers.OO <- 0.00572848909553674  #OOCAS
# 
# durs.main.YY <- 1/rates.main.YY
# durs.main.YO <- 1/rates.main.YO
# durs.main.OO <- 1/rates.main.OO
# durs.pers.YY <- 1/rates.pers.YY
# durs.pers.YO <- 1/rates.pers.YO
# durs.pers.OO <- 1/rates.pers.OO


# Age-sex-specific mortality rates 
ages <- 18:39 
asmr.B <- c(rep(0, 17), 
            1-(1-c(rep(0.00159, 7), 
                   rep(0.00225, 10), 
                   rep(0.00348, 5)))^(1/(365/time.unit)), 1) 


asmr.W <- c(rep(0, 17), 
            1-(1-c(rep(0.00103, 7), 
                   rep(0.00133, 10), 
                   rep(0.00214, 5)))^(1/(365/time.unit)), 1) 


# I, R, V role frequencies 
# role.B.prob <- role.W.prob <- 
# (c(0.242, 0.321, 0.437) + c(0.228, 0.228, 0.544))/2  #should be .543 not .544 but whatevs

################## PT and AGE ############################
#Y=18-24
#O=25-40
role.Y.prob <- c(0.13964, 0.32883, 0.53153)  #I, R, V
role.O.prob <- c(0.2941, 0.2594, 0.4465)  #I, R, V 

role.prob.overall <- (c(0.13964, 0.32883, 0.53153) + c(0.2941, 0.2594, 0.4465))/2

#done 11/1/2016



# Create meanstats 
st <- calc_nwstats_msm( 
  #method = 1, 
  time.unit = time.unit, 
  #num.B = num.B, 
  #num.W = num.W, 
  num.Y = num.Y, 
  num.O = num.O, 
  #deg.mp.B = deg.mp.B, 
  #deg.mp.W = deg.mp.W,
  deg.mp.Y = deg.mp.Y,
  deg.mp.O = deg.mp.O,
  #mdeg.inst.B = mdeg.inst.B, 
  #mdeg.inst.W = mdeg.inst.W, 
  mdeg.inst.Y = mdeg.inst.Y,
  mdeg.inst.O = mdeg.inst.O,
  qnts.B = qnts.B, 
  qnts.W = qnts.W, 
  prop.hom.mpi.B = prop.hom.mpi.B, 
  prop.hom.mpi.W = prop.hom.mpi.W, 
  balance = "mean", 
  sqrt.adiff.BB = sqrt.adiff.BB, 
  sqrt.adiff.WW = sqrt.adiff.WW, 
  sqrt.adiff.BW = sqrt.adiff.BW, 
  diss.main = ~offset(edges), 
  diss.pers = ~offset(edges), 
  durs.main = durs.main, 
  durs.pers = durs.pers, 
  #durs.main.age = durs.main.age,
  #durs.pers.age = durs.pers.age,
  ages = ages, 
  asmr.B = asmr.B, 
  asmr.W = asmr.W, 
  #role.B.prob = role.B.prob, 
  #role.W.prob = role.W.prob #,
  role.Y.prob = role.Y.prob,
  role.O.prob = role.O.prob
) 

st

#Save summary statistics data object "st"
save(st, file = "C:/Users/jchapi2/Documents/GitHub/EpiModelHIV/est/nwstats.rda") #("est/nwstats.rda")
rm(list = ls()) 

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
st$stats.m
st$totdeg.p.by.dm

st$asmr.B
st$durs.main
st$stats.m
st$stats.p
?calc_nwstats_msm

#ST, AGE
# Network Statistics Summary
# ==========================
#   Mean degree frequencies by age
# Young 15-24 (0/1) 0.732 0.268
# Older 25-40 (0/1): 0.7 0.3
# 
# Casual degree frequencies by age
# Young 15-24 (0/1/2) 0.711 0.214 0.075
# Older 25-40 (0/1/2): 0.683 0.215 0.102
# 
# Main network model target statistics:
#   1440.912 166.185 77.553 321.657 132.447 668.1029
# 
# Casual network model target statistics:
#   1993.442 321.291 586.551 276.975 643.314 1168.822
# 
# Instant network model target statistics:
#   346.6009 43.30001 45.34834 13.80831 18.01897 6.443879 206.5025 83.90994 100.9244 14.41924 42.49443 5.937334 0 14.05786 142.7391 441.2945 188.4354
# 
# Main Model Dissolution Coefficients
# =======================
#   Dissolution Model: ~offset(edges)
# Target Statistics: 34.12442
# Crude Coefficient: 3.500271
# Mortality/Exit Rate: 3.588345e-05
# Adjusted Coefficient: 3.502723
# 
# Casual Model Dissolution Coefficients
# =======================
#   Dissolution Model: ~offset(edges)
# Target Statistics: 20.49728
# Crude Coefficient: 2.970275
# Mortality/Exit Rate: 3.588345e-05
# Adjusted Coefficient: 2.971747













# RACE Looks like this:
# Network Statistics Summary
# ==========================
#   Mean degree frequencies by race
# Black (0/1) 0.712 0.288
# White (0/1): 0.712 0.288
# 
# Casual degree frequencies by race
# Black (0/1/2) 0.6905 0.2145 0.095
# White (0/1/2): 0.6905 0.2145 0.095
# 
# Main network model target statistics:
#   1440 470 210 667.68
# 
# Casual network model target statistics:
#   2022.5 890 950 1185.859
# 
# Instant network model target statistics:
#   338.5872 123.6851 146.3925 25.93342 63.47831 11.58727 0 14 142.8 441 184.0786

# IE.  With nothing specified in the dissolution models (just offset), based
#      on the main duration/casual duration, this is what we should expect?
# Main Model Dissolution Coefficients
# =======================
#   Dissolution Model: ~offset(edges)
# Target Statistics: 58.14292
# Crude Coefficient: 4.045555
# Mortality/Exit Rate: 3.588345e-05
# Adjusted Coefficient: 4.049737
# 
# Casual Model Dissolution Coefficients
# =======================
#   Dissolution Model: ~offset(edges)
# Target Statistics: 23.74488
# Crude Coefficient: 3.12434
# Mortality/Exit Rate: 3.588345e-05
# Adjusted Coefficient: 3.126046
