
#' @method print nwstats
#' @export
print.nwstats <- function(x, ...) {

  cat("Network Statistics Summary")
  cat("\n==========================")
  cat("\nMean degree frequencies by age")
  cat("\nYoung 15-24 (0/1)", rowSums(x$deg.mp.Y))
  cat("\nOlder 25-40 (0/1):", rowSums(x$deg.mp.O))
  cat("\n\nCasual degree frequencies by age")
  cat("\nYoung 15-24 (0/1/2)", colSums(x$deg.mp.Y))
  cat("\nOlder 25-40 (0/1/2):", colSums(x$deg.mp.O))

  cat("\n\nMain network model target statistics:\n")
  cat(x$stats.m)

  cat("\n\nCasual network model target statistics:\n")
  cat(x$stats.p)

  cat("\n\nInstant network model target statistics:\n")
  cat(x$stats.i)

  cat("\n\nMain Model ")
  print(x$coef.diss.m)
  cat("\n\nCasual Model ")
  print(x$coef.diss.p)

}
