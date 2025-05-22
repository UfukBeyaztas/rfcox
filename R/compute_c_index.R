compute_c_index <- function(time, status, lp) {
  n <- length(time)
  concordant <- 0
  discordant <- 0
  tied <- 0
  comparable <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (status[i] == 1 && time[i] < time[j]) {
        comparable <- comparable + 1
        if (lp[i] > lp[j]) {
          concordant <- concordant + 1
        } else if (lp[i] < lp[j]) {
          discordant <- discordant + 1
        } else {
          tied <- tied + 1
        }
      } else if (status[j] == 1 && time[j] < time[i]) {
        comparable <- comparable + 1
        if (lp[j] > lp[i]) {
          concordant <- concordant + 1
        } else if (lp[j] < lp[i]) {
          discordant <- discordant + 1
        } else {
          tied <- tied + 1
        }
      }
    }
  }
  
  c_index <- (concordant + 0.5 * tied) / comparable
  return(c_index)
}
