
exp_cif <- function(time_grid, event_times, alpha, beta, mu) {
  n <- length(time_grid)
  lambda <- rep(0, n)
  
  for (i in 1:length(time_grid)) {
    time_point <- time_grid[i]
    past_event_times <- event_times[event_times < time_point]
    if (length(past_event_times) == 0) lambda[i] <- mu
    else {
      lambda[i] <- mu + sum(alpha * exp(-beta * (time_point - past_event_times)))
    }
  }
  return(lambda)
}