nstk <- 2
nreg <- 3
n_yr <- 30
n_time <- 2
n_age <- 5
m <- 0.5
t_frac <- c(0.5, 0.5)
r0 <- c(5, 10)
n <- array(0, dim = c(nstk, nreg, n_yr, n_time, n_age))
move <- array(0, dim = c(nstk, nreg, nreg, n_time))

# Stock 1, feeding (time 1) - rows sum to 1
# FROM region 1: 60% stay, 20% each to regions 2&3
move[1, 1, 1, 1] <- 0.6
move[1, 1, 2, 1] <- 0.2
move[1, 1, 3, 1] <- 0.2
# FROM region 2: assume similar
move[1, 2, 1, 1] <- 0.2
move[1, 2, 2, 1] <- 0.6
move[1, 2, 3, 1] <- 0.2
# FROM region 3: assume similar
move[1, 3, 1, 1] <- 0.2
move[1, 3, 2, 1] <- 0.2
move[1, 3, 3, 1] <- 0.6

# Stock 1, natal (time 2) - more dispersal
move[1, 1, 1, 2] <- 1
move[1, 1, 2, 2] <- 0
move[1, 1, 3, 2] <- 0
move[1, 2, 1, 2] <- 1
move[1, 2, 2, 2] <- 0
move[1, 2, 3, 2] <- 0
move[1, 3, 1, 2] <- 1
move[1, 3, 2, 2] <- 0
move[1, 3, 3, 2] <- 0

# Stock 2, feeding (time 1)
move[2, 1, 1, 1] <- 0.65
move[2, 1, 2, 1] <- 0.175
move[2, 1, 3, 1] <- 0.175
move[2, 2, 1, 1] <- 0.175
move[2, 2, 2, 1] <- 0.65
move[2, 2, 3, 1] <- 0.175
move[2, 3, 1, 1] <- 0.175
move[2, 3, 2, 1] <- 0.175
move[2, 3, 3, 1] <- 0.65

# Stock 2, natal (time 2)
move[2, 1, 1, 2] <- 0
move[2, 1, 2, 2] <- 1
move[2, 1, 3, 2] <- 0
move[2, 2, 1, 2] <- 0
move[2, 2, 2, 2] <- 1
move[2, 2, 3, 2] <- 0
move[2, 3, 1, 2] <- 0
move[2, 3, 2, 2] <- 1
move[2, 3, 3, 2] <- 0

# pseudo code:
# 1) recruitment
# 2) movement by time step
# 3) mortality by time step
# 4) project movement and mortality in previous time step to next time step
# 5) if at the end of the year, then project population from last time step to the first time step
# in the next year, and age by one year

for(stk in 1:nstk) {
  for(y in 1:n_yr) {
    for(t in 1:n_time) {

      # Recruitment at age 1, time 1 only
      if(t == 1) {
        if(y == 1) n[stk, stk , y, 1, 1] <- r0[stk]
        else {
          spawners <- sum(n[stk, stk, y-1, n_time, ])
          recruits <- r0[stk] * 0.1 * spawners
          n[stk, stk, y, 1, 1] <- recruits
        }
      }

      # For each age, apply movement and mortality
      for(a in 1:n_age) {
        # Apply movement
        if(a != 1) n[stk, , y, t, a] <- (n[stk, , y, t, a] %*% move[stk, , , t])

        # Apply mortality
        n[stk, , y, t, a] <- n[stk, , y, t, a] * exp(-m * t_frac[t])

        # project stuff to next time step
        if(t != n_time) n[stk, , y, t+1, a] <- n[stk, , y, t, a]
      }

      # At end of year (t == n_time), age fish to next year
      if(t == n_time & y < n_yr) {
        # Age 1 from this year moved to age 2 next year
        n[stk, , y+1, 1, 2:n_age] <- n[stk, , y, t, 1:(n_age-1)]
        # Plus group
        n[stk, , y+1, 1, n_age] <- n[stk, , y+1, 1, n_age] + n[stk, , y, t, n_age]
      }
    }
  }
}

