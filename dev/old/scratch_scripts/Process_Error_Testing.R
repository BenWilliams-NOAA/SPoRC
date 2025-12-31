# purpose: to demonstrate how one might incorporate process error on movement
library(tidyverse)
# read in data
rep <- readRDS("/Users/matthewcheng/Desktop/PostDoc/Spatial Assessments and Sablefish/SPoRC/mlt_rg_sable_mod_results.RDS")

# new movement array to store values in
move_new <- array(0, dim = dim(rep$rep$Movement))

# Approach 1: Dirichlet Draws (Probably no appropriate given Laplace Approximation Assumptions) ---------------
for(r in 1:5) {
  for(y in 1:62) {
    for(a in 1:30) {
      for(s in 1:2) {
        move_new[r,,y,a,s] <- compResidual::rdirichlet(n = 1, alpha = rep$rep$Movement[r,,y,a,s] * 30)
      } # end s loop
    } # end a loop
  } # end y loop
} # end r loop

reshape2::melt(move_new) %>%
  filter(Var4 %in% c(2, 13, 30)) %>%
  ggplot(aes(x = Var3, y = value, color = factor(Var2))) +
  geom_line() +
  facet_grid(Var4~Var1)

# Approach 2: Logit Variables (iid across years and ages) ---------------------------------
constant_move_logit <- array(0, dim = c(5,5))
move_new <- array(0, dim = dim(rep$rep$Movement))

# Populate constant matrix first
ref <- 1 # reference region
for(r in 1:5) {
  logits_row <- log(rep$rep$Movement[r,,1,1,1][-ref] / rep$rep$Movement[r,,1,1,1][ref])
  constant_move_logit[r,] <- c(0, logits_row)
} # end r loop

iid_sigma <- 0.005 * kronecker(diag(1, 62), diag(1, 30))
sim <- mvnfast::rmvn(4, rep(0, 62 * 30), sigma = iid_sigma)
sim <- array(sim, dim = c(4,62,30))

# Generate process error terms (iid)
for(r in 1:5) {
  for(y in 1:62) {
    for(a in 1:30) {
      for(s in 1:2) {
        tmp_move <- constant_move_logit[r,]
        tmp_move[-ref] <- tmp_move[-ref] + sim[,y,a] # adding iid rnorm (bascially a mvn iid)
        move_new[r,,y,a,s] <- exp(tmp_move) / sum(exp(tmp_move))
      } # end s loop
    } # end a loop
  } # end y loop
} # end r loop

reshape2::melt(move_new) %>%
  filter(Var4 %in% c(2, 13, 30)) %>%
  ggplot(aes(x = Var3, y = value, color = factor(Var2))) +
  geom_line() +
  facet_grid(Var4~Var1)

# Approach 3: Logit Variables (iid across years only) ---------------------------------
constant_move_logit <- array(0, dim = c(5,5))
move_new <- array(0, dim = dim(rep$rep$Movement))

# Populate constant matrix first
ref <- 1 # reference region
for(r in 1:5) {
  logits_row <- log(rep$rep$Movement[r,,1,1,1][-ref] / rep$rep$Movement[r,,1,1,1][ref])
  constant_move_logit[r,] <- c(0, logits_row)
} # end r loop

iid_sigma <- 0.005 * diag(1, 62)
sim <- mvnfast::rmvn(4, rep(0, 62), sigma = iid_sigma)
sim <- array(sim, dim = c(4,62))

# Generate process error terms (iid)
for(r in 1:5) {
  for(y in 1:62) {
    for(a in 1:30) {
      for(s in 1:2) {
        tmp_move <- constant_move_logit[r,]
        tmp_move[-ref] <- tmp_move[-ref] + sim[,y] # adding iid rnorm (bascially a mvn iid)
        move_new[r,,y,a,s] <- exp(tmp_move) / sum(exp(tmp_move))
      } # end s loop
    } # end a loop
  } # end y loop
} # end r loop

reshape2::melt(move_new) %>%
  filter(Var4 %in% c(2, 13, 30)) %>%
  ggplot(aes(x = Var3, y = value, color = factor(Var2))) +
  geom_line() +
  facet_grid(Var4~Var1)

# Approach 4: Logit Variables (iid across ages only) ---------------------------------
constant_move_logit <- array(0, dim = c(5,5))
move_new <- array(0, dim = dim(rep$rep$Movement))

# Populate constant matrix first
ref <- 1 # reference region
for(r in 1:5) {
  logits_row <- log(rep$rep$Movement[r,,1,1,1][-ref] / rep$rep$Movement[r,,1,1,1][ref])
  constant_move_logit[r,] <- c(0, logits_row)
} # end r loop

iid_sigma <- 0.005 * diag(1, 30)
sim <- mvnfast::rmvn(4, rep(0, 30), sigma = iid_sigma)
sim <- array(sim, dim = c(4,30))

# Generate process error terms (iid)
for(r in 1:5) {
  for(y in 1:62) {
    for(a in 1:30) {
      for(s in 1:2) {
        tmp_move <- constant_move_logit[r,]
        tmp_move[-ref] <- tmp_move[-ref] + sim[,a] # adding iid rnorm (bascially a mvn iid)
        move_new[r,,y,a,s] <- exp(tmp_move) / sum(exp(tmp_move))
      } # end s loop
    } # end a loop
  } # end y loop
} # end r loop

reshape2::melt(move_new) %>%
  filter(Var3 %in% c(2, 13, 30)) %>%
  ggplot(aes(x = Var4, y = value, color = factor(Var2))) +
  geom_line() +
  facet_grid(Var3~Var1)

# Approach 5: Logit Variables (2dar1; year x age) ---------------------------------

# Define correlation parameters
age_cor <- 0.99  # Correlation between adjacent ages
year_cor <- 0.99  # Correlation between adjacent years

# Create correlation matrices for age and year
age_cormat <- matrix(0, nrow=30, ncol=30)
for(i in 1:30) for(j in 1:30) age_cormat[i,j] <- age_cor^abs(i-j)
year_cormat <- matrix(0, nrow=62, ncol=62)
for(i in 1:62) for(j in 1:62) year_cormat[i,j] <- year_cor^abs(i-j)
spatio_temporal_cor <- kronecker(year_cormat, age_cormat) # get kronecker of both
twodar1_sigma <- 0.25 * spatio_temporal_cor

# Generate correlated errors for all year-age combinations
n_total <- 62 * 30  # Total year-age combinations
sim_all <- mvtnorm::rmvnorm(4, mean=rep(0, 62 * 30), sigma=twodar1_sigma)
sim <- array(0, dim = c(4,62,30))
for(i in 1:4) sim[i,,] <- array(sim_all[i,], c(62,30))

# Generate process error terms (iid)
for(r in 1:5) {
  for(y in 1:62) {
    for(a in 1:30) {
      for(s in 1:2) {
        tmp_move <- constant_move_logit[r,]
        tmp_move[-ref] <- tmp_move[-ref] + sim[,y,a] # adding iid rnorm (bascially a mvn iid)
        move_new[r,,y,a,s] <- exp(tmp_move) / sum(exp(tmp_move))
      } # end s loop
    } # end a loop
  } # end y loop
} # end r loop

reshape2::melt(move_new) %>%
  filter(Var4 %in% c(2, 13, 30)) %>%
  ggplot(aes(x = Var3, y = value, color = factor(Var2))) +
  geom_line() +
  facet_grid(Var4~Var1)

