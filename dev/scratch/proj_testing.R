# Purpose: To demonstrate the collapsibility of SPoRC using Alaska sablefish
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 6/11/25


# Setup -------------------------------------------------------------------

library(here)
library(ggplot2)
library(tidyverse)
library(SPoRC)
library(geomtextpath)
devtools::load_all(here("R"))

# Read in model outputs
five_dat <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "data.RDS"))
three_dat <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "data.RDS"))
sgl_dat <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "data.RDS"))

five_rep <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "rep.RDS"))
three_rep <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "rep.RDS"))
sgl_rep <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "rep.RDS"))

five_sdrep <- readRDS(here("dev", "dev_output", "5_Region_Model_Sablefish", "sd_rep.RDS"))
three_sdrep <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "sd_rep.RDS"))
sgl_sdrep <- readRDS(here("dev", "dev_output", "1_Region_Model_Sablefish_SptComparison", "sd_rep.RDS"))

# Projections -------------------------------------------------------------

# Define HCR to use
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp # define stock status
  # If stock status is > 1
  if(stock_status >= 1) f <- frp
  # If stock status is between brp and alpha
  if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
  # If stock status is less than alpha
  if(stock_status < alpha) f <- 0
  return(f)
}

# single area model
sgl_rep$h_trans <- 0.8
sgl_ref_pt <- Get_Reference_Points(data = sgl_dat,
                                   rep = sgl_rep,
                                   SPR_x = 0.4,
                                   type = 'single_region',
                                   what = 'SPR',
                                   calc_rec_st_yr = 20,
                                   rec_age = 2)

# set up quantities to use in projection function
single_rg_st <- 62
n_sims <- 1
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 100
n_regions <- 1
n_ages <- length(sgl_dat$ages)
n_sexes <- sgl_dat$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
terminal_NAA <- array(sgl_rep$NAA[,single_rg_st,,], dim = c(n_regions, n_ages, n_sexes))
terminal_NAA0 <- array(sgl_rep$NAA0[,single_rg_st,,], dim = c(n_regions, n_ages, n_sexes))
WAA <- array(rep(sgl_dat$WAA[,single_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(sgl_dat$WAA[,single_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # weight at age
MatAA <- array(rep(sgl_dat$MatAA[,single_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(sgl_rep$fish_sel[,single_rg_st,,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- NULL
terminal_F <- array(sgl_rep$Fmort[,single_rg_st,], dim = c(n_regions, n_fish_fleets))
natmort <- array(sgl_rep$natmort[,single_rg_st,,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes))
recruitment <- array(sgl_rep$Rec[,20:(single_rg_st - 2)], dim = c(n_regions, length(20:(single_rg_st - 2))))
sexratio <- array(0.5, dim = c(n_regions, n_proj_yrs, n_sexes))

# storage
sgl_f_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
sgl_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
sgl_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims))

# do population projection
for(sim in 1:n_sims) {

  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  sexratio = sexratio,
                                  n_fish_fleets = n_fish_fleets,
                                  do_recruits_move = do_recruits_move,
                                  recruitment = recruitment,
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  terminal_F = terminal_F,
                                  natmort = natmort,
                                  WAA = WAA,
                                  WAA_fish = WAA_fish,
                                  MatAA = MatAA,
                                  fish_sel = fish_sel,
                                  Movement = Movement,
                                  f_ref_pt = array(sgl_ref_pt$f_ref_pt, dim = c(sgl_dat$n_regions, n_proj_yrs)),
                                  b_ref_pt = array(sgl_ref_pt$b_ref_pt, dim = c(sgl_dat$n_regions, n_proj_yrs)),
                                  HCR_function = HCR_function,
                                  recruitment_opt = "mean_rec",
                                  fmort_opt = "Input",
                                  t_spawn = t_spawn,
                                  bh_rec_opt = list(
                                    recruitment_dd = 1,
                                    rec_lag = 2,
                                    R0 = sgl_rep$R0,
                                    h = sgl_rep$h_trans,
                                    Rec_Prop = sgl_rep$Rec_trans_prop,
                                    WAA = array(WAA[,1,,1], c(1,30)),
                                    MatAA = array(MatAA[,1,,1], c(1,30)),
                                    natmort = array(natmort[,1,,1], dim = c(1, 30)),
                                    SSB = sgl_rep$SSB,
                                    Movement = sgl_rep$Movement[,,1,,1],
                                    do_recruits_move = sgl_dat$do_recruits_move,
                                    t_spawn = 0,
                                    sex_ratio_f = rep(0.5, 1)

                                  )
  )

  sgl_ssb_proj[,,sim] <- out$proj_SSB
  sgl_catch_proj[,,,sim] <- out$proj_Catch
  sgl_f_proj[,,sim] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

}


# three_rep$Movement <- array(diag(1, 3), dim = c(3,3,62,30,2))
# three_rep$Movement[,,,1:29,] <- diag(1,3)
# three_rep$Movement <- array(1/3, dim = c(3,3,62,30,2))
# Three area model
# three_rep <- readRDS("/Users/matthewcheng/Desktop/constant_move_rep.RDS")
three_rep$h_trans <- c(0.9, 0.9, 0.9)
three_ref_pt <- Get_Reference_Points(data = three_dat,
                                     rep = three_rep,
                                     SPR_x = 0.4,
                                     type = 'multi_region',
                                     what = 'global_BH_MSY',
                                     calc_rec_st_yr = 20,
                                     rec_age = 2
                                     # local_BH_MSY_newton_steps = 10
                                     )

# quantities to use in projection
three_rg_st <- 62
n_sims <- 1
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 100
n_regions <- 3
n_ages <- length(three_dat$ages)
n_sexes <- three_dat$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
terminal_NAA <- array(three_rep$NAA[,three_rg_st,,], dim = c(n_regions, n_ages, n_sexes))
terminal_NAA0 <- array(three_rep$NAA0[,three_rg_st,,], dim = c(n_regions, n_ages, n_sexes))
WAA <- array(rep(three_dat$WAA[,three_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(three_dat$WAA[,single_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # weight at age
MatAA <- array(rep(three_dat$MatAA[,three_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(three_rep$fish_sel[,three_rg_st,,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- aperm(abind::abind(replicate(n_proj_yrs, three_rep$Movement[,,three_rg_st,,], simplify = FALSE), along = 5), perm = c(1,2,5,3,4))
terminal_F <- array(three_rep$Fmort[,three_rg_st,], dim = c(n_regions, n_fish_fleets))
natmort <- array(three_rep$natmort[,three_rg_st,,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes))
recruitment <- array(three_rep$Rec[,20:(three_rg_st - 2)], dim = c(n_regions, length(20:(three_rg_st - 2))))
sexratio <- array(0.5, dim = c(n_regions, n_proj_yrs, n_sexes))

# storage
three_f_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
three_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
three_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims))

# do population projection
for(sim in 1:n_sims) {

  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  sexratio = sexratio,
                                  n_fish_fleets = n_fish_fleets,
                                  do_recruits_move = do_recruits_move,
                                  recruitment = recruitment,
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  terminal_F = terminal_F,
                                  natmort = natmort,
                                  WAA = WAA,
                                  WAA_fish = WAA_fish,
                                  MatAA = MatAA,
                                  fish_sel = fish_sel,
                                  Movement = Movement,
                                  f_ref_pt = array(three_ref_pt$f_ref_pt, dim = c(three_dat$n_regions, n_proj_yrs)),
                                  b_ref_pt = array(three_ref_pt$b_ref_pt, dim = c(three_dat$n_regions, n_proj_yrs)),
                                  HCR_function = HCR_function,
                                  recruitment_opt = "bh_rec",
                                  fmort_opt = "Input",
                                  t_spawn = t_spawn,
                                  bh_rec_opt =  list(
                                    recruitment_dd = 1,
                                    rec_lag = 2,
                                    R0 = three_rep$R0,
                                    h = three_rep$h_trans,
                                    Rec_Prop = three_rep$Rec_trans_prop,
                                    WAA = array(WAA[,1,,1], c(3,30)),
                                    MatAA = array(MatAA[,1,,1], c(3,30)),
                                    natmort = array(natmort[,1,,1], dim = c(3, 30)),
                                    SSB = three_rep$SSB,
                                    Movement = three_rep$Movement[,,1,,1],
                                    do_recruits_move = three_dat$do_recruits_move,
                                    t_spawn = 0,
                                    sex_ratio_f = rep(0.5, 3)
                                  )
  )

  three_ssb_proj[,,sim] <- out$proj_SSB
  three_catch_proj[,,,sim] <- out$proj_Catch
  three_f_proj[,,sim] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

}

# compare b0
sum(three_ref_pt$virgin_b_ref_pt)
out$proj_SSB[,n_proj_yrs]

(sum(three_ref_pt$virgin_b_ref_pt) - sum(out$proj_SSB[,n_proj_yrs])) / sum(out$proj_SSB[,n_proj_yrs]) * 100

# compare r0
three_rep$R0 * three_rep$Rec_trans_prop
rowSums(out$proj_NAA[,n_proj_yrs,1,])

three_ref_pt$virgin_b_ref_pt


local_R0 = three_rep$R0 * three_rep$Rec_trans_prop # get local R0 based on recruitment proportions
SSB = out$proj_SSB[,n_proj_yrs]
S0 = out$proj_SSB[,n_proj_yrs]
h = rep(0.7,3)
(4*h*local_R0*SSB) / ((1-h)*S0 + (5*h-1)*SSB)
local_R0

# three_ref_pt$virgin_b_ref_pt <- c(43.02781698, 115.63662478,  93.55357452)

cbind(three_rep$SSB, out$proj_SSB[])

five_rep$h_trans <- rep(0.9, 5)

# five area model
five_ref_pt <- Get_Reference_Points(data = five_dat,
                                    rep = five_rep,
                                    SPR_x = 0.4,
                                    type = 'multi_region',
                                    what = 'global_BH_MSY',
                                    calc_rec_st_yr = 20,
                                    rec_age = 2
                                    # local_BH_MSY_newton_steps = 10
                                    )

# quantities to use in projection
five_rg_st <- 62
n_sims <- 1
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 100
n_regions <- 5
n_ages <- length(five_dat$ages)
n_sexes <- five_dat$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
terminal_NAA <- array(five_rep$NAA[,length(five_dat$years),,], dim = c(n_regions, n_ages, n_sexes))
terminal_NAA0 <- array(five_rep$NAA0[,length(five_dat$years),,], dim = c(n_regions, n_ages, n_sexes))
WAA <- array(rep(five_dat$WAA[,length(five_dat$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(five_dat$WAA[,single_rg_st,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # weight at age
MatAA <- array(rep(five_dat$MatAA[,length(five_dat$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(five_rep$fish_sel[,length(five_dat$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- aperm(abind::abind(replicate(n_proj_yrs, five_rep$Movement[,,length(five_dat$years),,], simplify = FALSE), along = 5), perm = c(1,2,5,3,4))
terminal_F <- array(five_rep$Fmort[,length(five_dat$years),], dim = c(n_regions, n_fish_fleets))
natmort <- array(five_rep$natmort[,length(five_dat$years),,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes))
recruitment <- array(five_rep$Rec[,20:(length(five_dat$years) - 2)], dim = c(n_regions, length(20:(length(five_dat$years) - 2))))
sexratio <- array(0.5, dim = c(n_regions, n_proj_yrs, n_sexes))

# storage
five_f_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
five_ssb_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_sims))
five_catch_proj <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims))

# do population projection
for(sim in 1:n_sims) {

  # do projection
  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  sexratio = sexratio,
                                  n_fish_fleets = n_fish_fleets,
                                  do_recruits_move = do_recruits_move,
                                  recruitment = recruitment,
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  terminal_F = terminal_F,
                                  natmort = natmort,
                                  WAA = WAA,
                                  WAA_fish = WAA_fish,
                                  MatAA = MatAA,
                                  fish_sel = fish_sel,
                                  Movement = Movement,
                                  f_ref_pt = array(five_ref_pt$f_ref_pt, dim = c(five_dat$n_regions, n_proj_yrs)),
                                  b_ref_pt = array(five_ref_pt$b_ref_pt, dim = c(five_dat$n_regions, n_proj_yrs)),
                                  HCR_function = HCR_function,
                                  recruitment_opt = "bh_rec",
                                  fmort_opt = "Input",
                                  t_spawn = t_spawn,
                                  bh_rec_opt = list(
                                    recruitment_dd = 1,
                                    rec_lag = 2,
                                    R0 = five_rep$R0,
                                    h = five_rep$h_trans,
                                    Rec_Prop = five_rep$Rec_trans_prop,
                                    WAA = array(WAA[,1,,1], c(5,30)),
                                    MatAA = array(MatAA[,1,,1], c(5,30)),
                                    natmort = array(natmort[,1,,1], dim = c(5, 30)),
                                    SSB = five_rep$SSB,
                                    Movement = five_rep$Movement[,,1,,1],
                                    do_recruits_move = five_dat$do_recruits_move,
                                    t_spawn = 0,
                                    sex_ratio_f = rep(0.5, 5)
                                  )
  )

  five_ssb_proj[,,sim] <- out$proj_SSB
  five_catch_proj[,,,sim] <- out$proj_Catch
  five_f_proj[,,sim] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used

}

# plot(out$proj_F[1,])


# Combine Projections -----------------------------------------------------

# bind reports
all_rg_ssb <- reshape2::melt(sgl_rep$SSB) %>%
  mutate(Type = "Single-Area",
         log_se = sgl_sdrep$sd[names(sgl_sdrep$value) == 'log(Aggregated_SSB)']) %>%
  bind_rows(
    reshape2::melt(t(colSums(three_rep$SSB))) %>%
      mutate(Type = 'Three-Area',
             log_se = three_sdrep$sd[names(three_sdrep$value) == 'log(Aggregated_SSB)']),
    reshape2::melt(t(colSums(five_rep$SSB))) %>%
      mutate(Type = 'Five-Area',
             log_se = five_sdrep$sd[names(five_sdrep$value) == 'log(Aggregated_SSB)'])
  ) %>%
  rename(Year = Var2) %>%
  select(-Var1) %>%
  mutate(Type = factor(Type, levels = c("Single-Area", "Three-Area", "Five-Area")),
         Year = Year + 1959)

# Five region SSB data
five_rg_ssb <- reshape2::melt(five_rep$SSB) %>%
  rename(Region = Var1, Year = Var2) %>%
  left_join(data.frame(Region = 1:5, b40 = five_ref_pt$b_ref_pt), by = 'Region') %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS',
    Region == 2 ~ 'AI',
    Region == 3 ~ 'WGOA',
    Region == 4 ~ 'CGOA',
    Region == 5 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = five_sdrep$sd[names(five_sdrep$value) == 'log(SSB)'],
         Year = Year + 1959)

# Three region SSB data
three_rg_ssb <- reshape2::melt(three_rep$SSB) %>%
  rename(Region = Var1, Year = Var2) %>%
  left_join(data.frame(Region = 1:3, b40 = three_ref_pt$b_ref_pt), by = 'Region') %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS + AI + WGOA',
    Region == 2 ~ 'CGOA',
    Region == 3 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS + AI + WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = three_sdrep$sd[names(three_sdrep$value) == 'log(SSB)'],
         Year = Year + 1959)

# bind five region historical estimates
five_rg_ssb <- five_rg_ssb %>% mutate(Sim = 1, Type = 'Historical') %>%
  bind_rows(reshape2::melt(five_ssb_proj) %>%
              rename(Region = Var1, Year = Var2, Sim = Var3) %>%
              dplyr::mutate(Region = dplyr::case_when(
                Region == 1 ~ 'BS',
                Region == 2 ~ 'AI',
                Region == 3 ~ 'WGOA',
                Region == 4 ~ 'CGOA',
                Region == 5 ~ 'EGOA'
              ),
              Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
              log_se = NA,
              Year = Year + 2020,
              Type = 'Projection'))

# bind three region historical estimates
three_rg_ssb <- three_rg_ssb %>% mutate(Sim = 1, Type = 'Historical') %>%
  bind_rows(reshape2::melt(three_ssb_proj) %>%
              rename(Region = Var1, Year = Var2, Sim = Var3) %>%
              dplyr::mutate(Region = dplyr::case_when(
                Region == 1 ~ 'BS + AI + WGOA',
                Region == 2 ~ 'CGOA',
                Region == 3 ~ 'EGOA'
              ),
              Region = factor(Region, levels = c("BS + AI + WGOA", "CGOA", "EGOA")),
              log_se = NA,
              Year = Year + 2020,
              Type = 'Projection'))

# bind aggregated historical estimates
all_rg_ssb <- all_rg_ssb %>% mutate(Sim = 1, Type_Period = 'Historical') %>%
  bind_rows(

    reshape2::melt(sgl_ssb_proj) %>%
      rename(Region = Var1, Year = Var2, Sim = Var3) %>%
      dplyr::mutate(,
                    log_se = NA,
                    Year = Year + 2020,
                    Type_Period = 'Projection',
                    Type = 'Single-Area') %>%
      group_by(Year, Sim, Type_Period, Type) %>%
      summarize(value = sum(value)),

    reshape2::melt(three_ssb_proj) %>%
      rename(Region = Var1, Year = Var2, Sim = Var3) %>%
      dplyr::mutate(,
                    log_se = NA,
                    Year = Year + 2020,
                    Type_Period = 'Projection',
                    Type = 'Three-Area') %>%
      group_by(Year, Sim, Type_Period, Type) %>%
      summarize(value = sum(value)),

    reshape2::melt(five_ssb_proj) %>%
      rename(Region = Var1, Year = Var2, Sim = Var3) %>%
      dplyr::mutate(,
                    log_se = NA,
                    Year = Year + 2020,
                    Type_Period = 'Projection',
                    Type = 'Five-Area') %>%
      group_by(Year, Sim, Type_Period, Type) %>%
      summarize(value = sum(value))

  ) %>%
  left_join(data.frame(Type = c("Single-Area", "Three-Area", "Five-Area"),
                       b40 = c(sgl_ref_pt$b_ref_pt, sum(three_ref_pt$b_ref_pt), sum(five_ref_pt$b_ref_pt)))) %>%
  mutate(Type = factor(Type, levels = c("Single-Area", "Three-Area", "Five-Area")))

# Three region Rec data
three_rg_rec <- reshape2::melt(three_rep$Rec) %>%
  rename(Region = Var1, Year = Var2) %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS + AI + WGOA',
    Region == 2 ~ 'CGOA',
    Region == 3 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS + AI + WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = three_sdrep$sd[names(three_sdrep$value) == 'log(Rec)'],
         Year = Year + 1959)

# Five region rec data
five_rg_rec <- reshape2::melt(five_rep$Rec) %>%
  rename(Region = Var1, Year = Var2) %>%
  dplyr::mutate(Region = dplyr::case_when(
    Region == 1 ~ 'BS',
    Region == 2 ~ 'AI',
    Region == 3 ~ 'WGOA',
    Region == 4 ~ 'CGOA',
    Region == 5 ~ 'EGOA'
  ),
  Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  mutate(log_se = five_sdrep$sd[names(five_sdrep$value) == 'log(Rec)'],
         Year = Year + 1959)

# Plots -------------------------------------------------------------------

## Regional Biomass --------------------------------------------------------
# five region plot
five_ssb_plot <- ggplot() +
  geom_line(five_rg_ssb %>% filter(Type == 'Historical'),
            mapping = aes(x = Year, y = value), color = '#00a473', lty = 1, lwd = 1.3) +
  geom_ribbon(five_rg_ssb %>% filter(Type == 'Historical'),
              mapping = aes(x = Year, y = value, ymin = exp(log(value) - 1.96 * log_se),
                            ymax = exp(log(value) + 1.96 * log_se)), alpha = 0.25, color = NA, fill = '#00a473') +
  geom_line(five_rg_ssb %>% filter(Type == 'Projection'),
            mapping = aes(x = Year, y = value, group = Sim), color = '#00a473', lty = 1, lwd = 1.3, alpha = 1) +
  geom_hline(five_rg_ssb, mapping = aes(yintercept = b40), lty = 2, lwd = 1, color = '#00a473') +
  geom_vline(xintercept = 2021, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~Region, nrow = 5) +
  labs(x = 'Year', y = '') +
  theme_bw(base_size = 25)

# three region plot
three_ssb_plot <- ggplot() +
  geom_line(three_rg_ssb %>% filter(Type == 'Historical'),
            mapping = aes(x = Year, y = value), color = '#ef5f10', lty = 1, lwd = 1.3) +
  geom_ribbon(three_rg_ssb %>% filter(Type == 'Historical'),
              mapping = aes(x = Year, y = value, ymin = exp(log(value) - 1.96 * log_se),
                            ymax = exp(log(value) + 1.96 * log_se)), alpha = 0.25, color = NA, fill = '#ef5f10') +
  geom_line(three_rg_ssb %>% filter(Type == 'Projection'),
            mapping = aes(x = Year, y = value, group = Sim), color = '#ef5f10', lty = 1, lwd = 1.3, alpha = 1) +
  geom_hline(three_rg_ssb, mapping = aes(yintercept = b40), lty = 2, lwd = 1, color = '#ef5f10') +
  geom_vline(xintercept = 2021, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~Region, nrow = 3) +
  labs(x = 'Year', y = 'Spawning Stock Biomass ') +
  theme_bw(base_size = 25)

# aggregated plot
all_ssb_plot <- ggplot() +
  geom_line(all_rg_ssb %>% filter(Type_Period == 'Historical'),
            mapping = aes(x = Year, y = value, color = Type), lty = 1, lwd = 1.3) +
  geom_ribbon(all_rg_ssb %>% filter(Type_Period == 'Historical'),
              mapping = aes(x = Year, y = value, ymin = exp(log(value) - 1.96 * log_se),
                            ymax = exp(log(value) + 1.96 * log_se), fill = Type), alpha = 0.25, color = NA) +
  geom_line(all_rg_ssb %>% filter(Type_Period == 'Projection'),
            mapping = aes(x = Year, y = value, group = interaction(Sim, Type), color = Type), lty = 1, lwd = 1.3, alpha = 1) +
  geom_hline(all_rg_ssb, mapping = aes(yintercept = b40, color = Type), lty = 2, lwd = 1) +
  geom_vline(xintercept = 2021, lty = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = c("#7e79b0", "#ef5f10", "#00a473")) +
  scale_fill_manual(values = c("#7e79b0", "#ef5f10", "#00a473")) +
  labs(x = 'Year', y = 'Spawning Stock Biomass ', fill = 'Model', color = 'Model') +
  theme_bw(base_size = 25) +
  theme(legend.position = c(0.895, 0.125),
        legend.background = element_blank())


# combine plots
multi_ssb_plots <- cowplot::plot_grid(three_ssb_plot, five_ssb_plot, ncol = 2,
                                      labels = c("B", "C"), label_size = 30, label_x = 0.03)

all_ssb_plot
# three_ssb_plot
# five_ssb_plot

