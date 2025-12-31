# Purpose: To demonstrate the ability to parameterize flexible movement formulations in SPoRC
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 6/13/25

# Setup -------------------------------------------------------------------
library(here)
library(ggplot2)
library(tidyverse)
library(SPoRC)

# Read in 3 area model with 3 age blocks for movement
age_move <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "input_list.RDS"))
age_move_pars <- readRDS(here("dev", "dev_output", "3_Region_Model_Sablefish", "sd_rep.RDS"))
start_movement_vals <- array(age_move_pars$par.fixed[names(age_move_pars$par.fixed) == 'move_pars'], c(3, 2, 3))

# Constant Movement -------------------------------------------------------

# setting up movement parameterization
constant_move <- Setup_Mod_Movement(input_list = age_move,
                                    # Model options
                                    Movement_ageblk_spec = "constant", # estimating constant movement
                                    Movement_yearblk_spec = "constant", # time-invariant for movement
                                    Movement_sexblk_spec = "constant", # sex-invariant movement
                                    do_recruits_move = 0, # recruits do not move
                                    use_fixed_movement = 0, # estimating movement
                                    Use_Movement_Prior = 1, # priors used for movement
                                    Movement_prior = 3, # vague prior to penalize movement away from the extremes
                                    cont_vary_movement = 'none'
)

# Fit tag data by pooling all ages
constant_move$data$move_age_tag_pool <- list(1:length(constant_move$data$ages))

# Change starting values for movement
constant_move$par$move_pars[] <- start_movement_vals[,,1]

# use recruitment proportion prior
constant_move$data$Use_Rec_prop_Prior <- 1
constant_move$data$Rec_prop_prior <- 5

# Fit model
constant_move_obj <- fit_model(data = constant_move$data,
                               parameters = constant_move$par,
                               mapping = constant_move$map,
                               random = NULL,
                               newton_loops = 3,
                               silent = FALSE)

# get sdreport and report file
constant_move_sd_rep <- sdreport(constant_move_obj)
constant_move_rep <- constant_move_obj$rep

# Output files
saveRDS(constant_move_sd_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "constant_move_sd_rep.RDS"))
saveRDS(constant_move_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "constant_move_rep.RDS"))

#  Age Block Movement ---------------------------------------------

# setting up movement parameterization
age_move_block <- Setup_Mod_Movement(input_list = age_move,
                                          # Model options
                                          Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)), # estimating age-varying movement
                                          Movement_yearblk_spec = "constant", # constant movement
                                          Movement_sexblk_spec = "constant", # sex-invariant movement
                                          do_recruits_move = 0, # recruits do not move
                                          use_fixed_movement = 0, # estimating movement
                                          Use_Movement_Prior = 1, # priors used for movement
                                          Movement_prior = 3, # vague prior to penalize movement away from the extremes
                                          cont_vary_movement = 'none'
)

# Fit tag data by pooling by age blocks
age_move_block$data$move_age_tag_pool <- list(c(1:6), c(7:15), c(16:30))

# use recruitment proportion prior
age_move_block$data$Use_Rec_prop_Prior <- 1
age_move_block$data$Rec_prop_prior <- 5

# Fit model
age_move_block_obj <- fit_model(data = age_move_block$data,
                                     parameters = age_move_block$par,
                                     mapping = age_move_block$map,
                                     random = NULL,
                                     newton_loops = 5,
                                     silent = FALSE)

# get sdreport and report file
age_move_block_sd_rep <- sdreport(age_move_block_obj)
age_move_block_rep <- age_move_block_obj$rep

# Output files
saveRDS(age_move_block_sd_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "age_move_block_sd_rep.RDS"))
saveRDS(age_move_block_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "age_move_block_rep.RDS"))


# Time Block Movement -----------------------------------------------------

# setting up movement parameterization
time_move_block <- Setup_Mod_Movement(input_list = age_move,
                                      # Model options
                                      Movement_ageblk_spec = "constant", # estimating constant movement
                                      Movement_yearblk_spec = list(c(1:20), c(21:30), c(31:40), c(41:50), c(51:62)), # time-varying for movement
                                      Movement_sexblk_spec = "constant", # sex-invariant movement
                                      do_recruits_move = 0, # recruits do not move
                                      use_fixed_movement = 0, # estimating movement
                                      Use_Movement_Prior = 1, # priors used for movement
                                      Movement_prior = 3, # vague prior to penalize movement away from the extremes
                                      cont_vary_movement = 'none'
)

# Fit tag data by pooling all ages
time_move_block$data$move_age_tag_pool <- list(1:length(time_move_block$data$ages))

# use recruitment proportion prior
time_move_block$data$Use_Rec_prop_Prior <- 1
time_move_block$data$Rec_prop_prior <- 5

# Fit model
time_move_block_obj <- fit_model(data = time_move_block$data,
                                 parameters = time_move_block$par,
                                 mapping = time_move_block$map,
                                 random = NULL,
                                 newton_loops = 5,
                                 silent = FALSE)

# get sdreport and report file
time_move_block_sd_rep <- sdreport(time_move_block_obj)
time_move_block_rep <- time_move_block_obj$rep

# Output files
saveRDS(time_move_block_sd_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_move_block_sd_rep.RDS"))
saveRDS(time_move_block_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_move_block_rep.RDS"))

# Time and Age Block Movement ---------------------------------------------

# setting up movement parameterization
time_age_move_block <- Setup_Mod_Movement(input_list = age_move,
                                          # Model options
                                          Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)), # estimating age-varying movement
                                          Movement_yearblk_spec = list(c(1:20), c(21:30), c(31:40), c(41:50), c(51:62)), # time-varying for movement
                                          Movement_sexblk_spec = "constant", # sex-invariant movement
                                          do_recruits_move = 0, # recruits do not move
                                          use_fixed_movement = 0, # estimating movement
                                          Use_Movement_Prior = 1, # priors used for movement
                                          Movement_prior = 3, # vague prior to penalize movement away from the extremes
                                          cont_vary_movement = 'none'
)

# Fit tag data by pooling by age blocks
time_age_move_block$data$move_age_tag_pool <- list(c(1:6), c(7:15), c(16:30))

# use recruitment proportion prior
time_age_move_block$data$Use_Rec_prop_Prior <- 1
time_age_move_block$data$Rec_prop_prior <- 5

# Fit model
time_age_move_block_obj <- fit_model(data = time_age_move_block$data,
                                     parameters = time_age_move_block$par,
                                     mapping = time_age_move_block$map,
                                     random = NULL,
                                     newton_loops = 5,
                                     silent = FALSE)

# get sdreport and report file
time_age_move_block_sd_rep <- sdreport(time_age_move_block_obj)
time_age_move_block_rep <- time_age_move_block_obj$rep

# Output files
saveRDS(time_age_move_block_sd_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_age_move_block_sd_rep.RDS"))
saveRDS(time_age_move_block_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_age_move_block_rep.RDS"))

# iid_y_a -------------------------------------------------------------------

# setting up movement parameterization
iid_y_a <- Setup_Mod_Movement(input_list = age_move,
                                 # Model options
                                 Movement_ageblk_spec = "constant", # estimating constant movement
                                 Movement_yearblk_spec = "constant", # time-invariant for movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = 3, # vague prior to penalize movement away from the extremes
                                 cont_vary_movement = 'iid_y_a',
                                 Movement_cont_pe_pars_spec = "fix"
)

# Fit tag data by fitting all ages
iid_y_a$data$move_age_tag_pool <- as.list(1:length(iid_y_a$data$ages))

# Change starting values for movement
iid_y_a$par$move_pe_pars[] <- log(0.5) # fix sigma

# use recruitment proportion prior
iid_y_a$data$Use_Rec_prop_Prior <- 1
iid_y_a$data$Rec_prop_prior <- 5


# Fit model
iid_y_a_obj <- fit_model(data = iid_y_a$data,
                            parameters = iid_y_a$par,
                            mapping = iid_y_a$map,
                            random = NULL,
                            newton_loops = 3,
                            silent = FALSE)

# get sdreport and report file
iid_y_a_sd_rep <- sdreport(iid_y_a_obj)
iid_y_a_rep <- iid_y_a_obj$rep

# Output files
saveRDS(iid_y_a_sd_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "iid_y_a_sd_rep.RDS"))
saveRDS(iid_y_a_rep, here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "iid_y_a_rep.RDS"))



# Plot Comparisons --------------------------------------------------------
# Read in sd reports from all models
age_move_pars <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "age_move_block_sd_rep.RDS"))
constant_pars <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "constant_move_sd_rep.RDS"))
time_move_pars <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_move_block_sd_rep.RDS"))
time_age_move_pars <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_age_move_block_sd_rep.RDS"))
iid_y_a_move_pars <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "iid_y_a_sd_rep.RDS"))

# Read in reports from all models
age_move_rep <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "age_move_block_rep.RDS"))
constant_rep <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "constant_move_rep.RDS"))
time_move_rep <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_move_block_rep.RDS"))
time_age_move_rep <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "time_age_move_block_rep.RDS"))
iid_y_a_move_rep <- readRDS(here("dev", "sporc_manuscript_demonstrations", "movement_param_comparison", "iid_y_a_rep.RDS"))

constant_rep$Movement[,,1,1,1]
age_move_rep$Movement[,,1,1,1]
constant_rep$Movement[,,1,1,1]
time_move_rep$Movement[,,1,1,1]
time_age_move_rep$Movement[,,1,1,1]
iid_y_a_move_rep$Movement[,,1,1,1]

# turn model outputs into list
mod_sd_rep <- list(
  age_block = age_move_pars,
  constant = constant_pars,
  time_block = time_move_pars,
  time_age_block = time_age_move_pars,
  iid_y_a = iid_y_a_move_pars
)

# turn model outputs into list
mod_rep <- list(
  age_block = age_move_rep,
  constant = constant_rep,
  time_block = time_move_rep,
  time_age_block = time_age_move_rep,
  iid_y_a = iid_y_a_move_rep
)

ts_df <- data.frame()
idx_fits_df <- data.frame()
movement_df <- data.frame()

for(i in 1:length(mod_sd_rep)) {

  sd_rep <- mod_sd_rep[[i]] # get sdrep
  rep <- mod_rep[[i]] # get report

  # Regional Spawning Stock Biomass
  ssb_plot_df <- reshape2::melt(array(sd_rep$value[names(sd_rep$value) == 'SSB'], dim = c(3, 62))) %>%
    dplyr::rename(Region = Var1, Year = Var2) %>%
    dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) == "log(SSB)"]) %>%
    dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                  upr = exp(log(value) + 1.96 * se),
                  Region = paste("Region", Region),
                  Type = 'SSB',
                  Name = names(mod_sd_rep)[i])

  depletion_plot_df <- ssb_plot_df %>%
    group_by(Region, Name) %>%
    mutate(value = value / value[1],
           Type = 'Depletion')

  # Recruitment
  rec_plot_df <- reshape2::melt(array(sd_rep$value[names(sd_rep$value) == 'Rec'], dim = c(3, 62))) %>%
    dplyr::rename(Region = Var1, Year = Var2) %>%
    dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) == "log(Rec)"]) %>%
    dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                  upr = exp(log(value) + 1.96 * se),
                  Region = paste("Region", Region),
                  Type = 'Recruitment',
                  Name = names(mod_sd_rep)[i])

  # Get aggregated spawning stock biomass
  agg_ssb_df <- ssb_plot_df %>%
    group_by(Year) %>%
    summarize(value = sum(value)) %>%
    mutate(lwr = NA, upr = NA, se = NA,
           Region = 'Aggregated',
           Type = 'SSB',
           Name = names(mod_sd_rep)[i])

  # get aggregated depletion
  agg_dep_df <- agg_ssb_df %>%
    group_by(Region, Name) %>%
    mutate(value = value / value[1],
           Type = 'Depletion')

  # get aggregated recruitment
  agg_rec_df <- rec_plot_df %>%
    group_by(Year) %>%
    summarize(value = sum(value)) %>%
    mutate(lwr = NA, upr = NA, se = NA,
           Region = 'Aggregated',
           Type = 'Recruitment',
           Name = names(mod_sd_rep)[i])

  # bind together
  ts_df <- rbind(ts_df, ssb_plot_df, rec_plot_df, agg_ssb_df, agg_rec_df, depletion_plot_df, agg_dep_df)

  # get movement
  move_tmp <- reshape2::melt(rep$Movement) %>%
    dplyr::rename(Region_From = Var1, Region_To = Var2, Year = Var3, Age = Var4, Sex = Var5) %>%
    dplyr::mutate(Region_From = paste("From Region", Region_From),
                  Region_To = paste("To Region", Region_To),
                  Sex = paste("Sex", Sex),
                  Name = names(mod_sd_rep)[i]
    )

  movement_df <- rbind(move_tmp, movement_df)

}

# set up color pal
colors <- c("#117733", "black", "#88CCEE", "#44AA99", "#DDCC77")

# Do some more munging
ts_df <- ts_df %>%
  mutate(Region = case_when(
    Region == "Region 1" ~ "BS + AI + WGOA",
    Region == "Region 2" ~ "CGOA",
    Region == "Region 3" ~ "EGOA",
    Region == 'Aggregated' ~ "Aggregated"
  ),
  Name = factor(Name, levels = c("constant", "age_block", "time_block", "time_age_block", "iid_y_a"))) %>%
  group_by(Region, Year) %>%
  mutate(rd = (value - value[Name == 'constant']) / value[Name == 'constant'])

# Time Series Plot
png(here("dev", "sporc_manuscript_demonstrations", "figs", "movement_comp.png"), width = 1000, height = 750)
print(
  ggplot(ts_df %>% filter(Type != 'Depletion'), aes(x = Year + 1959, y = value, color = Name)) +
    geom_line(lwd = 1.3) +
    facet_grid(Type~Region, scales = 'free') +
    scale_color_manual(values = colors) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top') +
    labs(x = 'Year', y = 'Value', fill = 'Model', color = 'Model')
)
dev.off()
