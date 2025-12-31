library(here)
devtools::load_all(here("R"))
library(tidyverse)
library(SPoRC)
# Load in data
load('/Users/matthewcheng/Downloads/YFT_4area_observations_1_100_ESS_05.RData')

# Define dimensions
n_regions <- dat_4A_2$N_areas
yrs <- dat_4A_2$styr:dat_4A_2$endyr
n_ages <- dat_4A_2$Nages
n_fish_fleets <- 7
n_lens <- dat_4A_2$N_lbins
n_sexes <- 1

# Setup observed catch arrays

# Fleet 1: gi
# Fleet 2: hd
# Fleet 3: ll
# Fleet 4: other
# Fleet 5: bb
# Fleet 6: ps
# Fleet 7: trol

# Observed catches
ObsCatch <- array(0, dim = c(n_regions, length(yrs), n_fish_fleets))
ObsCatch[1,,1] <- dat_4A_1$catch$fishing_gi_1
ObsCatch[4,,1] <- dat_4A_1$catch$fishing_gi_4
ObsCatch[1,,2] <- dat_4A_1$catch$fishing_hd_1
ObsCatch[1,,3] <- dat_4A_1$catch$fishing_ll_1
ObsCatch[2,,3] <- dat_4A_1$catch$fishing_ll_2
ObsCatch[3,,3] <- dat_4A_1$catch$fishing_ll_3
ObsCatch[4,,3] <- dat_4A_1$catch$fishing_ll_4
ObsCatch[1,,4] <- dat_4A_1$catch$fishing_other_1
ObsCatch[4,,4] <- dat_4A_1$catch$fishing_other_4
ObsCatch[1,,5] <- dat_4A_1$catch$fishing_bb_1
ObsCatch[1,,6] <- dat_4A_1$catch$fishing_ps_1
ObsCatch[2,,6] <- dat_4A_1$catch$fishing_ps_2
ObsCatch[4,,6] <- dat_4A_1$catch$fishing_ps_4
ObsCatch[1,,7] <- dat_4A_1$catch$fishing_trol_1
ObsCatch[2,,7] <- dat_4A_1$catch$fishing_trol_2
ObsCatch[4,,7] <- dat_4A_1$catch$fishing_trol_4

# Use Catch Indicator
UseCatch <- array(0, dim = c(n_regions, length(yrs), n_fish_fleets))
UseCatch[which(ObsCatch != 0)] <- 1

# Setup fishery indices
dat_4A_1$CPUE <- dat_4A_1$CPUE %>% mutate(year = as.numeric(paste(year)))
ObsFishIdx <- array(0, dim = c(n_regions, length(yrs), n_fish_fleets))
ObsFishIdx[1,(dat_4A_1$CPUE %>% filter(index == 17))$year,3] <- (dat_4A_1$CPUE %>% filter(index == 17))$cpu
ObsFishIdx[2,(dat_4A_1$CPUE %>% filter(index == 18))$year,3] <- (dat_4A_1$CPUE %>% filter(index == 18))$cpu
ObsFishIdx[3,(dat_4A_1$CPUE %>% filter(index == 19))$year,3] <- (dat_4A_1$CPUE %>% filter(index == 19))$cpu
ObsFishIdx[4,(dat_4A_1$CPUE %>% filter(index == 20))$year,3] <- (dat_4A_1$CPUE %>% filter(index == 20))$cpu

# Use fishery index indicator
UseFishIdx <- array(0, dim = c(n_regions, length(yrs), n_fish_fleets))
UseFishIdx[which(ObsFishIdx != 0)] <- 1

# Setup fishery compositions
ObsFishLenComps <- array(0, dim = c(n_regions, length(yrs), n_lens, n_sexes, n_fish_fleets))
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 1))$Yr,,1,1] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 1))[,-c(1:6)])
ObsFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 2))$Yr,,1,1] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 2))[,-c(1:6)])
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 3))$Yr,,1,2] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 3))[,-c(1:6)])
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 4))$Yr,,1,3] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 4))[,-c(1:6)])
ObsFishLenComps[2,(dat_4A_1$lencomp %>% filter(FltSvy == 5))$Yr,,1,3] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 5))[,-c(1:6)])
ObsFishLenComps[3,(dat_4A_1$lencomp %>% filter(FltSvy == 6))$Yr,,1,3] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 6))[,-c(1:6)])
ObsFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 7))$Yr,,1,3] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 7))[,-c(1:6)])
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 8))$Yr,,1,4] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 8))[,-c(1:6)])
ObsFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 9))$Yr,,1,4] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 9))[,-c(1:6)])
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 10))$Yr,,1,5] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 10))[,-c(1:6)])
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 11))$Yr,,1,6] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 11))[,-c(1:6)])
ObsFishLenComps[2,(dat_4A_1$lencomp %>% filter(FltSvy == 12))$Yr,,1,6] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 12))[,-c(1:6)])
ObsFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 13))$Yr,,1,6] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 13))[,-c(1:6)])
ObsFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 14))$Yr,,1,7] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 14))[,-c(1:6)])
ObsFishLenComps[2,(dat_4A_1$lencomp %>% filter(FltSvy == 15))$Yr,,1,7] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 15))[,-c(1:6)])
ObsFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 16))$Yr,,1,7] <- as.matrix((dat_4A_1$lencomp %>% filter(FltSvy == 16))[,-c(1:6)])

# Use Fishery Length Indicators
fleet10_yrs <- (dat_4A_1$lencomp %>% filter(FltSvy == 10))$Yr
fleet12_yrs <- (dat_4A_1$lencomp %>% filter(FltSvy == 12))$Yr
fleet13_yrs <- (dat_4A_1$lencomp %>% filter(FltSvy == 13))$Yr

fleet12_yrs
fleet12_yrs[!fleet12_yrs %in% which(UseCatch[2,,6] == 0)]
fleet13_yrs[!fleet13_yrs %in% which(UseCatch[4,,6] == 0)]

UseFishLenComps <- array(0, dim = c(n_regions, length(yrs), n_fish_fleets))
UseFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 1))$Yr,1]  <- 1
UseFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 2))$Yr,1]  <- 1
UseFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 3))$Yr,2]  <- 1
UseFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 4))$Yr,3]  <- 1
UseFishLenComps[2,(dat_4A_1$lencomp %>% filter(FltSvy == 5))$Yr,3]  <- 1
UseFishLenComps[3,(dat_4A_1$lencomp %>% filter(FltSvy == 6))$Yr,3]  <- 1
UseFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 7))$Yr,3]  <- 1
UseFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 8))$Yr,4]  <- 1
UseFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 9))$Yr,4]  <- 1
UseFishLenComps[1,fleet10_yrs[!fleet10_yrs %in% which(UseCatch[1,,5] == 0)],5] <- 1
UseFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 11))$Yr,6] <- 1
UseFishLenComps[2,fleet12_yrs[!fleet12_yrs %in% which(UseCatch[2,,6] == 0)],6] <- 1
UseFishLenComps[4,fleet13_yrs[!fleet13_yrs %in% which(UseCatch[4,,6] == 0)],6] <- 1
UseFishLenComps[1,(dat_4A_1$lencomp %>% filter(FltSvy == 14))$Yr,7] <- 1
UseFishLenComps[2,(dat_4A_1$lencomp %>% filter(FltSvy == 15))$Yr,7] <- 1
UseFishLenComps[4,(dat_4A_1$lencomp %>% filter(FltSvy == 16))$Yr,7] <- 1


M <- c(0.34, 0.30, 0.26, 0.22, 0.17, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14, 0.16, 0.18, 0.19, 0.20,
       0.20, 0.20, 0.19, 0.17, 0.16, 0.15, 0.14, 0.14, 0.14, 0.14, 0.13, 0.13,
       0.13)

# Get WAA
laa <- c(
  22, 35.3, 41.3, 45.7, 49.9, 53.9, 57.7, 64.2, 74.7, 85.5, 96.4, 105.3,
  112.5, 118.5, 123.4, 127.4, 130.7, 133.4, 135.6, 137.4, 138.9, 140.1, 141.1,
  141.9, 142.5, 143.1, 143.5, 144.2
)

waa <- 0.000025 * laa^2.97
mataa <- c( 0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
)

# Function to construct ALK given mean length and sd
get_al_trans_matrix = function(age_bins, len_bins, mean_length, sd) {
  # container array
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))

  for(a in 1:length(age_bins)) {
    # Use actual bin lower limits as per SS3 specification
    # Assume len_bins contains the lower limits of each bin
    bin_lower_limits = len_bins

    # Calculate cumulative probabilities at bin lower limits
    AL = pnorm(bin_lower_limits, mean_length[a], sd[a])

    # Calculate bin probabilities according to the SS3 specification
    for(l in 1:length(len_bins)) {
      if(l == 1) {
        # First bin: from -infinity to lower limit of bin 2
        age_length[a, l] = pnorm(bin_lower_limits[2], mean_length[a], sd[a])
      } else if(l == length(len_bins)) {
        # Last bin: from lower limit of last bin to +infinity
        age_length[a, l] = 1 - pnorm(bin_lower_limits[l], mean_length[a], sd[a])
      } else {
        # Middle bins: from lower limit of bin l to lower limit of bin l+1
        age_length[a, l] = pnorm(bin_lower_limits[l+1], mean_length[a], sd[a]) -
          pnorm(bin_lower_limits[l], mean_length[a], sd[a])
      }
    }
  }
  return(age_length)
}


transmat <- get_al_trans_matrix(age_bins = 1:n_ages,
                                len_bins = log(seq(10,200,5)),
                                mean_length = log(laa),
                                sd = rep(0.2, n_ages))

image(transmat)

SizeAgeTrans <- array(0, dim = c(n_regions, length(yrs), n_lens, n_ages, n_sexes))

for(r in 1:n_regions) {
  for(y in 1:length(yrs)) {
    SizeAgeTrans[r,y,,,1] <- t(transmat)
  }
}

dat_4A_1$tag_releases

# Create the tag cohort identifier (region + year combination)
df <- dat_4A_1$tag_releases %>%
  mutate(cohort = paste(reg, yr, sep = "_"))

# Reshape to 3D structure using pivot_wider and arrays
# First, let's see what unique values we have
unique_cohorts <- unique(df$cohort)
unique_ages <- 1:n_ages
unique_genders <- sort(unique(df$gender))

cat("Unique cohorts:", paste(unique_cohorts, collapse = ", "), "\n")
cat("Unique ages:", paste(unique_ages, collapse = ", "), "\n")
cat("Unique genders:", paste(unique_genders, collapse = ", "), "\n")

# Initialize the array with appropriate dimensions
cohort_array <- array(
  data = NA,
  dim = c(length(unique_cohorts), length(unique_ages), length(unique_genders)),
  dimnames = list(
    cohort = unique_cohorts,
    age = unique_ages,
    gender = unique_genders
  )
)

# Fill the array with nrel values
for(i in 1:nrow(df)) {
  cohort_idx <- which(unique_cohorts == df$cohort[i])
  age_idx <- which(unique_ages == df$age[i])
  gender_idx <- which(unique_genders == df$gender[i])

  cohort_array[cohort_idx, age_idx, gender_idx] <- df$nrel[i]
}

tag_release_indicator_df <- df %>%
  select(reg, yr) %>%
  distinct() %>%
  arrange(reg, yr) %>%
  mutate(cohort = paste(reg, yr, sep = "_")) %>%
  column_to_rownames("cohort") %>%
  as.matrix()

combined_data <- dat_4A_1$tag_recaps %>%
  # Join with release data to get original release info
  left_join(dat_4A_1$tag_releases, by = "tg", suffix = c("_recap", "_release")) %>%
  # Calculate tag years (time since release)
  mutate(
    tag_years = yr_recap - yr_release,
    age_at_recap = age + tag_years,  # Age at recapture
    cohort = paste(reg, yr_release, sep = "_")  # Tag cohort identifier
  ) %>%
  # Filter to only include tag years <= 15
  filter(tag_years <= 15 & tag_years >= 0)

# Get unique dimensions
unique_tag_years <- 0:15  # 0 to 15 tag years
unique_cohorts <- sort(unique(combined_data$cohort))
unique_regions <- sort(unique(combined_data$reg))
unique_ages <- sort(unique(combined_data$age))
unique_sexes <- sort(unique(combined_data$gender))

cat("Dimensions:\n")
cat("Tag years: 0 to 15 (16 levels)\n")
cat("Cohorts:", length(unique_cohorts), "levels\n")
cat("Regions:", length(unique_regions), "levels\n")
cat("Ages:", length(unique_ages), "levels\n")
cat("Sexes:", length(unique_sexes), "levels\n")

# Create 5-dimensional array structure
# Dimensions: [tag_years, cohorts, regions, ages, sexes]
tag_array <- array(
  data = 0,
  dim = c(
    length(unique_tag_years),
    length(unique_cohorts),
    4,
    n_ages,
    1
  )
)

# Fill the array with recapture data
for(i in 1:nrow(combined_data)) {
  ty_idx <- which(unique_tag_years == combined_data$tag_years[i])
  cohort_idx <- which(unique_cohorts == combined_data$cohort[i])
  region_idx <- which(unique_regions == combined_data$reg[i])
  age_idx <- which(unique_ages == combined_data$age[i])
  sex_idx <- which(unique_sexes == combined_data$gender[i])

  if(length(ty_idx) > 0 && length(cohort_idx) > 0 &&
     length(region_idx) > 0 && length(age_idx) > 0 &&
     length(sex_idx) > 0) {
    tag_array[ty_idx, cohort_idx, region_idx, age_idx, sex_idx] <-
      tag_array[ty_idx, cohort_idx, region_idx, age_idx, sex_idx] + combined_data$recaps[i]
  }
}

cohort_array[is.na(cohort_array)] = 0

# Setup Model -------------------------------------------------------------

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = yrs,
                            # vector of years (1 - 62)
                            ages = 1:n_ages,
                            # vector of ages (1 - 30)
                            lens = 1:n_lens,
                            # number of lengths (41 - 99)
                            n_regions = n_regions,
                            # number of regions (5)
                            n_sexes = n_sexes,
                            # number of sexes (2)
                            n_fish_fleets = n_fish_fleets,
                            # number of fishery fleet (2)
                            n_srv_fleets = 1,
                            # number of survey fleets (2)
                            verbose = TRUE
)

inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)


# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0, # not using bias ramp
                            sigmaR_switch = 1, # switch to using late sigma in year 16
                            dont_est_recdev_last = 0, # don't estimate last rec dev

                            # Model options
                            rec_dd = 'global',
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix", # fixing
                            # h_spec = 'fix',
                            # steepness = inv_steepness(0.8),
                            # InitDevs_spec = "est_shared_r",
                            # RecDevs_spec = 'est_all',
                            # initial deviations are shared across regions,
                            # but recruitment deviations are region specific
                            ln_sigmaR = log(c(1, 1)),
                            # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(97300000),
                            # starting value for global R0
                            R0_prop = array(c(0.2, 0.2, 0.2),
                                            dim = c(input_list$data$n_regions - 1))
                            # starting value for R0 proportions in multinomial logit space
)

input_list$map$ln_global_R0 = factor(NA)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = array(rep(waa, each = length(yrs) * n_regions),
                                                dim = c(n_regions, length(yrs), n_ages, n_sexes)), # weight at age
                                    MatAA = array(rep(mataa, each = length(yrs) * n_regions),
                                                  dim = c(n_regions, length(yrs), n_ages, n_sexes)), # maturity at age
                                    AgeingError = NULL,
                                    # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = SizeAgeTrans,
                                    # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    Fixed_natmort = array(rep(M, each = length(yrs) * n_regions),
                                                          dim = c(n_regions, length(yrs), n_ages, n_sexes))
                                    # values to fix natural mortality at
)


# values to fix natural mortality at
# setting up movement parameterization
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = "constant",
                                 # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_yearblk_spec = "constant", # time-invariant movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = 3
                                 # vague prior to penalize movement away from the extremes
)

# setting up tagging parameterization
input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 0, # using tagging data
                                max_tag_liberty = 16, # maximum number of years to track a cohort

                                # Data Inputs
                                tag_release_indicator = tag_release_indicator_df,
                                # tag release indicator (first col = tag region,
                                # second col = tag year),
                                # total number of rows = number of tagged cohorts
                                Tagged_Fish = cohort_array, # Released fish
                                # dimensioned by total number of tagged cohorts, (implicitly
                                # tracks the release year and region), age, and sex
                                Obs_Tag_Recap = tag_array,
                                # dimensioned by max tag liberty, tagged cohorts, regions,
                                # ages, and sexes

                                # Model options
                                Tag_LikeType = "Poisson", # Poisson
                                mixing_period = 2, # Don't fit tagging until release year + 1
                                t_tagging = 1, # tagging happens midway through the year,
                                # movement does not occur within that year
                                tag_selex = "SexSp_AllFleet", # tagging recapture selectivity
                                # is a weighted average of fishery selectivity of two fleets
                                tag_natmort = "AgeSp_SexSp", # tagging natural mortality is
                                # age and sex-specific
                                Use_TagRep_Prior = 0, # tag reporting rate priors are used
                                move_age_tag_pool = "all", # whether or
                                # not to pool tagging data when fitting (for computational cost)
                                move_sex_tag_pool = "all", # whether or not to pool
                                # sex-specific data when fitting
                                Tag_Reporting_blocks =  paste("Block_1_Year_1-terminal_Region_",
                                                              c(1:input_list$data$n_regions), sep = ''),
                                Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                Tag_Shed_spec = "fix", # fixing chronic shedding
                                TagRep_spec = "fix", # tag reporting rates are
                                # not region specific
                                # Specify starting values or fixing values
                                ln_Init_Tag_Mort = log(1e-10), # fixing initial tag mortality
                                ln_Tag_Shed = log(1e-10),  # fixing tag shedding
                                ln_tag_theta = log(0.5),
                                # starting value for tagging overdispersion
                                Tag_Reporting_Pars = array(log(0.9 / (1-0.9)),
                                                           dim = c(input_list$data$n_regions, 1))
                                # starting values for tag reporting pars

)

# input_list <- Setup_Mod_Tagging(input_list, 0)

# setting up catch data
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = ObsCatch,
                                    Catch_Type = array(1, dim = c(length(yrs), n_fish_fleets)),
                                    catch_units = array("abd", dim = c(length(yrs), n_fish_fleets)),
                                    UseCatch = UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    ln_sigmaC =
                                      array(log(0.02), dim = c(input_list$data$n_regions,
                                                               length(input_list$data$years),
                                                               input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-5, dim = c(input_list$data$n_regions,
                                                                  input_list$data$n_fish_fleets))
                                    # some starting values for fishing mortality
)

# Fishery Indices and Compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = ObsFishIdx,
                                          ObsFishIdx_SE = array(0.2, dim(ObsFishIdx)),
                                          UseFishIdx =  UseFishIdx,
                                          ObsFishAgeComps = array(0, dim = c(n_regions, length(yrs), n_ages, n_sexes, n_fish_fleets)),
                                          UseFishAgeComps = array(0, dim(ObsFishIdx)),
                                          ObsFishLenComps = ObsFishLenComps,
                                          UseFishLenComps = UseFishLenComps,
                                          ISS_FishAgeComps = NULL,
                                          ISS_FishLenComps = array(5, dim = c(n_regions, length(yrs), n_sexes, n_fish_fleets)),

                                          # Model options
                                          fish_idx_type = rep("biom", n_fish_fleets),
                                          # fishery indices not used
                                          FishAgeComps_LikeType =
                                            rep("none", n_fish_fleets),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType =
                                            rep("Multinomial", n_fish_fleets),
                                          # length comp likelihoods for fishery fleet 1 and 2
                                          FishAgeComps_Type =
                                            paste("none_Year_1-terminal_Fleet", 1:n_fish_fleets, sep = '_'),
                                          # age comp structure for fishery fleet 1 and 2
                                          FishLenComps_Type =
                                            paste("spltRspltS_Year_1-terminal_Fleet", 1:n_fish_fleets, sep = '_'),
                                          # length comp structure for fishery fleet 1 and 2
                                          FishAge_comp_agg_type = rep(NA, n_fish_fleets),
                                          # ADMB aggregation quirks, ideally get rid of this
                                          FishLen_comp_agg_type = rep(0, n_fish_fleets)
                                          # ADMB aggregation quirks, ideally get rid of this
)

# Survey Indices and Compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,

                                         # data inputs
                                         ObsSrvIdx = array(0, dim = c(n_regions, length(yrs), input_list$data$n_srv_fleets)),
                                         ObsSrvIdx_SE = array(0, dim = c(n_regions, length(yrs), input_list$data$n_srv_fleets)),
                                         UseSrvIdx = array(0, dim = c(n_regions, length(yrs), input_list$data$n_srv_fleets)),
                                         ObsSrvAgeComps = array(0, dim = c(n_regions, length(yrs), n_ages, n_sexes, input_list$data$n_srv_fleets)),
                                         UseSrvAgeComps = array(0, dim = c(n_regions, length(yrs), input_list$data$n_srv_fleets)),
                                         ObsSrvLenComps = array(0, dim = c(n_regions, length(yrs), n_lens, n_sexes, input_list$data$n_srv_fleets)),
                                         UseSrvLenComps = array(0, dim = c(n_regions, length(yrs), input_list$data$n_srv_fleets)),
                                         ISS_SrvAgeComps = NULL,
                                         ISS_SrvLenComps = NULL,

                                         # Model options
                                         srv_idx_type = rep("none", input_list$data$n_srv_fleets),
                                         SrvAgeComps_LikeType = "none",
                                         SrvLenComps_LikeType = "none",
                                         SrvAgeComps_Type = "none_Year_1-terminal_Fleet_1",
                                         SrvLenComps_Type = "none_Year_1-terminal_Fleet_1",
                                         SrvAge_comp_agg_type = NA,
                                         SrvLen_comp_agg_type = NA
)


# Fishery Selectivity and Catchability
input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = paste('none_Fleet', 1:n_fish_fleets, sep = '_'),
                                      # fishery selectivity, whether continuous time-varying

                                      # fishery selectivity blocks
                                      fish_sel_blocks =
                                        paste('none_Fleet', 1:n_fish_fleets, sep = '_'),
                                      # no blocks for trawl fishery

                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1",
                                          "logist1_Fleet_2",
                                          "logist1_Fleet_3",
                                          "logist1_Fleet_4",
                                          "logist1_Fleet_5",
                                          "logist1_Fleet_6",
                                          "logist1_Fleet_7"),

                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        paste('none_Fleet', 1:n_fish_fleets, sep = '_'),
                                      # no blocks since q is not estimated

                                      # whether to estimate all fixed effects
                                      # for fishery selectivity and later modify
                                      # to fix and share parameters
                                      fish_fixed_sel_pars_spec =
                                        rep("est_shared_r", n_fish_fleets),

                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        rep("est_shared_r", n_fish_fleets)
                                      # fix fishery q since not used
)

input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel = "none_Fleet_1",

                                     # survey selectivity blocks
                                     srv_sel_blocks = "none_Fleet_1", # no blocks for jp and domestic survey

                                     # survey selectivity form
                                     srv_sel_model = "logist1_Fleet_1",

                                     # survey catchability blocks
                                     srv_q_blocks = "none_Fleet_1",

                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec = "fix",

                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     # spatially-invariant q
                                     srv_q_spec = "fix"
)

# set up model weighting stuff
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  Wt_Catch = 1,
                                  Wt_FishIdx = 1,
                                  Wt_SrvIdx = 1,
                                  Wt_Rec = 1,
                                  Wt_F = 1,
                                  Wt_Tagging = 1,
                                  # Composition model weighting
                                  Wt_FishAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map

data$ISS_FishLenComps[,,,] <- 2
# parameters$ln_fish_fixed_sel_pars[] <- log(10)

# make AD model function
obj <- RTMB::MakeADFun(SPoRC:::cmb(SPoRC:::SPoRC_rtmb, data),
                       parameters = parameters,
                       map = mapping,
                       random = NULL, silent = F)

# Now, optimize the function
optim <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

# newton steps
try_improve <- tryCatch(expr =
                          for(i in 1:4) {
                            g = as.numeric(obj$gr(optim$par))
                            h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                            optim$par = optim$par - solve(h,g)
                            optim$objective = obj$fn(optim$par)
                          }
                        , error = function(e){e}, warning = function(w){w})

# save report
obj$rep <- obj$report(obj$env$last.par.best)
plot(obj$rep$NAA[1,1,,1], type = 'l')

# sdrep <- sdreport(obj)

load("/Users/matthewcheng/Downloads/IOTC_SS_final.Rdata")
plot(IOTC_SS_spat_bio[1,,1], type = 'l')
lines(obj$rep$SSB[1,], type = 'l')

plot(IOTC_SS_spat_bio[2,,1], type = 'l')
lines(obj$rep$SSB[2,], type = 'l')

plot(IOTC_SS_spat_bio[3,,1], type = 'l')
lines(obj$rep$SSB[3,], type = 'l')

plot(IOTC_SS_spat_bio[4,,1], type = 'l')
lines(obj$rep$SSB[4,], type = 'l')

get_catch_fits_plot(list(data), list(obj$rep), '1')
get_idx_fits_plot(list(data), list(obj$rep), '1')

plot(obj$rep$SSB[1,], type = 'l')
lines(obj$rep$SSB[2,])
lines(obj$rep$SSB[3,])
lines(obj$rep$SSB[4,])

plot(obj$rep$Rec[1,], type = 'l')
plot(obj$rep$Rec[2,], type = 'l')
plot(obj$rep$Rec[3,], type = 'l')
plot(obj$rep$Rec[4,], type = 'l')

plot(obj$rep$PredCatch[1,,1])
lines(data$ObsCatch[1,,1])

obj$rep$jnLL
sum(obj$rep$Catch_nLL)
sum(obj$rep$Fmort_nLL)
sum(obj$rep$FishIdx_nLL)
sum(obj$rep$FishLenComps_nLL)

sum(obj$rep$FishLenComps_nLL[2,235,,6])

plot(data$ObsFishIdx[1,,3])
lines(obj$rep$PredFishIdx[1,,3])

plot(data$ObsFishIdx[2,,3])
lines(obj$rep$PredFishIdx[2,,3])

plot(data$ObsFishIdx[3,,3])
lines(obj$rep$PredFishIdx[3,,3])

plot(data$ObsFishIdx[4,,3])
lines(obj$rep$PredFishIdx[4,,3])

plot(data$ObsFishIdx[1,,3])
lines(obj$rep$PredFishIdx[1,,3])


plot(obj$rep$Fmort[1,,1], type = 'l')

