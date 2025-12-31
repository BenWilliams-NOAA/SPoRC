
single_region_SPR_int = function(log_F_x, data) {
  n_ages = dim(data$fish_sel)[1] # number of ages
  # exponentitate reference points to "estimate"
  F_x = exp(log_F_x)
  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_ages)) # 2 slots in rows, for unfished, and fished at F_x
  # Set up the initial recruits
  Nspr[,1] = 1
  # Loop through and decrement recruit
  for(j in 2:(n_ages - 1)) {
    Nspr[1,j] = Nspr[1,j-1] * exp(-1 * data$natmort[j-1]) # unfished
    Nspr[2,j] = Nspr[2,j-1] * exp(-1 * (data$natmort[j-1] + sum(data$F_fract_flt * F_x * data$fish_sel[j-1,]))) # fished
  }
  # Accumulate plus group
  Nspr[1,n_ages] = Nspr[1,n_ages-1] * exp(-1 * data$natmort[n_ages-1])/(1-exp(-1*data$natmort[n_ages])) # unfished
  Nspr[2,n_ages] = Nspr[2,n_ages-1] * exp(-1 * (data$natmort[n_ages-1] + sum(data$F_fract_flt * F_x * data$fish_sel[n_ages-1,])))/
    (1 - exp(-1 * (data$natmort[n_ages] + sum(data$F_fract_flt * F_x * data$fish_sel[n_ages,])))) # fished
  # Convert numbers at age to spawning biomass at age (data$t_spwn accounts for mortality up until spawning)
  for(j in 1:n_ages) {
    SB_age[1,j] = Nspr[1,j] * data$WAA[j] * data$MatAA[j] * exp(-data$t_spwn * data$natmort[j]) # unfished
    SB_age[2,j] = Nspr[2,j] * data$WAA[j] * data$MatAA[j] * exp(-data$t_spwn * (data$natmort[j] + sum(data$F_fract_flt * F_x * data$fish_sel[j,]))) # fished
  }
  # Get spawning biomass per recruit to get spawning potential ratio
  SB0 = sum(SB_age[1,])
  SB_F_x = sum(SB_age[2,])
  SPR = SB_F_x / SB0
  # compute objective function to get F_x
  sprpen = 100 * (SPR - data$SPR_x)^2
  return(sprpen)
}

data_fx = list()
# Extract out relevant elements
n_ages <- length(ages) # number of ages
n_years <- length(years) # number of years
data_fx$t_spwn <- t_spwn # specified mortality time up until spawning
# fishing mortality fraction
data_fx$F_fract_flt <- Fmort[1,n_years,] / sum(Fmort[1,n_years,]) # get fleet F fraction to derive population level selectivity
# fishery selectivity
data_fx$fish_sel <- array(fish_sel[,1,,1,,drop = FALSE], dim = c(n_ages, n_fish_fleets)) # get female selectivity for all fleets
# natural mortality
data_fx$natmort <- as.vector(natmort[,1,,1,drop = FALSE]) # get female natural mortality
# weight at age
data_fx$WAA <- WAA[,1,,1,drop = FALSE] # weight at age for females
# maturity at age
data_fx$MatAA <- MatAA[,1,,1,drop = FALSE] # maturity at age for females
data_fx$SPR_x = 0.4

log_fx = log(0.1) # init guess for log(F_x)
# Newton-Raphson optimization loop
# Newton-Raphson optimization loop with numerical gradients/Hessians
h = 1e-3 # step size for numerical differentiation

for(i in 1:4) {
  objective = function(x) single_region_SPR_int(x, data_fx)
  grad = (objective(log_fx + h) - objective(log_fx - h)) / (2 * h)
  hess = (objective(log_fx + h) - 2 * objective(log_fx) + objective(log_fx - h)) / (h^2)
  log_fx = log_fx - as.numeric(grad / hess)
}

RTMB::ADREPORT(log_fx)
RTMB::REPORT(log_fx)
