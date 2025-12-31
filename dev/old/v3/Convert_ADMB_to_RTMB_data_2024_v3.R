# Purpose: To convert ADMB data to RTMB format for the 2024 sablefish assessment model

library(here)
library(R2admb)
library(tidyverse)

# Read in data
tem_dat <- dget(here("dev", '2024 Base (23.5)_final model_v3', 'tem.rdat'))
ageing_dat <- dget(here("dev",'2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error
tem_admb_dat <- readLines(here("dev","2024 Base (23.5)_final model_v3", "tem_2024_na_wh.dat"))
tem_rep <- readLines(here("dev",'2024 Base (23.5)_final model_v3', 'sable.rep'))
tem_par <- read_pars(here("dev",'2024 Base (23.5)_final model_v3', 'tem'))

# Dimensions
ages <- 1:30 # ages
years <- rownames(tem_dat$t.series) # years
n_regions <- 1
lens <- seq(41, 99, 2)
n_sexes <- 2
n_fish_fleets <- 2
n_srv_fleets <- 3

# Biologicals -------------------------------------------------------------
# Setup biological arrays
WAA <- array(NA, dim = c(n_regions, length(years), length(ages), n_sexes))
MatAA <- array(NA, dim = c(n_regions, length(years), length(ages), n_sexes))
WAA[1,,,1] <- matrix(tem_dat$growthmat$wt.f.block1, nrow = length(years), ncol = length(ages), byrow = TRUE)
WAA[1,,,2] <- matrix(tem_dat$growthmat$wt.m.block1, nrow = length(years), ncol = length(ages), byrow = TRUE)
MatAA[1,1:length(years),,1] <- matrix(tem_dat$growthmat$mage.block2, nrow = length(years), ncol = length(ages), byrow = TRUE)
MatAA[1,1:length(years),,2] <- matrix(tem_dat$growthmat$mage.block2, nrow = length(years), ncol = length(ages), byrow = TRUE)

# Size Transition Matrix (Growth)
SizeAgeTrans <- array(0, dim = c(n_regions, length(years), length(lens), length(ages), n_sexes)) # size age transition matrix
# first time block
SizeAgeTrans[1,1:35,,,1] <- aperm(replicate(length(1:35), tem_dat$sizeage.f.block1), perm = c(3,2,1)) # female
SizeAgeTrans[1,1:35,,,2] <- aperm(replicate(length(1:35), tem_dat$sizeage.m.block1), perm = c(3,2,1)) # male
# Size age transition, second time block
SizeAgeTrans[1,36:length(years),,,1] <- aperm(replicate(length(36:length(years)), tem_dat$sizeage.f.block2), perm = c(3,2,1)) # female
SizeAgeTrans[1,36:length(years),,,2] <- aperm(replicate(length(36:length(years)), tem_dat$sizeage.m.block2), perm = c(3,2,1)) # male


# Catch -------------------------------------------------------------------
ObsCatch <- array(NA, c(n_regions, length(1:length(years)), n_fish_fleets))
ObsCatch[1,,1] <- as.numeric(strsplit(tem_admb_dat[46], split = " ")[[1]])[1:length(years)] # fixed gear catches
ObsCatch[1,,2] <- as.numeric(strsplit(tem_admb_dat[48], split = " ")[[1]])[1:length(years)] # trawl gear catches
ObsCatch[ObsCatch == 0] <- NA # set 0 catches to NA so we aren't fitting
UseCatch <- array(1, dim = c(n_regions, length(years), n_fish_fleets))
UseCatch[1,1:3,2] <- 0 # don't use first 3 years of trawl catches (0)


# Fishery Indices ---------------------------------------------------------
ObsFishIdx <- array(NA, c(n_regions, length(1:length(years)), n_fish_fleets))
ObsFishIdx_SE <- array(NA, c(n_regions, length(1:length(years)), n_fish_fleets))
UseFishIdx <- array(1, c(n_regions, length(1:length(years)), n_fish_fleets))
colnames(ObsFishIdx) <- years # define row years
colnames(ObsFishIdx_SE) <- years # define row years

ObsFishIdx[1,colnames(ObsFishIdx) %in% rownames(tem_dat$obssrv6),1] <- as.numeric(strsplit(tem_admb_dat[144], split = " ")[[1]])
ObsFishIdx_SE[1,colnames(ObsFishIdx_SE) %in% rownames(tem_dat$obssrv6),1] <- as.numeric(strsplit(tem_admb_dat[146], split = " ")[[1]])

# Domestic LL Fishery (after 1995)
ObsFishIdx[1,colnames(ObsFishIdx) %in% rownames(tem_dat$obssrv5),1] <- as.numeric(strsplit(tem_admb_dat[127], split = " ")[[1]])
ObsFishIdx_SE[1,colnames(ObsFishIdx_SE) %in% rownames(tem_dat$obssrv5),1] <- as.numeric(strsplit(tem_admb_dat[129], split = " ")[[1]])
UseFishIdx[is.na(ObsFishIdx)] <- 0 # don't fit if missing data


# Fishery Age Comps -------------------------------------------------------
# Note that NA is in trawl fishery slot so it doesn't fit
ObsFishAgeComps <- array(NA, dim = c(n_regions, length(years), length(ages), n_sexes, n_fish_fleets))
colnames(ObsFishAgeComps) <- years # define row years

# Fishery Age Comps (Sex Specific)
ObsFishAgeComps[1,colnames(ObsFishAgeComps) %in% rownames(tem_dat$oac.fish1),,1,1] <- tem_dat$oac.fish1 # Observed fixed gear fishery age comps (aggregated)
UseFishAgeComps <- array(0, dim = c(n_regions, length(years), n_fish_fleets))
colnames(UseFishAgeComps) <- years # define row years
UseFishAgeComps[1,colnames(UseFishAgeComps) %in% rownames(tem_dat$oac.fish1),1] <- 1 # only fit if have age comp data

# Data weighting for fishery age compositions
ISS_FishAgeComps <- array(0, dim = c(n_regions, length(years), n_sexes, n_fish_fleets))
colnames(ISS_FishAgeComps) <- years # define row years
ISS_FishAgeComps[1,colnames(ISS_FishAgeComps) %in% rownames(tem_dat$oac.fish1),1,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery ages (aggregated)


# Fishery Length Comps ----------------------------------------------------
ObsFishLenComps <- array(NA, dim = c(n_regions, length(years), length(lens), n_sexes, n_fish_fleets))
colnames(ObsFishLenComps) <- years # define row years

# observed fixed and trawl gear fishery length comps (joint)
# observed fixed gear fishery length comps
ObsFishLenComps[1,colnames(ObsFishLenComps) %in% rownames(tem_dat$olc.fish1.f),,1,1] <- tem_dat$olc.fish1.f # females
ObsFishLenComps[1,colnames(ObsFishLenComps) %in% rownames(tem_dat$olc.fish1.m),,2,1] <- tem_dat$olc.fish1.m # males

# observed trawl gear fishery length comps
ObsFishLenComps[1,colnames(ObsFishLenComps) %in% rownames(tem_dat$olc.fish3.f),,1,2] <- tem_dat$olc.fish3.f # females
ObsFishLenComps[1,colnames(ObsFishLenComps) %in% rownames(tem_dat$olc.fish3.m),,2,2] <- tem_dat$olc.fish3.m # males

# Use fishery length comp controls
UseFishLenComps <- array(0, dim = c(n_regions, length(years), n_fish_fleets))
colnames(UseFishLenComps) <- years # define row years
UseFishLenComps[1,colnames(UseFishLenComps) %in% rownames(tem_dat$olc.fish1.f),1] <- 1 # only fit if have len comp data
UseFishLenComps[1,colnames(UseFishLenComps) %in% rownames(tem_dat$olc.fish3.m),2] <- 1 # only fit if have len comp data

# Data weighting for fishery length compositions
ISS_FishLenComps <- array(0, dim = c(n_regions, length(years), n_sexes, n_fish_fleets))
colnames(ISS_FishLenComps) <- years # define row years

# split iss
ISS_FishLenComps[1,colnames(ISS_FishLenComps) %in% rownames(tem_dat$olc.fish1.f),1,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery females
ISS_FishLenComps[1,colnames(ISS_FishLenComps) %in% rownames(tem_dat$olc.fish1.m),2,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery males
ISS_FishLenComps[1,colnames(ISS_FishLenComps) %in% rownames(tem_dat$olc.fish3.f),1,2] <- 20 # Assuming constant ISS of 20 for trawl gear fishery females
ISS_FishLenComps[1,colnames(ISS_FishLenComps) %in% rownames(tem_dat$olc.fish3.m),2,2] <- 20 # Assuming constant ISS of 20 for trawl gear fishery males


# Survey Indices ----------------------------------------------------------
ObsSrvIdx <- array(NA, dim = c(n_regions, length(1:length(years)), n_srv_fleets))
ObsSrvIdx_SE <- array(NA, dim = c(n_regions, length(1:length(years)), n_srv_fleets))
UseSrvIdx <- array(1, dim = c(n_regions, length(1:length(years)), n_srv_fleets))
colnames(ObsSrvIdx) <- years # define row years
colnames(ObsSrvIdx_SE) <- years # define row years

# Domestic LL Survey post 1995
ObsSrvIdx[1,colnames(ObsSrvIdx) %in% rownames(tem_dat$obssrv3),1] <- as.numeric(strsplit(tem_admb_dat[93], split = " ")[[1]])
ObsSrvIdx_SE[1,colnames(ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv3),1] <- as.numeric(strsplit(tem_admb_dat[95], split = " ")[[1]])

# Domestic Trawl Survey
ObsSrvIdx[1,colnames(ObsSrvIdx) %in% rownames(tem_dat$obssrv7),2] <- as.numeric(strsplit(tem_admb_dat[161], split = " ")[[1]])
ObsSrvIdx_SE[1,colnames(ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv7),2] <- as.numeric(strsplit(tem_admb_dat[163], split = " ")[[1]])

# Coop LL Survey pre 1995
ObsSrvIdx[1,colnames(ObsSrvIdx) %in% rownames(tem_dat$obssrv4),3] <- as.numeric(strsplit(tem_admb_dat[110], split = " ")[[1]])
ObsSrvIdx_SE[1,colnames(ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv4),3] <- as.numeric(strsplit(tem_admb_dat[112], split = " ")[[1]])
UseSrvIdx[is.na(ObsSrvIdx)] <- 0 # don't fit if missing data


# Survey Age Comps --------------------------------------------------------
ObsSrvAgeComps <- array(NA, dim = c(n_regions, length(years), length(ages), n_sexes, n_srv_fleets))
colnames(ObsSrvAgeComps) <- years # define row years

# agg ages
ObsSrvAgeComps[1,colnames(ObsSrvAgeComps) %in% rownames(tem_dat$oac.srv1),,1,1] <- tem_dat$oac.srv1 # Observed domestic ll survey age comps (aggregated)
ObsSrvAgeComps[1,colnames(ObsSrvAgeComps) %in% rownames(tem_dat$oac.srv2),,1,3] <- tem_dat$oac.srv2 # Observed coop jp ll survey age comps (aggregated)

UseSrvAgeComps <- array(0, dim = c(n_regions, length(1:length(years)), n_srv_fleets))
colnames(UseSrvAgeComps) <- years # define row years
UseSrvAgeComps[1,colnames(UseSrvAgeComps) %in% rownames(tem_dat$oac.srv1),1] <- 1 # only fit if have age comp data
UseSrvAgeComps[1,colnames(UseSrvAgeComps) %in% rownames(tem_dat$oac.srv2),3] <- 1 # only fit if have age comp data

# Data weighting for fishery age compositions
ISS_SrvAgeComps <- array(0, dim = c(n_regions, length(years), n_sexes, n_srv_fleets))
colnames(ISS_SrvAgeComps) <- years # define row years

ISS_SrvAgeComps[1,colnames(ISS_SrvAgeComps) %in% rownames(tem_dat$oac.srv1),1,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey ages (aggregated)
ISS_SrvAgeComps[1,colnames(ISS_SrvAgeComps) %in% rownames(tem_dat$oac.srv2),1,3] <- 20 # Assuming constant ISS of 20 for coop jp survey ages (aggregated)


# Survey Length Comps -----------------------------------------------------
ObsSrvLenComps <- array(NA, dim = c(n_regions,length(years), length(lens), n_sexes, n_srv_fleets))
colnames(ObsSrvLenComps) <- years # define row years

# observed domestic survey ll length comps
ObsSrvLenComps[1,colnames(ObsSrvLenComps) %in% rownames(tem_dat$olc.srv1.f),,1,1] <- tem_dat$olc.srv1.f # females
ObsSrvLenComps[1,colnames(ObsSrvLenComps) %in% rownames(tem_dat$olc.srv1.m),,2,1] <- tem_dat$olc.srv1.m # males

# observed domestic trawl survey length comps
ObsSrvLenComps[1,colnames(ObsSrvLenComps) %in% rownames(tem_dat$olc.srv7.f),,1,2] <- tem_dat$olc.srv7.f # females
ObsSrvLenComps[1,colnames(ObsSrvLenComps) %in% rownames(tem_dat$olc.srv7.m),,2,2] <- tem_dat$olc.srv7.m # males

# observed coop jp ll survey length comps
srv_trawl_iss <- as.numeric(strsplit(tem_admb_dat[586], split = " ")[[1]]) # get ISS from trawl survey
ObsSrvLenComps[1,colnames(ObsSrvLenComps) %in% rownames(tem_dat$olc.srv2.f),,1,3] <- tem_dat$olc.srv2.f # females
ObsSrvLenComps[1,colnames(ObsSrvLenComps) %in% rownames(tem_dat$olc.srv2.m),,2,3] <- tem_dat$olc.srv2.m # males

# Use survey length comp controls
UseSrvLenComps <- array(0, dim = c(n_regions, length(years), n_srv_fleets))
colnames(UseSrvLenComps) <- years # define row years
UseSrvLenComps[1,colnames(UseSrvLenComps) %in% rownames(tem_dat$olc.srv1.f),1] <- 1 # only fit if have len comp data
UseSrvLenComps[1,colnames(UseSrvLenComps) %in% rownames(tem_dat$olc.srv7.f),2] <- 1 # only fit if have len comp data
UseSrvLenComps[1,colnames(UseSrvLenComps) %in% rownames(tem_dat$olc.srv2.f),3] <- 1 # only fit if have len comp data

# Data weighting for survey length compositions
ISS_SrvLenComps <- array(0, dim = c(n_regions,length(years), n_sexes, n_srv_fleets))
colnames(ISS_SrvLenComps) <- years # define row years
ISS_SrvLenComps[1,colnames(ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv1.f),1,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey females
ISS_SrvLenComps[1,colnames(ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv1.m),2,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey males
ISS_SrvLenComps[1,colnames(ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv7.f),1,2] <- srv_trawl_iss # Assuming constant ISS of 20 for domestic trawl survey females
ISS_SrvLenComps[1,colnames(ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv7.m),2,2] <- srv_trawl_iss # Assuming constant ISS of 20 for domestic trawl survey males
ISS_SrvLenComps[1,colnames(ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv2.f),1,3] <- 20 # Assuming constant ISS of 20 for coop jp ll survey females
ISS_SrvLenComps[1,colnames(ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv2.m),2,3] <- 20 # Assuming constant ISS of 20 for coop jp ll survey males

# output into list
RTMB_data <- list(WAA = WAA,
                  MatAA = MatAA,
                  SizeAgeTrans = SizeAgeTrans,
                  ObsCatch = ObsCatch,
                  UseCatch = UseCatch,
                  ObsFishIdx = ObsFishIdx,
                  ObsFishIdx_SE = ObsFishIdx_SE,
                  UseFishIdx = UseFishIdx,
                  ObsFishAgeComps = ObsFishAgeComps,
                  UseFishAgeComps = UseFishAgeComps,
                  ISS_FishAgeComps = ISS_FishAgeComps,
                  ObsFishLenComps = ObsFishLenComps,
                  UseFishLenComps = UseFishLenComps,
                  ISS_FishLenComps = ISS_FishLenComps,
                  ObsSrvIdx = ObsSrvIdx,
                  ObsSrvIdx_SE = ObsSrvIdx_SE,
                  UseSrvIdx = UseSrvIdx,
                  ObsSrvAgeComps = ObsSrvAgeComps,
                  UseSrvAgeComps = UseSrvAgeComps,
                  ISS_SrvAgeComps = ISS_SrvAgeComps,
                  ObsSrvLenComps = ObsSrvLenComps,
                  UseSrvLenComps = UseSrvLenComps,
                  ISS_SrvLenComps = ISS_SrvLenComps
                  )


write_rds(RTMB_data, file = here("dev", '2024 Base (23.5)_final model', 'RTMB.RDS'))
