

rep <- three_rep <- readRDS("/Users/matthewcheng/Desktop/PostDoc/Spatial Assessments and Sablefish/SPoRC/dev/dev_output/3_Region_Model_Sablefish/rep.RDS")

R0 <- 10
sexratio <- 1
init_iter <- 1e4
n_ages <- 30
n_sexes <- 1
n_regions <- 3
do_recruits_move <- 0
Movement <- three_rep$Movement
# Movement <- array(diag(1,n_regions), c(n_regions,n_regions,1,30,2))
Init_NAA_next_year <- Init_NAA <- array(0, dim = c(n_regions, n_ages, n_sexes))
natmort <- 0.1 * exp(-1 * 1:n_ages)
# but this does
natmort[1:5] <- 0.1
natmort[6:14] <- 0.2
natmort[27:29] <- 0.3
natmort[30] <- 0.3
# natmort[] <- seq(0.1, 0.05, length.out = n_ages)

# Set up initial equilibrium age structure
for(r in 1:n_regions) {
  for(s in 1:n_sexes) {
    tmp_cumsum_Z <- cumsum(natmort[1:(n_ages-1)]) # cumulative sum of total mortality
    Init_NAA[r,,s] <- c(R0 * sexratio, R0 * sexratio * exp(-tmp_cumsum_Z)) # exponential mortality model
  } # end s loop
} # end r loop

# Apply annual cycle and iterate to equilibrium
for(i in 1:30) {
  for(s in 1:n_sexes) {
    Init_NAA_next_year[,1,s] <- R0 * sexratio # recruitment
    # movement
    if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] <- t(Init_NAA[,a,s]) %*% Movement[,,1,a,1] # recruits dont move
    if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] <- t(Init_NAA[,a,s]) %*% Movement[,,1,a,1] # recruits move
    for(r in 1:n_regions) {
      # ageing and mortality
      Init_NAA_next_year[r,2:n_ages,s] <- Init_NAA[r,1:(n_ages-1),s] * exp(-(natmort[1:(n_ages-1)] ))
      # accumulate plus group
      Init_NAA_next_year[r,n_ages,s] <- Init_NAA_next_year[r,n_ages,s]  + Init_NAA[r,n_ages,s] * exp(-(natmort[n_ages] ))
    }
    Init_NAA <- Init_NAA_next_year # iterate to next cycle
  } # end s loop

  if(sum(Init_NAA_next_year[,29,] - Init_NAA[,29,]) != 0) print(Init_NAA_next_year[,29,] - Init_NAA[,29,])
} # end i loop

lines(Init_NAA_next_year[1,,])

# movement matrices
M29 <- Movement[,,1,29,1]
M30 <- Movement[,,1,30,1]

# age-specific survival from your natmort vector:
s29 <- exp(-natmort[29])
s30 <- exp(-natmort[30])

# prev-age vector (before you move/apply survival in the analytic expr)
N29 <- Init_NAA[,29,1]

# Source (new entrants): move age-29 fish then apply age-29 survival
# Using movement orientation consistent with t(N) %*% Movement in the loop -> use t(M29) %*% N29
source <- as.vector(t(M29) %*% N29) * s29
# (equivalently: source <- ( (t(N29) %*% M29) ) * s29 )

# Transition operator for plus-group survivors:
# survivors are (M30^T %*% N30) * s30, so the LHS operator is diag(s30) %*% t(M30)
T_mat_correct <- s30 * t(M30)  # NOT diag(s30)Retry
I_mat <- diag(n_regions)

# solve linear system
N30_equil <- solve(I_mat - T_mat_correct, source)

# compare
cbind(iterative = round(Init_NAA[,30,1], 8), analytic = round(as.numeric(N30_equil), 8), diff = round(Init_NAA[,30,1] - as.numeric(N30_equil), 8))


eigen(T_mat_correct)$values



# projection initial abundance forward
Rec_trans_prop <- rep$Rec_trans_prop
sexratio <- rep$sexratio
R0 <- rep$R0
Movement <- rep$Movement
init_F <- rep$init_F
fish_sel <- rep$fish_sel
natmort <- rep$natmort
n_regions <- 3
n_ages <- 30
n_sexes <- 2
Init_NAA_next_year <- Init_NAA <- array(0, dim = c(n_regions, n_ages, n_sexes))


for(i in 1:c(n_ages * 5)) {
for(s in 1:n_sexes) {
  Init_NAA[,1,s] = R0 * sexratio[,1,s] * Rec_trans_prop # initialize recruitment
  # movement
  if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,1,a,s] # recruits don't move
  if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,1,a,s] # recruits move
  for(r in 1:n_regions) {
    # ageing and mortality
    Init_NAA[r,2:n_ages,s] = Init_NAA[r,1:(n_ages-1),s] * exp(-(natmort[r,1,1:(n_ages-1),s] + (init_F * fish_sel[r,1,1:(n_ages-1),s,1])))
    # accumulate plus group
    Init_NAA[r,n_ages,s] = (Init_NAA[r,n_ages,s]) + (Init_NAA[r,n_ages,s] * exp(-(natmort[r,1,n_ages,s] + (init_F * fish_sel[r,1,n_ages,s,1]))))
  } # end r loop
} # end s loop
} # end a loop

lines(Init_NAA[1,,1])

# [1] 5.51052012 4.96183193 3.98721046 3.25425304 2.69550992
# [6] 2.26320242 1.92342411 1.27216032 0.89233339 0.66408577
# [11] 0.52120628 0.42704611 0.36125424 0.31247109 0.27430349
# [16] 0.24310674 0.22232675 0.20218026 0.18320543 0.16563514
# [21] 0.14953267 0.13486995 0.12157221 0.10954338 0.09868026
# [26] 0.08888021 0.08004519 0.07208363 0.06491118 0.11108159

# [1] 5.51052012 4.96183193 3.98721046 3.25425304 2.69550992
# [6] 2.26320242 1.92342411 1.27216032 0.89233339 0.66408577
# [11] 0.52120628 0.42704611 0.36125424 0.31247109 0.27430349
# [16] 0.24310674 0.22232675 0.20218026 0.18320543 0.16563514
# [21] 0.14953267 0.13486995 0.12157221 0.10954338 0.09868026
# [26] 0.08888021 0.08004519 0.07208363 0.06491118 0.11108159
