if(n_regions > 1) {
  # Apply continuous movement + mortality
  for(s in 1:n_sexes) {

    # Age-1 (recruits) - age forward with mortality only, no movement
    NAA[,y+1,2,s] = NAA[,y,1,s] * exp(-ZAA[,y,1,s])
    NAA0[,y+1,2,s] = NAA0[,y,1,s] * exp(-natmort[,y,1,s])

    # Ages 2+ with continuous movement + mortality
    for(a in 2:(n_ages-1)) {
      # Fished
      Q_fished = Mrate[,,y,a,s]
      diag(Q_fished) = diag(Q_fished) - ZAA[,y,a,s]
      NAA[,y+1,a+1,s] = as.vector(Matrix::expm(Q_fished) %*% NAA[,y,a,s])

      # Unfished
      Q_unfished = Mrate[,,y,a,s]
      diag(Q_unfished) = diag(Q_unfished) - natmort[,y,a,s]
      NAA0[,y+1,a+1,s] = as.vector(Matrix::expm(Q_unfished) %*% NAA0[,y,a,s])
    }

    # Plus group
    Q_fished_plus = Mrate[,,y,n_ages,s]
    diag(Q_fished_plus) = diag(Q_fished_plus) - ZAA[,y,n_ages,s]
    NAA[,y+1,n_ages,s] = NAA[,y+1,n_ages,s] + as.vector(Matrix::expm(Q_fished_plus) %*% NAA[,y,n_ages,s])

    Q_unfished_plus = Mrate[,,y,n_ages,s]
    diag(Q_unfished_plus) = diag(Q_unfished_plus) - natmort[,y,n_ages,s]
    NAA0[,y+1,n_ages,s] = NAA0[,y+1,n_ages,s] + as.vector(Matrix::expm(Q_unfished_plus) %*% NAA0[,y,n_ages,s])
  } # end s loop
} else {
  # No movement - just mortality and ageing
  NAA[,y+1,2:n_ages,] = NAA[,y,1:(n_ages-1),] * exp(-ZAA[,y,1:(n_ages-1),])
  NAA[,y+1,n_ages,] = NAA[,y+1,n_ages,] + NAA[,y,n_ages,] * exp(-ZAA[,y,n_ages,])

  NAA0[,y+1,2:n_ages,] = NAA0[,y,1:(n_ages-1),] * exp(-natmort[,y,1:(n_ages-1),])
  NAA0[,y+1,n_ages,] = NAA0[,y+1,n_ages,] + NAA0[,y,n_ages,] * exp(-natmort[,y,n_ages,])
}

### Compute Biomass Quantities ----------------------------------------------
Total_Biom[,y] = apply(NAA[,y,,,drop = FALSE] * WAA[,y,,,drop = FALSE], 1, sum) # Total biomass
SSB[,y] = apply(NAA[,y,,1,drop = FALSE] * WAA[,y,,1,drop = FALSE] * MatAA[,y,,1,drop = FALSE] * exp(-ZAA[,y,,1,drop = FALSE] * t_spawn), 1, sum)
Dynamic_SSB0[,y] = apply(NAA0[,y,,1,drop = FALSE] * WAA[,y,,1,drop = FALSE] * MatAA[,y,,1,drop = FALSE] * exp(-natmort[,y,,1,drop = FALSE] * t_spawn), 1, sum)

# If single sex model, multiply SSB calculations by 0.5
if(n_sexes == 1) {
  SSB[,y] = SSB[,y] * 0.5
  Dynamic_SSB0[,y] = Dynamic_SSB0[,y] * 0.5
}
} # end y loop
