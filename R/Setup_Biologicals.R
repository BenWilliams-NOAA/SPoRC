#' Set up simulation stuff for biologicals
#'
#' @param base_M_value Base natural mortality value dimensioned by regions, ages, sexes
#' @param M_pattern Natural mortality pattern. Options include: constant
#' @param base_WAA_values Base weight-at-age values dimensioned by regions, ages, sexes
#' @param base_WAA_fish_values Base weight-at-age values for the fishery dimensioned by regions, ages, sexes, fishery fleets
#' @param WAA_pattern Weight-at-age pattern. Options include: constant
#' @param base_Maturity_AA_values Base maturity values dimensioned by regions, ages, sexes
#' @param Maturity_AA_pattern Maturity pattern. Options include: constant
#' @param sim_list Simulation list objects
#'
#' @export Setup_Sim_Biologicals
Setup_Sim_Biologicals <- function(
                                  base_M_value,
                                  M_pattern,
                                  base_WAA_values,
                                  base_WAA_fish_values,
                                  WAA_pattern,
                                  base_Maturity_AA_values,
                                  Maturity_AA_pattern,
                                  sim_list
                                  ) {

  # Create containers to store values
  M <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  WAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  WAA_fish <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  Maturity_AA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

  for(sim in 1:sim_list$n_sims) {
    for(r in 1:sim_list$n_regions) {
      for(y in 1:sim_list$n_yrs) {
        for(s in 1:sim_list$n_sexes) {

          # Natural mortality constant
          if(M_pattern == "constant") M[r,y,,s,sim] <- base_M_value[r,,s]

          # WAA constant
          if(WAA_pattern == "constant") WAA[r,y,,s,sim] <- base_WAA_values[r,,s]
          if(WAA_pattern == "constant") for(f in 1:sim_list$n_fish_fleets) WAA_fish[r,y,,s,f,sim] <- base_WAA_fish_values[r,,s,f]

          # Maturity constant
          if(Maturity_AA_pattern == "constant") Maturity_AA[r,y,,s,sim] <- base_Maturity_AA_values[r,,s]

        } # end s loop
      } # end y loop
    } # end r loop
  } # end sim loop

  # output into list
  sim_list$M <- M
  sim_list$WAA <- WAA
  sim_list$WAA_fish <- WAA_fish
  sim_list$Maturity_AA <- Maturity_AA

  return(sim_list)

}

#' Helper function to map natural mortality blocks
#'
#' This function maps natural mortality (\code{ln_M}) to a block structure across region, year, age, and sex dimensions.
#' It assigns unique integer identifiers to each block defined by the user's specifications, which is stored in the \code{M_blocks}
#' array in the input list.
#'
#' @param input_list A named list object containing model data, parameters, and mapping structures.
#' @param M_spec Character string indicating whether to estimate or fix natural mortality.
#'   Options are:
#'   \itemize{
#'     \item \code{"est_ln_M"}: Estimate natural mortality parameters.
#'     \item \code{"fix"}: Fix all natural mortality parameters (requires them to be passed in the input list).
#'   }
#' @param M_regionblk_spec_vals A list of numeric vectors specifying the region indices grouped in each block.
#' @param M_yearblk_spec_vals A list of numeric vectors specifying the year indices grouped in each block.
#' @param M_ageblk_spec_vals A list of numeric vectors specifying the age indices grouped in each block.
#' @param M_sexblk_spec_vals A list of numeric vectors specifying the sex indices grouped in each block.
#'
#' @return An updated \code{input_list} with mapped natural mortality blocks:
#' \itemize{
#'   \item \code{input_list$map$ln_M} is a factor indicating fixed or estimated mortality parameters.
#'   \item \code{input_list$data$M_blocks} is a 4D array (region × year × age × sex) with unique block IDs.
#' }
#'
#' @keywords internal
do_M_mapping <- function(input_list,
                         M_spec,
                         M_regionblk_spec_vals,
                         M_yearblk_spec_vals,
                         M_ageblk_spec_vals,
                         M_sexblk_spec_vals) {

  # Validate options
  if(!M_spec %in% c('est_ln_M', 'fix')) stop("M_spec needs to be specified as either est_ln_M (only for a single sex) or fix")

  # set up whether fixing M or estimating
  if(M_spec == 'est_ln_M') input_list$map$ln_M <- factor(1:length(input_list$par$ln_M))
  if(M_spec == 'fix') input_list$map$ln_M <- factor(rep(NA, length(input_list$par$ln_M)))

  # create array for blocks
  M_blocks <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))

  # loop through to get counters for blocking structure for indexing
  counter <- 1
  for (regionblk in 1:length(M_regionblk_spec_vals)) {
    map_r <- M_regionblk_spec_vals[[regionblk]]

    for (yearblk in 1:length(M_yearblk_spec_vals)) {
      map_y <- M_yearblk_spec_vals[[yearblk]]

      for (ageblk in 1:length(M_ageblk_spec_vals)) {
        map_a <- M_ageblk_spec_vals[[ageblk]]

        for (sexblk in 1:length(M_sexblk_spec_vals)) {
          map_s <- M_sexblk_spec_vals[[sexblk]]

          # Assign the current counter to this block
          M_blocks[map_r, map_y, map_a, map_s] <- counter
          counter <- counter + 1

        } # end sexblk
      } # end ageblk
    } # end yearblk
  } # end regionblk

  collect_message("Natural Mortality specified as: ", M_spec)
  collect_message("Natural Mortality Region Blocks is specified as: ", length(M_regionblk_spec_vals))
  collect_message("Natural Mortality Year Blocks is specified as: ", length(M_yearblk_spec_vals))
  collect_message("Natural Mortality Age Blocks is specified as: ", length(M_ageblk_spec_vals))
  collect_message("Natural Mortality Sex Blocks is specified as: ", length(M_sexblk_spec_vals))

  input_list$data$M_blocks <- M_blocks

  return(input_list)
}

#' Setup biological inputs for estimation model
#'
#' @param input_list List containing data, parameter, and map lists for the model.
#' @param WAA Numeric array of weight-at-age (spawning), dimensioned \code{[n_regions, n_years, n_ages, n_sexes]}.
#' @param MatAA Numeric array of maturity-at-age, dimensioned \code{[n_regions, n_years, n_ages, n_sexes]}.
#' @param AgeingError Numeric matrix or array representing the ageing error transition matrix.
#'   If a matrix (2D), dimensions should be \code{[number of modeled ages, number of observed composition ages]}
#'   and the ageing error is assumed to be constant over time.
#'   If an array (3D), dimensions should be \code{[number of years, number of modeled ages, number of observed composition ages]}
#'   allowing ageing error to vary by year.
#'   Defaults to an identity matrix (no ageing error) if not specified, assuming observed age bins exactly match modeled age bins.
#'
#'   **Note:** If the observed age composition bins differ from the modeled age bins
#'   (e.g., observed ages 2–10 while modeled ages are 1–10), the default identity matrix will cause a dimensional mismatch
#'   and misalignment. In such cases, users should provide a custom ageing error matrix mapping modeled to observed ages.
#'   For example, to drop the first modeled age bin, supply a matrix like \code{diag(1, 10)[, 2:10]}.
#'   This ensures proper alignment of age bins for likelihood calculations.
#' @param Use_M_prior Integer flag indicating whether to apply a natural mortality prior (\code{0} = no, \code{1} = yes).
#' @param M_prior Numeric vector of length two giving the mean (in normal space) and standard deviation of the natural mortality prior.
#' @param fit_lengths Integer flag indicating whether to fit length data (\code{0} = no, \code{1} = yes).
#' @param SizeAgeTrans Numeric array of size-at-age transition probabilities, dimensioned \code{[n_regions, n_years, n_lens, n_ages, n_sexes]}.
#' @param M_spec Character string specifying natural mortality estimation approach. Defaults to \code{est_ln_M}, which estimates mortality to be invariant, if blocks are not specified. Options:
#' \itemize{
#'   \item \code{"est_ln_M"}: Estimates natural mortality across the defined natural mortality blocks.
#'   \item \code{"fix"}: Fix all natural mortality parameters using the provided array.
#' }
#' @param Fixed_natmort Numeric array of fixed natural mortality values, dimensioned \code{[n_regions, n_years, n_ages, n_sexes]}. Required if \code{M_spec = "fix"}.
#' @param Selex_Type Character string specifying whether selectivity is age or length-based. Default is age-based
#' \itemize{
#'   \item \code{"length"}: Length-based selectivity.
#'   \item \code{"age"}: Age-based selectivity
#' }
#' @param WAA_fish Numeric array of weight-at-age (fishery), dimensioned \code{[n_regions, n_years, n_ages, n_sexes, n_fish_fleets]}.
#' @param WAA_srv Numeric array of weight-at-age (survey), dimensioned \code{[n_regions, n_years, n_ages, n_sexes, n_srv_fleets]}.
#' @param addtocomp Numeric value for a constant to add to composition data. Default is 1e-3.
#' @param M_ageblk_spec Specification of age blocking for natural mortality estimation.
#'   Either a character string ("constant") or a list of index vectors, e.g., \code{list(1:10, 11:30)}, which specifies 2 age blocks for M.
#' @param M_regionblk_spec Specification of regional blocking for natural mortality.
#'   Either a character string ("constant") or a list of index vectors, e.g., \code{list(1:3, 4:5)}, which specifies 2 region blocks for M.
#' @param M_yearblk_spec Specification of year blocking for natural mortality.
#'   Either a character string ("constant") or a list of index vectors, e.g., \code{list(1:10, 11:30)}, which specifies 2 year blocks for M.
#' @param M_sexblk_spec Specification of sex blocking for natural mortality.
#'   Either a character string ("constant") or a list of index vectors, e.g., \code{list(1:2)}, which specifies sex-invariant M.
#' @param ... Additional arguments for starting values such as \code{ln_M} and \code{M_offset.} These are ignored if \code{M_spec = fix}.
#'
#' @export Setup_Mod_Biologicals
Setup_Mod_Biologicals <- function(input_list,
                                  WAA,
                                  WAA_fish = NULL,
                                  WAA_srv = NULL,
                                  MatAA,
                                  addtocomp = 1e-3,
                                  AgeingError = NULL,
                                  Use_M_prior = 0,
                                  M_prior = NA,
                                  fit_lengths = 0,
                                  SizeAgeTrans = NA,
                                  Selex_Type = 'age',
                                  M_spec = "est_ln_M",
                                  M_ageblk_spec = 'constant',
                                  M_regionblk_spec = 'constant',
                                  M_yearblk_spec = 'constant',
                                  M_sexblk_spec = 'constant',
                                  Fixed_natmort = NULL,
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Weight at age checking
  check_data_dimensions(WAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'WAA')
  if(!is.null(WAA_fish)) check_data_dimensions(WAA_fish, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'WAA_fish')
  if(!is.null(WAA_srv)) check_data_dimensions(WAA_srv, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'WAA_srv')

  # Maturity at age checking
  check_data_dimensions(MatAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'MatAA')

  # Length checking
  if(!fit_lengths %in% c(0,1)) stop("Values for fit_lengths are not valid. They are == 0 (not used), or == 1 (used)")
  collect_message("Length Composition data are: ", ifelse(fit_lengths == 0, "Not Used", "Used"))

  # Size Age Transition checking
  if(fit_lengths == 1) check_data_dimensions(SizeAgeTrans, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'SizeAgeTrans')
  if(fit_lengths == 1 & is.na(sum(SizeAgeTrans))) stop("Length composition are fit to, but the size-age transition matrix is NA")

  # Natural Mortality checking
  if(!is.null(M_spec)) {
    if(M_spec == 'fix') {
      if(is.null(Fixed_natmort)) stop("Please provide a fixed natural mortality array dimensioned by n_regions, n_years, n_ages, and n_sexes!")
      check_data_dimensions(Fixed_natmort, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'Fixed_natmort')
    }
  }

  # Natural Mortality prior checking
  if(!Use_M_prior %in% c(0,1)) stop("Values for Use_M_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  collect_message("Natural Mortality priors are: ", ifelse(Use_M_prior == 0, "Not Used", "Used"))

  # Checking ageing error dimensions
  if(!is.null(AgeingError)) {
    if(length(dim(AgeingError)) == 2) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), what = 'AgeingError') # user supplied ageing error is not time-varying
    if(length(dim(AgeingError)) == 3) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), n_years = length(input_list$data$years), what = 'AgeingError_t') # user supplied ageing error is time-varying
  }

  # Selectivity Options -----------------------------------------------------

  # Age based selectivity
  if(Selex_Type == 'age') {
    Selex_Type <- 0
    collect_message("Selectivity is aged-based.")
  } # if age based

  # Length based selectivity
  if(Selex_Type == 'length') {
    if(fit_lengths == 0) stop("Length composition data are not fit, but selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    Selex_Type <- 1
    collect_message("Selectivity is length-based")
  } # if length based


  # Weight at Age Options ---------------------------------------------------

  # setup fishery and survey specific weight at age (if not specified - just uses the WAA (spawning) already supplied)
  if(is.null(WAA_fish)) { # if no fishery WAA provided, use spawning WAA supplied
    WAA_fish <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(f in 1:input_list$data$n_fish_fleets) WAA_fish[,,,,f] <- WAA
    collect_message("WAA_fish was specified at NULL. Using the spawning WAA for WAA_fish")
  }

  # if no survey WAA provided, use spawning WAA supplied
  if(is.null(WAA_srv)) {
    WAA_srv <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(f in 1:input_list$data$n_srv_fleets) WAA_srv[,,,,f] <- WAA
    collect_message("WAA_srv was specified at NULL. Using the spawning WAA for WAA_srv")
  }

  # Ageing Error Options ----------------------------------------------------

  # setup ageing error if not provided
  if(is.null(AgeingError)) {
    AgeingError <- diag(1, length(input_list$data$ages)) # if no inputs for ageing error, then create identity matrix
    AgeingError_t <- array(0, dim = c(length(input_list$data$years), dim(AgeingError)))
    for(i in 1:length(input_list$data$years)) AgeingError_t[i,,] <- AgeingError
    warning("No ageing error matrix was provided. A default identity matrix was used, which assumes that the number and structure of modelled age bins exactly match the observed age bins. If the observed age composition data includes fewer age bins than the model (e.g., observed ages 2-10 while modelled ages are 1-10), this default assumption will cause a dimensional mismatch and potentially misalign the modelled and observed compositions. To avoid this, please provide an ageing error matrix of dimension n_model_ages x n_obs_ages that correctly maps modelled ages to observed age bins. For example, if observed ages are 2-10, supply a matrix that drops the first model age by using a shifted identity matrix: diag(1, 10)[, 2:10]. This will ensure the age bins are correctly aligned for likelihood calculations.")
  }

  # setup ageing error if user-supplied is not year specific
  if(!is.null(AgeingError) && length(dim(AgeingError)) == 2) {
    AgeingError_t <- array(0, dim = c(length(input_list$data$years), dim(AgeingError)))
    for(i in 1:length(input_list$data$years)) AgeingError_t[i,,] <- AgeingError
    collect_message("Ageing Error is specified to be time-invariant")
  }

  # ageing error if it is year specific (just reassigning)
  if(!is.null(AgeingError) && length(dim(AgeingError)) == 3) {
    AgeingError_t <- AgeingError
    collect_message("Ageing Error is specified to be time-varying")
  }


  # Natural Mortality Options -----------------------------------------------
  # Input indicator for estimating or not estimating M
  if(is.null(M_spec) || M_spec == "est_ln_M") input_list$data$use_fixed_natmort <- 0
  else if(M_spec == "fix") input_list$data$use_fixed_natmort <- 1

  # Populate Data List ------------------------------------------------------

  input_list$data$WAA <- WAA
  input_list$data$WAA_fish <- WAA_fish
  input_list$data$WAA_srv <- WAA_srv
  input_list$data$MatAA <- MatAA
  input_list$data$AgeingError <- AgeingError_t
  input_list$data$fit_lengths <- fit_lengths
  input_list$data$SizeAgeTrans <- SizeAgeTrans
  input_list$data$Use_M_prior <- Use_M_prior
  input_list$data$M_prior <- M_prior
  input_list$data$Fixed_natmort <- Fixed_natmort
  input_list$data$Selex_Type <- Selex_Type
  input_list$data$addtocomp <- addtocomp

  # Populate Parameter List -------------------------------------------------

  # If M is constant for ages
  if(is.character(M_ageblk_spec)){
    if(M_ageblk_spec == "constant") M_ageblk_spec_vals <- list(1:length(input_list$data$ages))
  } else M_ageblk_spec_vals <- M_ageblk_spec

  # If M is constant across years
  if(is.character(M_yearblk_spec)){
    if(M_yearblk_spec == "constant") M_yearblk_spec_vals = list(1:length(input_list$data$years))
  } else M_yearblk_spec_vals = M_yearblk_spec

  # If M is constant across sexes
  if(is.character(M_sexblk_spec)){
    if(M_sexblk_spec == "constant") M_sexblk_spec_vals <- list(1:input_list$data$n_sexes)
  } else M_sexblk_spec_vals <- M_sexblk_spec

  # If M is constant across regions
  if(is.character(M_regionblk_spec)){
    if(M_regionblk_spec == "constant") M_regionblk_spec_vals <- list(1:input_list$data$n_regions)
  } else M_regionblk_spec_vals <- M_regionblk_spec

  if("ln_M" %in% names(starting_values)) input_list$par$ln_M <- starting_values$ln_M
  else input_list$par$ln_M <- array(log(0.5), dim = c(length(M_regionblk_spec_vals), length(M_yearblk_spec_vals),
                                                      length(M_ageblk_spec_vals), length(M_sexblk_spec_vals)))

  # Mapping Options ---------------------------------------------------------
  input_list <- do_M_mapping(input_list, M_spec, M_regionblk_spec_vals,
                             M_yearblk_spec_vals, M_ageblk_spec_vals, M_sexblk_spec_vals) # natural mortality mapping

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

