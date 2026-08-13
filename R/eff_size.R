#' Calculate effect sizes
#'
#' Blah blah blah
#'
#' @import data.table
#' @import crayon
#' @importFrom tools R_user_dir
#' @importFrom parallel mclapply detectCores
#' @param ... variables required by effect size calculation (e.g., x, x1, x2, sd1, sd2, see Description)
#' @param effect_type Single length character specifying effect size. Currently supported effect sizes include: SMD, SMDH, lnRoM, lnCVR, lnVR, lnM, Zr, lnHWE
#' @param data An optional data frame containing the variables necessary for calculating effect size.
#' @param bind Should the function bind the effect sizes to the inputted data frame? Logical
#' @param paired Use paired version of effect size? If so, then inputs will differ from unpaired (e.g., `n` instead of `n1` and `n2`).
#' @param default_formulas Logical. Whether to return default estimators. If FALSE (recommended and the function's default), will return first and second derivative estimators. Note that second refers to bias-corrected estimates (which are not always less biased). See Details.
#' @param SAFE Logical. Whether to conduct SAFE bootstrapping of effect size. This is essential to calculate lnM in boundary scenarios
#' @param SAFE_boots Numeric. Number of bootstraps. Default is 1e6.
#' @param SAFE_max_secs Numeric. Number of seconds to conduct SAFE bootstrapping before timing out. See Details.
#' @param n_cores Numeric. Number of cores for parallel processing of SAFE calculation. See Details.
#' @param SAFE_distribution Character. Distribution to use for SAFE bootstrapping. If unspecified, will use default distribution. See Details.
#' @param sigma_matrix Matrix. Optional custom sigma_matrix for SAFE bootstrapping (experimental)
#' @param verbose How chatty do you want the function to be? Logical
#' @return A data.table with effect sizes and sample variances.
#' @export
eff_size <- function(...,
                     effect_type = NULL,
                     data = NULL,
                     bind = TRUE,
                     paired = FALSE,
                     default_formulas = TRUE,
                     formula_path = "data/effect_size_formulas.csv",
                     SAFE = FALSE,
                     SAFE_boots = 1e6,
                     SAFE_max_secs = 15,
                     n_cores = 1,
                     SAFE_distribution = NULL,
                     sigma_matrix = NULL,
                     verbose = T
                     ){

    env <- parent.frame()

    # Set up parallel for SAFE
    pbop <- pbapply::pboptions(type = "txt")
    pbapply::pboptions(pbop)

    #
    if(!is.null(data)) dat <- copy(data)

    # Load formulas -----------------------------------------------------------
    # load("data/effect_formulas.rda")

    data.table::setorder(effect_formulas, effect_size, calc_type)

    # Check that effect type is specified and filter formulas
    if(is.null(effect_type)){

      cat(blue(("\nEffect type name must be specified with 'effect_type' argument and provide necessary variables (named in arguments to function call) to match formula equations.\n")),
          blue("\nReturning effect size names & required variables for reference.\n\n"))
      return(unique(effect_formulas[, .(effect_size, paired_design, vars_required)]))

    }else if(effect_type %in% effect_formulas$effect_size){
      # filter to desired effect_type  and calculation
      effect_formulas.sub <- effect_formulas[effect_size == effect_type, ]

      if(any(paired %in% effect_formulas.sub$paired_design)) effect_formulas.sub <- effect_formulas.sub[paired_design == paired, ]

    }else if(!effect_type %in% effect_formulas$effect_size){

      cat(blue("Effect type name misspecified"),
          blue("\nReturning effect size names & required variables for reference.\n\n"))
      return(unique(effect_formulas[, .(effect_size, vars_required)]))

    }

    # If SAFE is unsupported let the user know and set SAFE to FALSE
    if(all(is.na(effect_formulas.sub$SAFE_family))){
      SAFE <- FALSE
      cat(blue("SAFE is currently not implemented for this effect size"))
    }

    if(SAFE == FALSE){
      # Drop extra rows for multiple SAFE methods:
      effect_formulas.sub <- effect_formulas.sub[default_safe_family %in% c("yes"), ]
    }

    # Get the required variables:
    vars <- strsplit(unique(effect_formulas.sub$vars_required), split = ", ") |>
      unlist()

    # >>> Parse inputs ----------------------------------------------------
    # If debugging, skip this section

    # Function can now accept vectors or NSE inputs (unquoted column names)
    call_expr <- match.call(expand.dots = TRUE)
    args <- as.list(call_expr)[-1]
    if("" %in% names(args)) stop("Numeric arguments must be named")

    args <- args[vars[vars %in% names(args)]] # because of optionally specified r
    input_vars <- list()

    #
    if(!is.null(data)){
      setDT(dat)
      for(i in 1:length(args)){
        input_vars[[i]] <- dat[, eval(args[[i]], envir = env)]
      }

    }else if(is.vector(eval(args[[1]], envir = env))){
      for(i in 1:length(args)){
        input_vars[[i]] <- eval(args[[i]], envir = env)
      }
    }
    names(input_vars) <- names(args)

    if(length(unique(lengths(input_vars))) > 1){ stop(cat("Input vectors", "(", red(paste(names(input_vars), collapse = ", ")), ")",  "are different lengths. Please double check inputs.")) }

    # >>> Preliminary checks and filtering --------------------------------------------------
    # Deal with missing 'r'
    if(paired == TRUE &
       !"r" %in% names(input_vars)){

      cat("Paired design selected", red("but 'r' not specified."), "Setting 'r' to 0.5\n")
      input_vars$r <- rep(0.5, max(lengths(input_vars)))


    }else if(paired == FALSE &
             !"r" %in% names(input_vars)){

      input_vars$r <- rep(0, max(lengths(input_vars))) # This is necessary for the shared sigma_matrices of some effect sizes

    }

    # Check for missing variables.
    if(!all(vars %in% names(input_vars))){
      return(cat("Missing the following variables:",
                 red(paste(setdiff(vars, names(input_vars)), collapse=", ")), "\n"))
    }

    # Print effect size specific warnings, e.g., 0 in lnOR and lnRR
    if(!is.na(unique(effect_formulas.sub$special_warnings)) &
       verbose == TRUE){
      cat(unique(effect_formulas.sub$special_warnings), "\n",
          "Leaving it to user's discretion to check prior to execution.\n\n")
    }

    # Deal with alternative SAFE distributions.
    if(is.null(SAFE_distribution) &
       "yes" %in% effect_formulas.sub$default_safe_family &
       SAFE == TRUE){
      # If unspecified (SAFE_distribution == NULL & there are multiple options for default, then choose default
      effect_formulas.sub <- effect_formulas.sub[default_safe_family == "yes", ]
    }else if(!is.null(SAFE_distribution)){
      # If SAFE_distribution is specified, subset to SAFE_distribution
      effect_formulas.sub <- effect_formulas.sub[SAFE_family == SAFE_distribution, ]
    }
    # If unspecified (SAFE_distribution == NULL & effect_formulas.sub$default is all NA then do nothing)

    if(nrow(effect_formulas.sub) == 0){
      return(cat(red("\nEffect size not available after filtering to type."),
                 "\n\nEffect sizes currently supported include:", paste(sort(unique(effect_formulas$effect_size)), collapse = "; "),
                 blue("\n\nTo add custom effect sizes please see XXXX")) )
    }


    # >>> Filter defaults ---------------------------------------------------------
    # For SAFE:
    if(SAFE == TRUE) definition_formula <- effect_formulas.sub[derivative == "first" & calc_type == "point_estimate", ]

    if(default_formulas == TRUE) effect_formulas.sub <- effect_formulas.sub[default == "yes", ]

    # >>> Calculate plugin effect size: -------------------------------------------------
    if(verbose){
      cat("Using the formulas:\n\t", blue(paste(
              unique(unlist(strsplit(effect_formulas.sub$formula, "; "))),
              collapse = "\n\t ")),
          "\nBe sure that all variables in formula are correctly named.\n\n")
    }

    plugins <- suppressWarnings(.calc_effect(effect_formulas.sub, input_vars))

    if(any(is.na(unlist(plugins)))) cat(magenta("Plugin effect sizes could not be calculated\n\n"))

    if(default_formulas == TRUE) setnames(plugins, names(plugins), gsub("_first|_second", "", names(plugins)))

    if(SAFE == TRUE){

      # >>> SAFE calculation ----------------------------------------------------------------
      # Extract reference plugin effect size. First order.
      definition <- suppressWarnings(.calc_effect(definition_formula, input_vars))
      plugin_effect_size <- definition$yi_first

      index <- seq(1:max(lengths(input_vars)))
      k <- 1

      if(length(plugin_effect_size) != max(index)){ return(cat("Shit.")) }

      # Run SAFE function for each element of input_vars:
      safe_out <- pbapply::pblapply(index, function(k){
        return(.SAFE_calc(formulas = definition_formula, # Changed to this from `effect_formulas.sub`
                          input = lapply(input_vars, "[[", k), # select the first element in each element...
                          plugin_effect = plugin_effect_size[k],
                          custom_sigma = sigma_matrix[[k]], # submit custom sigma_matrix if it exists.
                          SAFE_boots = 1e6,
                          SAFE_max_secs = SAFE_max_secs))
      },
      cl = n_cores) |>
        rbindlist()

      safe_out[, SAFE_complete := ifelse(number_SAFE_bootstraps == SAFE_boots, "yes", "no")]
      if(any(safe_out$number_SAFE_bootstraps < SAFE_boots) &
         !any(is.na(safe_out$yi_safe))){

        cat(magenta("\n\nBoundary issues prevented full number of SAFE bootstraps. Try reducing `SAFE_boots` (1e6 is the default) or increasing time limit `SAFE_max_secs` (15 seconds is default).\n\n"))

      }else if(any(is.na(safe_out$yi_safe))){

        cat(magenta("SAFE could not be calculated for at least one of the inputs\n\n"))

      }

      out <- cbind(plugins, safe_out)

    }else if(SAFE == FALSE){
      out <- plugins
    }

    # >>> Return objects ------------------------------------------------------
    if(bind == TRUE & !is.null(data)){
      out <- cbind(data, out)
    }

    return(out)
}
