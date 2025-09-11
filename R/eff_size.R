#' Calculate effect sizes
#'
#' Blah blah blah
#'
#' @import data.table
#' @import crayon
#' @importFrom tools R_user_dir
#' @importFrom parallel mclapply detectCores
#' @param ... variables required by effect size calculation (e.g., x, x1, x2, sd1, sd2, see Description)
#' @param effect_type Character of effect size. Currently supported effect sizes include SMD, lnRoM, lnCVR, OR, RR, blah blah
#' @param SAFE Logical. Whether to conduct SAFE bootstrapping of effect size
#' @param SAFE_boots Numeric. Number of bootstraps (default is 1e6)
#' @param SAFE_distribution Character. distribution to use for SAFE bootstrapping. See description
#' @param sigma_matrix Matrix. Optional custom sigma_matrix for SAFE bootstrapping.
#' @param verbose Verbose or not? Logical
#' @param parallelize Logical.
#' @param use_custom_formulas Logical. Load user-customized formulas? These are saved to a local path (tools::R_user_dir(package = "effsize", which = "data")) with the add_custom_effect_sizes() function.
#' @return A data.table with effect sizes and sample variances
#' @export
eff_size <- function(..., # This is where the input variables are passed in.
                     effect_type = NULL,
                     paired_design = FALSE,
                     SAFE = TRUE,
                     SAFE_boots = 1e6,
                     SAFE_distribution = NULL, # Aug 2025: Looks like some of these should provide a choice...
                     sigma_matrix = NULL, # Custom sigma matrix. Needs to be a list calculated off of data of hte same length as input_vars. Maybe down the road this could be a custom function
                     verbose = T,
                     parallelize = T,
                     use_custom_formula_table = FALSE,
                     # So users can add formula here instead of to table:
                     # Which is maybe a better plan anyways...
                     definition_formula = NULL,
                     cloud_filtering_rules = NULL
                     ){

    input_vars <- list(...)

    # Prepare formula table ---------------------------------------------------

    if(!is.null(definition_formula)){
      if(is.null(cloud_filtering_rules) | is.null(SAFE_distribution)) return(cat(crayon::red("'definition_formula' specified. Must specify cloud_filtering rules and SAFE_distribution type.")))
      effect_formulas.sub <- data.table(formula = definition_formula,
                                        SAFE_family = SAFE_distribution,
                                        cloud_filtering_rules = cloud_filtering_rules,
                                        formula_type = "yi_definition")
      effect_formulas.sub[, exec_formula := .create_exec_formula(formula, names(input_vars))]

    }else{
      # >>> Prepare premade formula table -----------------------------------------------
      path <- file.path(tools::R_user_dir(package = "effsize", which = "data"),
                        "data",
                        "modified_effect_formulas.rda")

      if(use_custom_formula_table == TRUE & file.exists(path)){
        load(path)
      }else{
        load("data/effect_formulas.rda")
        if(use_custom_formulas == TRUE) cat("No modified formula table present. Using stable formulas.")
      }
      setorder(effect_formulas, name, calc_type)

      # >>> Preliminary checks and filtering --------------------------------------------------

      # Could build some tutorial information into this:
      if(is.null(effect_type)){
        cat(crayon::red(("\nMust specify an effect size type ('effect_type') and necessary variables (named in arguments to function call) to match formula equations.\n")),
            crayon::blue("\nReturning effect size names & required variables for reference.\n\n"))
        return(unique(effect_formulas[, .(name, vars_required)]))
      }else{
        # filter to desired effect_type  and calculation
        effect_formulas.sub <- effect_formulas[name == effect_type, ]
      }

      if(length(unique(lengths(input_vars))) > 1){ return(cat("Input vectors", "(", red(paste(names(input_vars), collapse = ", ")), ")",  "are different lengths. Please double check inputs.")) }

      # Deal with missing 'r'
      if(paired_design == TRUE & !"r" %in% names(input_vars)){
        cat("Paired design selected", crayon::red("but 'r' not specified."), "Setting 'r' to 0.5.\n\nBe sure that n1 == n2 for all observations.")
        input_vars$r <- rep(0.5, max(lengths(input_vars)))
      }else if(paired_design == FALSE & "yes" %in% effect_formulas.sub$accepts_paired_design){
        input_vars$r <- rep(0, max(lengths(input_vars))) # This is necessary for the shared sigma_matrices of some effect size
      }

      # Check for missing variables.
      vars <- strsplit(unique(effect_formulas.sub$vars_required), split = ", ") |> unlist()
      if(length(setdiff(vars, names(input_vars))) > 0){
        return(cat("Missing the following variables:",
                   crayon::red(paste(setdiff(vars, names(input_vars)), collapse=", ")), "\n"))
      }

      # Print effect size specific warnings, e.g., 0 in lnOR and lnRR
      if(!is.na(unique(effect_formulas.sub$special_warnings)) & verbose == TRUE){
        cat(crayon::magenta(unique(effect_formulas.sub$special_warnings)),
            "\nLeaving it to user's discretion to check prior to execution.\n")
      }

      # Deal with alternative SAFE distributions.
      if(is.null(SAFE_distribution) & "yes" %in% effect_formulas.sub$default_SAFE_family){
        # If unspecified (SAFE_distribution == NULL & there are multiple options for default, then choose default
        effect_formulas.sub <- effect_formulas.sub[default_SAFE_family == "yes", ]
      }else if(!is.null(SAFE_distribution)){
        # If SAFE_distribution is specified, subset to SAFE_distribution
        effect_formulas.sub <- effect_formulas.sub[SAFE_family == SAFE_distribution, ]
      }

      # If unspecified (SAFE_distribution == NULL & effect_formulas.sub$default is all NA then do nothing)
      if(nrow(effect_formulas.sub) == 0){
        return(cat(crayon::red("\nEffect size not available after filtering to type."),
                   "\n\nEffect sizes currently supported include:", paste(sort(unique(effect_formulas$name)), collapse = "; "),
                   crayon::blue("\n\nTo add custom effect sizes please see add_custom_effect_sizes")) )
      }

    }

  # Calculate plugin effect size: -------------------------------------------------
  if(verbose){
    effect_formulas.sub[, to_console := paste0(formula_type, " <- ", formula)]
    cat("Using the formulas:\n\t", crayon::blue(paste(effect_formulas.sub$to_console, collapse = "\n\t ")),
        "\nBe sure that all variables in formula are correctly named.\n\n")
  }

  # effsize::
  plugins <- calc_effect(effect_formulas.sub,
                         input_vars)

  if(SAFE == FALSE){
    return(plugins)
  }else if(!"yi_definition" %in% effect_formulas.sub$formula_type){
    return(cat("'yi_definition' formula required to conduct SAFE calculation. If using custom formulas, please read vignette(). If this error is experienced with stable effect sizes as published, please file an Issue at github.com/ejlundgren/XXXXX"))
  }

  # SAFE calculation ----------------------------------------------------------------

  # Extract definition effect size.
  plugin_effect_size <- plugins$yi_definition
  definition_formula <- effect_formulas.sub[formula_type == "yi_definition", ]

  # Run SAFE function for each element of input_vars:
  SAFE_out <- parallel::mclapply(1:length(plugin_effect_size), function(k){
      input_k <- lapply(input_vars, "[[", k)

      # effsize::
      cloud <- parameter_cloud(family = unique(effect_formulas.sub$SAFE_family),
                               input = input_k,
                               cloud_filtering_rule <- unique(effect_formulas.sub$cloud_filtering_rules),
                               sigma_matrix = sigma_matrix[[k]], #' if specified by user. Otherwise calculated based on sim_family
                               SAFE_boots = SAFE_boots)

      # Add missing inputs (e.g., n)
      cloud <- data.table::data.table(cloud,
                                      input_k[!names(input_k) %in% names(cloud)] |>
                                        unlist() |>
                                        t() |>
                                        data.table())

      # Convert cloud
      #effsize::
      cloud_trans <- calc_effect(formulas = definition_formula,
                                 input = cloud)

      # bias corrected estimate of sampling variance and SE:
      safe_SE <- stats::sd(cloud_trans$yi_definition)
      safe_vi <- safe_SE^2

      bias_SAFE <- mean(cloud_trans$yi_definition) - plugin_effect_size[k]

      safe_yi <- plugin_effect_k - bias_SAFE

      return(data.table(yi_safe = safe_yi,
                               vi_safe = safe_vi,
                               SE_safe = safe_SE))
    },
    mc.cores = ifelse(parallelize == TRUE,
                         (parallel::detectCores()-1),
                         1),
      mc.allow.recursive = TRUE)

  out <- cbind(plugins, data.table::rbindlist(SAFE_out))
  return(out)
}
