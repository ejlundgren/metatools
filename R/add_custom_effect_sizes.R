#' Calculate effect sizes
#'
#' Blah blah blah
#'
#' @import data.table
#' @import crayon
#' @importFrom tools R_user_dir
#' @param derived_formula Vector of formulas as character strings. If mathematical functions are used (e.g., xxxxx()) they must be namespaced.
#' @param name Name of effect size
#' @param formula_type # at least
#' @param calculation_type Character vector of equal length as 'formula' specifying whether formula is for 'effect_size' or 'sampling_variance'. Required
#' @param vars_required xxxx Required. Character string of required input variables. E.g., x1, x2, n1, n2, etc.
#' @param label xxxx Optional character string to label results
#' @param filtering_rules xxxxx (NA or a character string that could be executed as R code. See details)
#' @param special_warnings xxxx (optional)
#' @param SAFE_family Options include: "1_normal", "2_multivariate_normal", "4_multivariate_normal", "4_multivariate_normal_wishart", "2_binomial", "2_multinomial"
#' @param add_to_stable_formulas Add to already modified formulas or to user-defined formulas?
#' @return A data.table with effect sizes and sample variances
#' @export
add_custom_effect_sizes <- function(formula = NULL,
                                    name,
                                    derivative,
                                    calculation_type,
                                    vars_required,
                                    label = NA,
                                    filtering_rules = NA,
                                    special_warnings = NA,
                                    SAFE_family,
                                    add_to_stable_formulas = TRUE){

  # Various checks...

  # Load formulas:
  if(add_to_stable_formulas == FALSE){
    path <- file.path(tools::R_user_dir(package = "effsize", which = "data"),
                      "data",
                      "modified_effect_formulas.rda")
    if(file.exists(path)){
      load(path)
    }else{
      cat("No modified formula table present. Using stable formulas")
    }
  }

  # Create exec_formula
  # Need to gsub 'var' from vars_required with 'input$' . Need for loops because gsub is not vectorized.
  vars <- strsplit(unique(vars_required), split = ", ") |> unlist()
  to <- paste0("input$", vars)
  from <- paste0("\\b", vars, "\\b")

  exec_formula <- formula
  for(k in 1:length(from)){
    exec_formula <- gsub(from[k], to[k], exec_formula)
  }

  new_effect <- data.table::data.table(name,
                                       derivative,
                                       sim_family = SAFE_family,
                                       vars_required,
                                       cloud_filtering_rules = filtering_rules,
                                       special_warnings,
                                       default_safe_family = NA,
                                       calc_type = calculation_type,
                                       formula,
                                       label,
                                       exec_formula,
                                       source = paste("User-specified custom effect", Sys.Date()))
  stopifnot(length(formula) == nrow(new_effect))

  if(all(is.na(new_effect$label))){
    new_effect[, label := ifelse(calculation_type == "effect_size", "yi", "vi")]
    new_effect[!is.na(derivative_order), label := ifelse(derivative_order == "first", paste0(label, "_first"), paste0(label, "_second"))]
  }

  # Testing -------------------------------------------------------------
  # Can probably check that formula works by doing:
  input <- lapply(vars, function(x){
    return(sample(10:100, 1))
  })
  names(input) <- vars

  test <- lapply(exec_formula, function(x){
    tryCatch(expr={
      t <- (eval(parse(text = x)))
    },
    error = function(e){
      print(e)
    })
  }) |> unlist()
  test

  if(any(is.null(test)) | is.character(test) |
     length(test) != length(exec_formula)){
    stop(cat(crayon::red("Novel effect size does not run at index position:",
                         which(lengths(test) == 0),
                         "\nCheck that all variables in formula are in vars_required. And check that formula runs in your local environment.")))
  }


  effect_formulas.mod <- rbind(effect_formulas,
                               new_effect,
                               fill = TRUE)

# >>> Save to user directory ----------------------------------------------
  path <- file.path(tools::R_user_dir(package = "effsize", which = "data"), "data")

  if(!file.exists(path)){
    dir.create(path, recursive = TRUE)
  }

  effect_formulas <- data.table::copy(effect_formulas.mod)
  save(effect_formulas,
       file = file.path(path, "modified_effect_formulas.rda"))
  cat(crayon::blue("Custom formulas saved to"), file.path(path, "modified_effect_formulas.rda"),
      crayon::red("\n\nTo reset, do XYZ..."))
}
