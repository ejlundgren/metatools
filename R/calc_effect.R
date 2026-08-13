#' Evaluate formulas
#'
#' Used internally to evaluate mathematical formulas.
#'
#' @param formulas Formula to evalulate
#' @param input List of named input variables to evaluate formula on
#' @return A data.table with effect sizes and sample variances
.calc_effect <- function(formulas,
                        input){
  # Concatenate the formulas into a single formula, separated with ';'
  exec <- paste(formulas$exec_formula, collapse = "; ")

  # This adds the effects/variances to the local env but with name assignation:
  eval(parse(text = exec), envir = environment())

  res_list <- lapply(unique(formulas$label), function(x)
    get(x, envir = environment()))
  names(res_list) <- unique(formulas$label)

  return(data.table::as.data.table(res_list))
}

