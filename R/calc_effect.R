#' Evaluate formulas
#'
#' Used internally
#'
#' @param formulas Formula to evalulate
#' @param input List of named input variables to evaluate formula on
#' @return A data.table with effect sizes and sample variances
#' @export
calc_effect <- function(formulas,
                        input){
  # Concatenate the formulas into a single formula, separated with ';'
  exec <- paste(paste(formulas$formula_type, "<-", formulas$exec_formula), collapse = "; ")

  # This adds the effects/variances to the local env but with name assignation:
  eval(parse(text = exec))

  # This gathers them:
  out <- eval(parse(text = paste0("data.table(", paste(unique(formulas$formula_type), collapse = ", "), ")")))

  return(out)
}

