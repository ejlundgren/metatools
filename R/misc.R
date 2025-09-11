#
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -----------------------------------------
#
# # TESTING -----------------------------------------------------------------

formula <- "(sd1^2 / (n1 * x1^2)) + (sd2^2 / (n2 * x2^2)) - ((2 * r * sd1 * sd2) / (x1 * x2 * sqrt(n1 * n2)))"
variable_names <- c("x1", "x2", "sd1", "sd2", "n1", "n2")

.create_exec_formula <- function(formula, variable_names){
  to <- paste0("input$", variable_names)
  from <- paste0("\\b", variable_names, "\\b")

  exec_formula <- formula

  for(k in 1:length(vars)){
    exec_formula <- gsub(vars[k], to[k], exec_formula)
  }
  return(exec_formula)
}

