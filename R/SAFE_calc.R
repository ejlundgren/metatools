#' DEPRECATED. Calculate effect sizes
#'
#' Used internally to calculate SAFE effect sizes. This function calls parameter_cloud and then transforms hyperparameters
#'
#' @import data.table
#' @importFrom stats sd
#' @param formulas blah
#' @param input_k blah
#' @param plugin_effect_k blah
#' @param sigma_matrix_k blah
#' @param SAFE_boots blah
#' @return A data.table with effect sizes and sample variances
#' @export
# SAFE_calc <- function(formula,
#                       cloud,
#                       plugin_effect_k,
#                       sigma_matrix_k,
#                       SAFE_boots = 1e6){
#
#
#   # Convert cloud
#   cloud_trans <- effsize::calc_effect(formulas = definition_formula,
#                                       input = cloud)$yi_definition
#
#   # bias corrected estimate of sampling variance and SE:
#   safe_SE <- stats::sd(cloud_trans)
#   safe_vi <- safe_SE^2
#
#   bias_SAFE <- mean(cloud_trans) - plugin_effect_k
#
#   safe_yi <- plugin_effect_k - bias_SAFE
#
#   return(data.table(yi_safe = safe_yi,
#                     vi_safe = safe_vi,
#                     SE_safe = safe_SE))
# }
