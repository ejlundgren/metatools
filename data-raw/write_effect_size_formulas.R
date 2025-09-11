#
#
# Write and test all effect size formulas
#
#
# dev branch:
# AIM: Redesign so that user formula modifications are easier / more flexible
# AND so that SAFE and plugin calculations are more independent.
# This is because SAFE always requires DEFINITION formula (and only effect size not samplign variance)
# I think the table should thus be WIDE.
#
# name = "SMD"
# yi_definition_formula = "xxxxx"
# yi_bias_corrected_formula = "xxxxxxx"
# vi_first_order_formula = "xxxxxx"
# vi_second_order_formula = "XXXXX"
# SAFE_families_accepted = "4_multivariate_normal, 4_multivariate_normal_wishart"
# SAFE_family_default = "4_multivariate_normal_wishart"
# accepts_paired_designs = TRUE or FALSE
#
# For the new lnRoM formulas:
# name = "lnRoM_lognormal"
# yi_definition = "log(x1/x2)"
# yi_bias_corrected_formula = "log(x) + XXXX - log(y) - XXXXX
# vi_formula = "XXXX"
# SAFE_families_accepted = "2_multivariate_lognormal"
# SAFE_family_default = "2_multivariate_lognormal"
# accepts_paired_designs = TRUE or FALSE
#

# OR:
# THis might be better actually (because 1st and 2nd order variance formulas makes this dataset a monstrosity)
# name = "SMD"
# yi_definition_formula = "xxxxx"
# derived_formula = "xxxxxxx"
# calculation = "vi" or "yi"
# SAFE_families_accepted = "4_multivariate_normal, 4_multivariate_normal_wishart"
# SAFE_family_default = "4_multivariate_normal_wishart"
# accepts_paired_designs = TRUE or FALSE
# HMMMM.I don't like the weird long/wide combo...

# With different cloud_filtering_rules per sim_family, there's no way to make this wide.
# Going to keep it long.
# But going to adjust the names to be clearer.
#
# 0. Prepare environment --------------------------------------------------

rm(list = ls())

library("groundhog")
groundhog.library(pkg = c("data.table",
                          "crayon", #"pryr",
                          "usethis",
                          "MASS", "metafor"),
                  date = "2025-04-15")


# 1. Create effect size table ---------------------------------------------
#' Get rid of _paired
effect_formulas <- data.table(name = c("reciprocal",
                                       "lnRoM",
                                       "lnRoM",
                                       "SMD", # Cohen's d
                                       "SMD",# Hedges' g
                                       "SMD_paired", # Cohen's d
                                       "SMD_paired", # Hedges' g
                                       "SMD", # Cohen's d
                                       "SMD",# Hedges' g
                                       "SMD_paired", # Cohen's d
                                       "SMD_paired", # Hedges' g
                                       "lnRoM_paired",
                                       "lnOR",
                                       "lnRR",
                                       "lnCVR",
                                       "lnCVR",
                                       "lnCVR_paired",
                                       "lnCVR_paired",
                                       "lnCVR",
                                       "lnCVR",
                                       "lnCVR_paired",
                                       "lnCVR_paired",
                                       "lnHWE_A",
                                       "lnM",
                                       "lnM_paired"
                                       ),
                              derivative = c("first", # reciprocal
                                        "first", #lnRoM first
                                        "second",# lnRoM second
                                        "first",# Cohen's d
                                        "second", # Hedgse' g
                                        "first",# Cohen's d paired
                                        "second", # Hedgse' g paired

                                        "first",# Cohen's d WISHART
                                        "second", # Hedgse' g WISHART
                                        "first",# Cohen's d paired WISHART
                                        "second", # Hedgse' g paired WISHART

                                        "first", # lnRoM_paired
                                        "first", # 'lnOR
                                        "first", #lnRR
                                        "first", # lnCVR first order
                                        "second", # lnCVR second order
                                        "first", # lnCVR paired first order
                                        "second", # lnCVR paired second order

                                        "first", # lnCVR first order WISHART
                                        "second", # lnCVR second order WISHART
                                        "first", # lnCVR paired first order WISHART
                                        "second", # lnCVR paired second order WISHART

                                        "first", # lnHWE_A
                                        "first", # lnM
                                        "first" # lnM_paired
                              ),
                              effect_size = c("1 / x",# reciprocal
                                              "log(x1 / x2)",#lnRox first
                                              "log(x1 / x2) + 0.5 * (sd1^2/(n1 * x1^2) - sd2^2/(n2 * x2^2))", # lnRox second
                                              "(x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) )", # Cohen's d
                                              "(1 - 3 / (4 * (n1 + n2 - 2) - 1) ) * ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ) )", # Hedges' g
                                              "(x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) )", # Cohen's d PAIRED
                                              "(1 - 3 / (4 * (n1 + n2 - 2) - 1) ) * ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ) )", # Hedges' g PAIRED
                                              #' [WISHART:]
                                              "(x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) )", # Cohen's d
                                              "(1 - 3 / (4 * (n1 + n2 - 2) - 1) ) * ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ) )", # Hedges' g
                                              "(x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) )", # Cohen's d PAIRED
                                              "(1 - 3 / (4 * (n1 + n2 - 2) - 1) ) * ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ) )", # Hedges' g PAIRED
                                              #
                                              "log(x1 / x2)", # lnRox_paired
                                              "log(((a) * (d)) / ((b) * (c)))", # lnOR
                                              "log( ((a) / n1) / ((c)/ n2) )", # lnRR
                                              "log(sd1 / x1) - log(sd2 / x2)", # lnCVR first order
                                              "log((sd1 / x1) / (sd2 / x2)) + 1/2 * (1 / (n1 - 1) - 1 / (n2 - 1)) + 1/2 * ((sd2^2/(n2 * x2^2)) - (sd1^2 / (n1 * x1^2)))", # lnCVR second order
                                              "log(sd1 / x1) - log(sd2 / x2)", # lnCVR paired first order
                                              "log((sd1 / x1) / (sd2 / x2)) + 1/2 * (1 / (n1 - 1) - 1 / (n2 - 1)) + 1/2 * ((sd2^2/(n2 * x2^2)) - (sd1^2 / (n1 * x1^2)))", #lnCVR paired second order
                                              # WISHART:
                                              "log(sd1 / x1) - log(sd2 / x2)", # lnCVR first order
                                              "log((sd1 / x1) / (sd2 / x2)) + 1/2 * (1 / (n1 - 1) - 1 / (n2 - 1)) + 1/2 * ((sd2^2/(n2 * x2^2)) - (sd1^2 / (n1 * x1^2)))", # lnCVR second order
                                              "log(sd1 / x1) - log(sd2 / x2)", # lnCVR paired first order
                                              "log((sd1 / x1) / (sd2 / x2)) + 1/2 * (1 / (n1 - 1) - 1 / (n2 - 1)) + 1/2 * ((sd2^2/(n2 * x2^2)) - (sd1^2 / (n1 * x1^2)))", #lnCVR paired second order
                                              #
                                              "log((n_Aa/(n_AA + n_Aa + n_aa)) / (2 * (n_AA/(n_AA + n_Aa + n_aa)) * (n_aa/(n_AA + n_Aa + n_aa))))", # lnHWE_A
                                              "log(sqrt(((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))) / ((2 * n1 * n2) / (n1 + n2)))) - log(sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))))", #lnM
                                              "log(sqrt(((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))) / ((2 * n1 * n2) / (n1 + n2)))) - log(sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))))" #lnM_paired
                              ),
                              sampling_variance = c("(sd^2 / n) * (-1 * x ^ -2)^2", # reciprocal
                                                    "sd1^2 / (n1 * x1^2) + sd2^2 / (n2 * x2^2)", #lnRox first
                                                    "sd1^2 / (n1 * x1^2) + sd2^2 / (n2 * x2^2) + 0.5 * ( (sd1^4 / (n1^2 * x1^4)) + (sd2^4 / (n2^2 * x2^4)))", # lnRoM second
                                                    "((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Cohen's d
                                                    "(1 - 3 / (4 * (n1 + n2 - 2) - 1) )^2 * ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Hedges' g
                                                    "((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Cohen's d paired
                                                    "(1 - 3 / (4 * (n1 + n2 - 2) - 1) )^2 * ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Hedges' g paired
                                                    #' *WISHART:*
                                                    "((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Cohen's d
                                                    "(1 - 3 / (4 * (n1 + n2 - 2) - 1) )^2 * ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Hedges' g
                                                    "((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Cohen's d paired
                                                    "(1 - 3 / (4 * (n1 + n2 - 2) - 1) )^2 * ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))", # Hedges' g paired
                                                    #
                                                    "(sd1^2 / (n1 * x1^2)) + (sd2^2 / (n2 * x2^2)) - ((2 * r * sd1 * sd2) / (x1 * x2 * sqrt(n1 * n2)))", # lnRox_paired
                                                    "(1 / (a)) + (1 / (b)) + (1 / (c)) + (1 / (d))", # lnOR
                                                    "(1 - (a / n1)) / a + (1 - (c / n2)) / c",# lnRR
                                                    "sd1^2/(n1 * x1^2) + sd2^2/(n2 * x2^2) + 1/(2*(n1 - 1)) + 1/(2*(n2 - 1))", # lnCVR first order
                                                    "sd1^2/(n1 * x1^2) + sd1^4/(2 * n1^2 * x1^4) + n1/(2*(n1 - 1)^2) + sd2^2/(n2 * x2^2) + sd2^4/(2 * n2^2 * x2^4) + n2/(2*(n2 - 1)^2)", # lnCVR second order
                                                    "sd1^2/(n1 * x1^2) + sd2^2/(n1 * x2^2) - 2*r*sd1*sd2/(n1 * x1 * x2) + 1/(n1 - 1) - r^2/(n1 - 1)", # lnCVR_paired first order
                                                    "sd1^2/(n1 * x1^2) + sd1^4/(2 * n1^2 * x1^4) + sd2^2/(n1 * x2^2) + sd2^4/(2 * n1^2 * x2^4) - 2*r*sd1*sd2/(n1 * x1 * x2) + r^2 * sd1^2 * sd2^2 * (x1^4 + x2^4) / (2 * n1^2 * x1^4 * x2^4) + n1/(n1 - 1)^2 - r^2/(n1 - 1) + r^4 * (sd1^8 + sd2^8) / (2 * (n1 - 1)^2 * sd1^4 * sd2^4)", # lnCVR paired second order
                                                    # WISHART:
                                                    "sd1^2/(n1 * x1^2) + sd2^2/(n2 * x2^2) + 1/(2*(n1 - 1)) + 1/(2*(n2 - 1))", # lnCVR first order
                                                    "sd1^2/(n1 * x1^2) + sd1^4/(2 * n1^2 * x1^4) + n1/(2*(n1 - 1)^2) + sd2^2/(n2 * x2^2) + sd2^4/(2 * n2^2 * x2^4) + n2/(2*(n2 - 1)^2)", # lnCVR second order
                                                    "sd1^2/(n1 * x1^2) + sd2^2/(n1 * x2^2) - 2*r*sd1*sd2/(n1 * x1 * x2) + 1/(n1 - 1) - r^2/(n1 - 1)", # lnCVR_paired first order
                                                    "sd1^2/(n1 * x1^2) + sd1^4/(2 * n1^2 * x1^4) + sd2^2/(n1 * x2^2) + sd2^4/(2 * n1^2 * x2^4) - 2*r*sd1*sd2/(n1 * x1 * x2) + r^2 * sd1^2 * sd2^2 * (x1^4 + x2^4) / (2 * n1^2 * x1^4 * x2^4) + n1/(n1 - 1)^2 - r^2/(n1 - 1) + r^4 * (sd1^8 + sd2^8) / (2 * (n1 - 1)^2 * sd1^4 * sd2^4)", # lnCVR paired second order
                                                    #
                                                    "((1 - (n_Aa/(n_AA + n_Aa + n_aa))) / (4 * (n_AA/(n_AA + n_Aa + n_aa)) * (n_aa/(n_AA + n_Aa + n_aa)) * (n_Aa/(n_AA + n_Aa + n_aa)))) / (n_AA + n_Aa + n_aa)", # lnHWE_A
                                                    "(1 / 4 * ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)))^2) * ( (n0/2)^2 * (2 * (sd1^2 / n1 + sd2^2 / n2)^4 + 4 * (sd1^2 / n1 + sd2^2 / n2)^2 * (x1 - x2)^2) + ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2)^2 / (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))^2) * ((2 * (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))^2) / n1 + n2 - 2) )", #lnM
                                                    "(1 / 4 * ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)))^2) * ( (n1/2)^2 * ( ((2 * (sd1^2 / n1 + sd2^2 / n2)^4)/n1^2) + ((4 * (sd1^2 / n1 + sd2^2 / n2)^2 * (x1 - x2)^2)/n1) ) ) + ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2)^2 / (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))^2) * ((sd1^4 + sd2^4 + 2*r^2*sd1^2*sd2^2) / (2 * (n-1)))" # lnM_paired
                                                    ),
                              sim_family = c( # for now, number of simulated terms + distribution type
                                "1_normal",# "reciprocal",
                                "2_multivariate_normal",# "lnRoM",
                                "2_multivariate_normal",# "lnRoM",
                                "4_multivariate_normal",# "Cohens_d",
                                "4_multivariate_normal",# "Hedges_g",
                                "4_multivariate_normal",# "Cohens_d",
                                "4_multivariate_normal",# "Hedges_g",
                                # WISHART:
                                "4_multivariate_normal_wishart",# "Cohens_d",
                                "4_multivariate_normal_wishart",# "Hedges_g",
                                "4_multivariate_normal_wishart",# "Cohens_d",
                                "4_multivariate_normal_wishart",# "Hedges_g",
                                #
                                "2_multivariate_normal",# "lnRoM_paired",
                                "4_binomial",# "lnOR",
                                "2_binomial",# "lnRR",
                                "4_multivariate_normal",# "lnCVR",
                                "4_multivariate_normal",# "lnCVR",
                                "4_multivariate_normal",# "lnCVR_paired",
                                "4_multivariate_normal",# "lnCVR_paired",
                                #  WISHART:
                                "4_multivariate_normal_wishart",# "lnCVR",
                                "4_multivariate_normal_wishart",# "lnCVR",
                                "4_multivariate_normal_wishart",# "lnCVR_paired",
                                "4_multivariate_normal_wishart",# "lnCVR_paired",
                                #
                                "3_multinomial",# "lnHWE_A"
                                "4_multivariate_normal", # lnM
                                "4_multivariate_normal" # lnM_paired
                              ),
                              vars_required = c(
                                "n, x, sd",# "reciprocal",
                                "x1, x2, sd1, sd2, n1, n2",# "lnRoM",
                                "x1, x2, sd1, sd2, n1, n2",# "lnRoM",
                                "x1, x2, sd1, sd2, n1, n2",# "Cohens_d",
                                "x1, x2, sd1, sd2, n1, n2",# "Hedges_g",
                                "x1, x2, sd1, sd2, r, n1, n2",# "Cohens_d_paired",
                                "x1, x2, sd1, sd2, r, n1, n2",# "Hedges_g_paired",
                                # WISHART:
                                "x1, x2, sd1, sd2, n1, n2",# "Cohens_d",
                                "x1, x2, sd1, sd2, n1, n2",# "Hedges_g",
                                "x1, x2, sd1, sd2, r, n1, n2",# "Cohens_d_paired",
                                "x1, x2, sd1, sd2, r, n1, n2",# "Hedges_g_paired",
                                #
                                "x1, x2, sd1, sd2, r, n1, n2",# "lnRoM_paired",
                                "a, b, c, d",# "lnOR",
                                "a, c, n1, n2",# "lnRR",
                                "x1, x2, sd1, sd2, n1, n2",# "lnCVR",
                                "x1, x2, sd1, sd2, n1, n2",# "lnCVR",
                                "x1, x2, sd1, sd2, r, n1, n2",# "lnCVR_paired",
                                "x1, x2, sd1, sd2, r, n1, n2",# "lnCVR_paired",
                                # WISHART:
                                "x1, x2, sd1, sd2, n1, n2",# "lnCVR",
                                "x1, x2, sd1, sd2, n1, n2",# "lnCVR",
                                "x1, x2, sd1, sd2, r, n1, n2",# "lnCVR_paired",
                                "x1, x2, sd1, sd2, r, n1, n2",# "lnCVR_paired",
                                #
                                "n_AA, n_Aa, n_aa",# "lnHWE_A")
                                "x1, x2, sd1, sd2, n1, n2", # lnM
                                "x1, x2, sd1, sd2, n1, n2, r" # lnM_paired
                              ),
                              cloud_filtering_rules = c(
                                NA,# "reciprocal",
                                "x1 > 0 & x2 > 0",# "lnRoM",
                                "x1 > 0 & x2 > 0",# "lnRoM",

                                "v1 > 0 & v2 > 0",# "Cohens_d",
                                "v1 > 0 & v2 > 0",# "Hedges_g",
                                "v1 > 0 & v2 > 0",# "Cohens_d_paired",
                                "v1 > 0 & v2 > 0",# "Hedges_g_paired",

                                # Wishart:
                                NA, #"v1 > 0 & v2 > 0",# "Cohens_d",
                                NA, #"v1 > 0 & v2 > 0",# "Hedges_g",
                                NA, #"v1 > 0 & v2 > 0",# "Cohens_d_paired",
                                NA, #"v1 > 0 & v2 > 0",# "Hedges_g_paired",
                                #
                                "x1 > 0 & x2 > 0",# "lnRoM_paired",
                                NA,# "lnOR",
                                NA,# "lnRR",
                                "x1 > 0 & x2 > 0 & v1 > 0 & v2 > 0",# "lnCVR",
                                "x1 > 0 & x2 > 0 & v1 > 0 & v2 > 0",# "lnCVR",
                                "x1 > 0 & x2 > 0 & v1 > 0 & v2 > 0",# "lnCVR_paired",
                                "x1 > 0 & x2 > 0 & v1 > 0 & v2 > 0",# "lnCVR_paired",
                                # Wishart:
                                "x1 > 0 & x2 > 0", #"v1 > 0 & v2 > 0",# "lnCVR",
                                "x1 > 0 & x2 > 0", #"v1 > 0 & v2 > 0",# "lnCVR",
                                "x1 > 0 & x2 > 0", #"v1 > 0 & v2 > 0",# "lnCVR_paired",
                                "x1 > 0 & x2 > 0", #"v1 > 0 & v2 > 0",# "lnCVR_paired",
                                #
                                NA, # "lnHWE_A"
                                "v1 > 0 & v2 > 0", # lnM
                                "v1 > 0 & v2 > 0" # lnM

                              ),
                              special_warnings = c(
                                                NA, # "reciprocal",
                                                "lnRoM cannot accept x1 or x2 ≤ 0.", # "lnRoM",
                                                "lnRoM cannot accept x1 or x2 ≤ 0.", # "lnRoM",
                                                #
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD", # Cohen's d
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD",# Hedges' g
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD_paired", # Cohen's d
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD_paired", # Hedges' g
                                                # Wishart:
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD", # Cohen's d
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD",# Hedges' g
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD_paired", # Cohen's d
                                                "SMD can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "SMD_paired", # Hedges' g
                                                #
                                                "lnRoM cannot accept x1 or x2 ≤ 0.", # "lnRoM_paired",
                                                "If any group (a, b, c, d) is 0, lnOR will be incalculable. If this is the case with your data, add 0.5 to all groups and rerun.", # "lnOR",
                                                "If any group (a, c) is 0, lnRR will be incalculable. If this is the case with your data, add 0.5 to all groups and rerun.",# "lnRR",
                                                #
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR",
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR",
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR_paired",
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR_paired",
                                                # Wishart:
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR",
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR",
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR_paired",
                                                "lnCVR can use normal or Wishart distribution to model variances. If unspecified with 'SAFE_distribution' argument, defaulting to Wishart", # "lnCVR_paired",
                                                #
                                                "Cannot accept values of 0 for n_AA, n_Aa, or n_aa", # "lnHWE_A"
                                                "Results are only reliable with SAFE method", #'lnM [Is this true?]
                                                "Results are only reliable with SAFE method" #'lnM [Is this true?]

                              ),
                              default_safe_family = c(
                                          "yes",#' "reciprocal",
                                          "yes",#' "lnRoM",
                                          "yes", #' "lnRoM",
                                          "no",#' "SMD", # Cohen's d
                                          "no",#' "SMD",# Hedges' g
                                          "no",#' "SMD_paired", # Cohen's d
                                          "no",#' "SMD_paired", # Hedges' g
                                          #' #' [Wishart SMD:]
                                          "yes",#' "SMD", # Cohen's d
                                          "yes",#' "SMD",# Hedges' g
                                          "yes",#' "SMD_paired", # Cohen's d
                                          "yes",#' "SMD_paired", # Hedges' g
                                          #' #
                                         "yes", #' "lnRoM_paired",
                                         "yes",#' "lnOR",
                                         "yes", #' "lnRR",
                                         #
                                         "no", #' "lnCVR",
                                         "no", #' "lnCVR",
                                         "no",#' "lnCVR_paired",
                                         "no",#' "lnCVR_paired",
                                          #' #' [Wishart lnCVR:]
                                          "yes",#' "lnCVR",
                                          "yes",#' "lnCVR",
                                          "yes",#' "lnCVR_paired",
                                          "yes",#' "lnCVR_paired",
                                          #' #
                                          "yes",#' "lnHWE_A",
                                          "yes",#' "lnM",
                                          "yes"#' "lnM_paired"
                                          )
                            )

effect_formulas

effect_formulas <- melt(effect_formulas,
                        id.vars = c("name", "derivative", "sim_family", "vars_required",
                                    "cloud_filtering_rules", "special_warnings",
                                    "default_safe_family"),
                        variable.name = "calc_type",
                        value.name = "formula")
effect_formulas <- effect_formulas[!is.na(formula)]


# >>> Filter out the unpaired lnROM and SMD and CVR -----------------------
effect_formulas[name %in% c("SMD", "lnRoM", "lnCVR")]
effect_formulas[name %in% c("SMD_paired", "lnRoM_paired", "lnCVR_paired")]
effect_formulas[grepl("wishart", sim_family)]

effect_formulas <- effect_formulas[!name %in% c("SMD", "lnRoM", "lnCVR")]
effect_formulas[grepl("wishart", sim_family)]

effect_formulas[name %in% c("SMD_paired", "lnRoM_paired", "lnCVR_paired"), accepts_paired_design := "yes"]
effect_formulas[name %in% c("SMD_paired", "lnRoM_paired", "lnCVR_paired"), name := gsub("_paired", "", name)]
effect_formulas[is.na(accepts_paired_design), accepts_paired_design := "no"]

effect_formulas[grepl("wishart", sim_family)]



# >>> Reclassify formulas -------------------------------------------------

effect_formulas[derivative == "first" & calc_type == "effect_size",
                formula_type := "yi_definition"]

effect_formulas[derivative == "second" & calc_type == "effect_size",
                formula_type := "yi_bias_corrected"]

effect_formulas[derivative == "first" & calc_type == "sampling_variance",]
effect_formulas[derivative == "first" & calc_type == "sampling_variance",
                formula_type := "vi_first_order"]

effect_formulas[derivative == "second" & calc_type == "sampling_variance",]
effect_formulas[derivative == "second" & calc_type == "sampling_variance",
                formula_type := "vi_second_order"]

setnames(effect_formulas, c("sim_family", "default_safe_family"),
                          c("SAFE_family", "default_SAFE_family"))

effect_formulas

# >>> Add an execution string with explicit environment/object control ----------------
#' [feel like the previous approach with list2env() is a bit risky. Going to make this explicit.]
# This is actually going to be slightly tricky... Cause we don't want 'a' from ln_OR gsub'ed with n_aa from HWE
# This is risky with 'r'. Because of sqrt....but also if i can figure out the regex on this, it'll make the function safer..

# test <- c("sqrt(r)", "2*r*n1*sqrt(n1)", "r * 2 + n1 * sqrt(n1)", "rnorm(5, 100, 10) * r + 10")
# gsub("\\br\\b", "HAH", test)
# "sqrt(HAH)" "2*HAH*n1*sqrt(n1)" "HAH * 2 + n1 * sqrt(n1)" "rnorm(5, 100, 10) * HAH + 10"

i <- 1
k <- 1
effect_formulas[, exec_formula := formula]

for(i in 1:nrow(effect_formulas)){
  vars <- strsplit(unique(effect_formulas[i, ]$vars_required), split = ", ") |> unlist()
  to <- paste0("input$", vars)
  vars <- paste0("\\b", vars, "\\b")

  for(k in 1:length(vars)){
    effect_formulas[i, exec_formula := gsub(vars[k], to[k], exec_formula)]
  }

}

effect_formulas

effect_formulas[grepl("paired", name) & grepl("sqrt", exec_formula)]

effect_formulas[, source := "Stable Release 2025 v1.1"]

# >>> Save table ----------------------------------------------------------
# Stable_Formula_Backup <- copy(effect_formulas)
# save(Stable_Formula_Backup, file = "data/Stable_Formula_Backup.rda")

usethis::use_data(effect_formulas,
                  overwrite = TRUE)

# Active_Formulas <- copy(effect_formulas)
# save(Active_Formulas, file = "data/Active_Formulas.rda")

# fwrite(effect_formulas, "data/Stable_Formula_Backup.csv", na = "NA")
# fwrite(effect_formulas, "data/Active_Formulas.csv", na = "NA")

unique(effect_formulas$cloud_filtering_rules)

# 2. SCRATCH SPACE FOR WRITING OUT FORMULAS: ------------------------------------------------------

# >>> Reciprocal ----------------------------------------------------------
x    <- mean(c(11.3, 9.7, 10.4, 12.0, 9.1))    # five maze completion times for five insects
sd <-  sd(c(11.3, 9.7, 10.4, 12.0, 9.1))
n    <- length(c(11.3, 9.7, 10.4, 12.0, 9.1))    # number of observations in vector x

# Sampling variance:
(sd^2 / n) * (-1 * x ^ -2)^2



# >>> Cohen's d -----------------------------------------------------------
x1 <- 15.2; sd1 <- 5.3; n1 <- 20 # group 1
x2 <- 12.7; sd2 <- 4.9; n2 <- 18 # group 2
r <- 0.5

# Effect:
# simple
(x1 - x2) / sqrt( Sp )

# Sp
((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)

# Combined
(x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) )

# Variance:
((n1 + n2) / (n1 * n2)) + d^2 / (2 * (n1 + n2 - 2))

# combined:
((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))

# testing against manuscript:
sqrt( ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2)) )
# good

# >>> Hedges' g: ------------------------------------------------------
# J =
(1 - 3 / (4 * (n1 + n2 - 2) - 1) )

# Sp =
(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

# Simple =
J * ((x1 - x2) / sqrt( Sp ) )

# Complete =
(1 - 3 / (4 * (n1 + n2 - 2) - 1) ) * ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ) )

# Hedges' g variance:
# simple =
J^2 * Var_d


# Complete =
(1 - 3 / (4 * (n1 + n2 - 2) - 1) )^2 * ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2))

sqrt((1 - 3 / (4 * (n1 + n2 - 2) - 1) )^2 * ((n1 + n2) / (n1 * n2)) + ((x1 - x2) / sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) ))^2 / (2 * (n1 + n2 - 2)))
# matches!



#' #' *testing with code from ms*
#' x1 <- 15.2;  s1 <- 5.3;  n1 <- 20   # group 1
#' x2 <- 12.7;  s2 <- 4.9;  n2 <- 18   # group 2
#' df <- n1 + n2 - 2
#' ## 1. Plug‐in Cohen d & delta‐method SE
#' Sp2    <- ((n1-1)*s1^2 + (n2-1)*s2^2) / df # This matches what I wrote
#' d_hat  <- (x1 - x2) / sqrt(Sp2) # This matches what I wrote
#' Var_d  <- (n1+n2)/(n1*n2) + d_hat^2/(2*df) # This matches what I wrote
#' se_delta <- sqrt(Var_d) # This matches as well
#'
#' ## 2. Hedges' g & its SE --------------------------------------------------
#' J      <- 1 - 3/(4*df - 1); J # This matches
#' g_hat  <- J * d_hat; g_hat # This matches
#' Var_g  <- J^2 * Var_d; Var_g
#' se_g   <- sqrt(Var_g); se_g

# >>> Paired lnRoM ------------------------------------------------------

# variance only:
(sd1^2 / (n1 * x1^2)) + (sd2^2 / (n2 * x2^2)) - ((2 * r * sd1 * sd2) / (x1 * x2 * sqrt(n1 * n2)))

# >>> lnOR ----------------------------------------------------------------

a <- 10
b <- 20
c  <- 20
d <- 5

# Effect:
log(((a+0.5) * (d+0.5)) / ((b+0.5) * (c+0.5)))

# Variance:
(1 / (a+0.5)) + (1 / (b+0.5)) + (1 / (c+0.5)) + (1 / (d+0.5))


# >>> lnRR ----------------------------------------------------------------
a = 2; b = 20; c = 7; d = 13
n1 <- a + b
n2 <- c + d

# Effect:
log( ((a+0.5) / n1) / ((c + 0.5)/ n2) )

# Variance:
# Simple:
(1 - p1) / a + (1 - p2) / c

#p1:
a / n1

#p2:
c / n2

# Full:
(1 - (a / n1)) / a + (1 - (c / n2)) / c


# >>> lnCVR ---------------------------------------------------------------
x1 = 17; sd1 = 2; n1 = 23;
x2 =  12; sd2 = 3; n2 = 27;
r  = 0

# first order:
# Effect:
log(sd1 / x1) - log(sd2 / x2)

# Unpaired:
# Variance:
# First:
(sd1^2 / (n1 * x1^2)) + (sd2^2 / (n2 * x2^2)) + 1/(2 * (n1 - 1)) + 1 / (2*(n2 - 1))

# Second order:
# Effect:
sd1^2/(n1 * x1^2) + sd1^4/(2 * n1^2 * x1^4) + n1/(2*(n1 - 1)^2) + sd2^2/(n2 * x2^2) + sd2^4/(2 * n2^2 * x2^4) + n2/(2*(n2 - 1)^2)



# >>> lnCVR with r --------------------------------------------------------

# Variance
#First
sd1^2/(n1 * x1^2) + sd2^2/(n1 * x2^2) - 2*r*sd1*sd2/(n1 * x1 * x2) + 1/(n1 - 1) - r^2/(n1 - 1)

# Second
sd1^2/(n1 * x1^2) + sd1^4/(2 * n1^2 * x1^4) + sd2^2/(n1 * x2^2) + sd2^4/(2 * n1^2 * x2^4) - 2*r*sd1*sd2/(n1 * x1 * x2) + r^2 * sd1^2 * sd2^2 * (x1^4 + x2^4) / (2 * n1^2 * x1^4 * x2^4) + n1/(n1 - 1)^2 - r^2/(n1 - 1) + r^4 * (sd1^8 + sd2^8) / (2 * (n1 - 1)^2 * sd1^4 * sd2^4)


# >>> Hardy Weinberg ------------------------------------------------------
n_AA <- 40
n_Aa <- 25
n_aa <- 50

# I'x guessing you'd want one for A and one for a?
# Effect
# n =
(n_AA + n_Aa + n_aa)

# p1
n_AA/n
(n_AA/(n_AA + n_Aa + n_aa))

#p2
n_Aa/n
(n_Aa/(n_AA + n_Aa + n_aa))

#p3
n_aa/n
(n_aa/(n_AA + n_Aa + n_aa))

# Short:
log(p2 / (2 * p1 * p3))

# Coxplete:
log((n_Aa/(n_AA + n_Aa + n_aa)) / (2 * (n_AA/(n_AA + n_Aa + n_aa)) * (n_aa/(n_AA + n_Aa + n_aa))))

# Variance:
# Short:
((1 - p2) / (4 * p1 * p3 * p2)) / n

# Complete
((1 - (n_Aa/(n_AA + n_Aa + n_aa))) / (4 * (n_AA/(n_AA + n_Aa + n_aa)) * (n_aa/(n_AA + n_Aa + n_aa)) * (n_Aa/(n_AA + n_Aa + n_aa)))) / (n_AA + n_Aa + n_aa)

sqrt(((1 - (n_Aa/(n_AA + n_Aa + n_aa))) / (4 * (n_AA/(n_AA + n_Aa + n_aa)) * (n_aa/(n_AA + n_Aa + n_aa)) * (n_Aa/(n_AA + n_Aa + n_aa)))) / (n_AA + n_Aa + n_aa))



# >>> lnM -----------------------------------------------------------------
x1 <- 15.2; sd1 <- 5.3; n1 <- 20 # group 1
x2 <- 12.7; sd2 <- 4.9; n2 <- 18 # group 2
r <- 0.5


# EFFECT:
# MSb
(((n1 * n2) / (n1 + n2)) * (x1 - x2)^2)


# MSw == s2w
(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))


# n0
((2 * n1 * n2) / (n1 + n2))

# s2b simple
(MSb - MSw) / n0


# s2b full
((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))) / ((2 * n1 * n2) / (n1 + n2))


# Simple
log(sqrt(s2b)) - log(sqrt(s2w))

# Full
log(sqrt(((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))) / ((2 * n1 * n2) / (n1 + n2)))) - log(sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))))

# Variance:

# Simple:
(1 / 4 * delta^2) * ( (n0/2)^2 * (2 * sdD^4 + 4 * sdD^2 * δ^2) + (MSb^2 / MSw^2) * ((2 * MSw^2) / n1 + n2 - 2) )

# delta:
MSb - MSw
((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)))

# n0
((2 * n1 * n2) / (n1 + n2))

# δ
(x1 - x2)

# s2D
(sd1^2 / n1 + sd2^2 / n2)

# Complete:
(1 / 4 * ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)))^2) * ( (n0/2)^2 * (2 * (sd1^2 / n1 + sd2^2 / n2)^4 + 4 * (sd1^2 / n1 + sd2^2 / n2)^2 * (x1 - x2)^2) + ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2)^2 / (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))^2) * ((2 * (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))^2) / n1 + n2 - 2) )


# PAIRED variance, simple:
(1 / 4 * delta^2) * ( (n1/2)^2 * ( ((2 * sdD^4)/n1^2) + ((4 * sdD^2 * δ^2)/n1) ) ) +
                        (MSb^2 / MSw^2) *
                        ((sd1^4 + sd2^4 + 2*r^2*sd1^2*sd2^2) / (2 * (n-1)))


(1 / 4 * ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2) - (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)))^2) * ( (n1/2)^2 * ( ((2 * (sd1^2 / n1 + sd2^2 / n2)^4)/n1^2) + ((4 * (sd1^2 / n1 + sd2^2 / n2)^2 * (x1 - x2)^2)/n1) ) ) + ((((n1 * n2) / (n1 + n2)) * (x1 - x2)^2)^2 / (((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))^2) * ((sd1^4 + sd2^4 + 2*r^2*sd1^2*sd2^2) / (2 * (n-1)))


