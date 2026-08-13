dat <- data.table(x1 = c(15, 20), x2 = c(14.5, 15),
                  sd1 = c(1.5, 2), sd2 = c(1.7, 1.3),
                  n1 = c(15, 20), n2 = c(15, 20))
input_vars <- as.list(dat)

effect_type = "lnM"
data = dat
bind = TRUE
paired = FALSE
default_formulas = TRUE
SAFE = TRUE
SAFE_boots = 1e6
SAFE_max_secs = 15
n_cores = 1
SAFE_distribution = NULL
sigma_matrix = NULL
verbose = T
load("data/effect_formulas.rda")


formulas
input
plugin_effect
custom_sigma
SAFE_boots = 1e6
SAFE_max_secs
