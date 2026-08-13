
rm(list = ls())

# source('scripts/DEV_one_ring_to_bind_them/1_universal_SAFE.R')

library("data.table")
library("MASS")

#' [FUNCTION DESCRIPTION]
# This function will return point and variance estimates based on plugins and SAFE method
# It's been designed so that users can add new effect sizes (by adding to the formula csv that underlies calculations)
# It also can accept a custom sigma_matrix argument for these new effect sizes, if they do not
# fall into already developed categories

input_vars = list(x1 = c(15, 8, 14, 30), # means
                  x2 = c(2, 5, 15, 3),
                  sd1 = c(1, 2, 0.5, 1.2), # sds
                  sd2 = c(0.7, 1, .75, 1),
                  n1 = c(15, 10, 10, 13), # ns
                  n2 = c(14, 11, 9, 12))
effect_type = "lnRoM" # designate effect size
paired_design <- FALSE
SAFE = TRUE
SAFE_boots = 1e6
SAFE_distribution = NULL; sigma_matrix = NULL; verbose = T; use_custom_formulas = FALSE

# FOR DEBUGGING SAFE_calc:
k <- 1
formulas = effect_formulas.sub
input_k = lapply(input_vars, "[[", k) # select the first element in each element...
plugin_effect_k = plugin_effect_size[k]
sigma_matrix_k = sigma_matrix[[k]] # submit custom sigma_matrix if it exists.
SAFE_boots = 1e6

# For parameter cloud:



#' [Quick examples:]
eff_size(x1 = c(15, 8, 14, 30), # means
         x2 = c(2, 5, 15, 3),
         sd1 = c(1, 2, 0.5, 1.2), # sds
         sd2 = c(0.7, 1, .75, 1),
         n1 = c(15, 10, 10, 13), # ns
         n2 = c(14, 11, 9, 12),
         effect_type = "lnRoM", # designate effect size
         SAFE = TRUE,
         SAFE_boots = 1e6)

eff_size(x1 = c(15, 8, 14, 30), # means
         x2 = c(2, 5, 15, 3),
         sd1 = c(1, 2, 0.5, 1.2), # sds
         sd2 = c(0.7, 1, .75, 1),
         n1 = c(15, 10, 10, 13), # ns
         n2 = c(14, 11, 9, 12),
         effect_type = "SMD", # designate effect size
         SAFE = TRUE,
         SAFE_boots = 1e6)

#' [The function has some error checking built in]
eff_size()


eff_size(effect_type = "lnRoM",
         x1 = 5,
         x2 = c(5, 10),
         n1 = c(4),
         n2 = c(4),
         sd1 = c(0.5),
         sd2 = c(0.1))

eff_size(effect_type = "baloney!")

#' [Adding new effect sizes simply requires adding a row to this csv with no other changes to underlying function code]
#' As long as the new effect size uses one of the already implemented sigma matrices.
#'
#' [Formulas are here:]
guide <- fread("data/effect_size_formulas.csv")
guide

#' [NOTE:] The names/labels will be simplified and there are some slight inaccuracies:
#' Cohen's d is being described as a first order effect and Hedges' g a second order,
#' even though these aren't derived from the delta method as far as I know. We can change terminology

#' [TO DO:]
#' *1.* Create functionality to add effect sizes.
#' *2.* What about zcor?
#' *3.* Diagnostic functions based on CV. :)


#  TEST OF EFF_SIZE FUNCTION -------------------------------------

# 1. Speed -------------------------------------------------------------------

x <- mean(c(11.3, 9.7, 10.4, 12.0, 9.1))
sd <- sd(c(11.3, 9.7, 10.4, 12.0, 9.1))
n <- length(c(11.3, 9.7, 10.4, 12.0, 9.1))

eff_size(x = x, sd = sd, n = n,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "reciprocal")

#' [SAFE is correct but plugin SE/Var is not...]

# This matches the correct results:
#' [Correct results:]
# plugin point:  plugin SE   plugin var      SE_SAFE     SAFE var
# 0.0952        0.0950      0.009025        0.0048      2.304e-05


# 2. lnRoM ----------------------------------------------------------------
set.seed(1)
x1 <- 13.4; sd1 <- 4.6; n1 <- 18 # group 1
x2 <- 16.1; sd2 <- 3.9; n2 <- 17 # group 2
source('scripts/DEV_one_ring_to_bind_them/1_universal_SAFE_v3.R')

out <- eff_size(x1 = x1, sd1 = sd1, n1 = n1,
         x2 = x2, sd2 = sd2, n2 = n2,
         r = 0,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "lnRoM")
out
#         First  Second SAFE_BC
# Point -0.1836 -0.1820 -0.1820
# SE     0.1000  0.1001  0.1007
# Var    0.01    0.010   0.0101
#' [MATCHES!]

# >>> paired ----------------------------------------------------------
# Now paired;
source('scripts/DEV_one_ring_to_bind_them/1_universal_SAFE_v3.R')
# x1 <- 13.4; sd1 <- 4.6; n1 <- 18 # group 1
# x2 <- 16.1; sd2 <- 3.9; n2 <- 17 # group 2
# r = 0.5
# input_vars <- list(x1=x1, x2=x2, sd1=sd1, sd2=sd2, n1=n1, n2=n2, r=r)
out <- eff_size(x1 = 13.4, sd1 = 4.6, n1 = 18,
                x2 = 16.1, sd2 = 3.9, n2 = 17,
                r = 0.5,
                SAFE = TRUE,
                SAFE_boots = 1e6,
                verbose = T,
                effect_type = "lnRoM_paired")
out

# 3. Hedges' g / Cohen d ------------------------------------------------------------

x1 <- 15.2;  sd1 <- 5.3;  n1 <- 20   # group 1
x2 <- 12.7;  sd2 <- 4.9;  n2 <- 18   # group 2

out <- eff_size(x1 = x1, sd1 = sd1, n1 = n1,
         x2 = x2, sd2 = sd2, n2 = n2,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "SMD")
out
#' *Second order refers to g, first order to d*
#' [MATCHES]
#       d         g         d_BC
# Point 0.4888    0.4785    0.4776
# SE    0.3300    0.3230    0.3408
# Var   0.1089    0.1043    0.1161


# 4. lnOR -----------------------------------------------------------------
a <- 2 ; b <- 20            #  2 / 22 events in Group 1
c <- 10; d <- 12            # 10 / 22 events in Group 2

out <- eff_size(a = a, b = b, c = c,
         d = d,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "lnOR")
out
#' [This matches. :)]
#       Plug_in       SAFE_BC
# Point -2.1203       -1.9515
# SE     0.8563       0.8714
# Var    0.7332       0.7593


# 5. lnRR -----------------------------------------------------------------

# a = 2; b = 20; c = 7; d = 13
# n1 <- a + b
# n2 <- c + d

out <- eff_size(a = 2, n1 = 22,
         c = 7, n2 = 20,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "lnRR")
out

#' [This matches :)]
#        Plug_in    SAFE_BC
# Point -1.3481     -1.2167
# SE     0.7399     0.7690
# Var    0.5475     0.5913


# 6. lnCVR -------------------------------------------------------------------
# x1 = 17, s1 = 2, n1 = 23,
# x2 =  12, s2 = 3, n2 = 27,
# r  = 0

out <- eff_size(x1 = 17, sd1 = 2, n1 = 23,
         x2 =  12, sd2 = 3, n2 = 27,
         r  = 0,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "lnCVR")
out
#' [Matches]
#         First       Second    SAFE_BC
# Point   -0.7538     -0.7494     -0.7479
# SE      0.2118      0.2160      0.2453
# Var     0.04485     0.0466      0.0601



# 6.5 lnCVR paired --------------------------------------------------------
out <- eff_size(x1 = 15, sd1 = 2, n1 = 25,
                x2 = 10, sd2 = 2, n2 = 25,
                r  = 0.5,
                SAFE = TRUE,
                SAFE_boots = 1e6,
                verbose = T,
                effect_type = "lnCVR_paired")
out

#' [INCORRECT]
#           First       Second      SAFE_BC
# Point     -0.4055     -0.4050     -0.4052
# SE        0.1803      0.1853      0.2112
# Var       0.0325      0.0343      0.0446


# 7. HWE ------------------------------------------------------------------
n_AA <- 40; n_Aa <- 25; n_aa <- 50

eff_size(n_AA = 40, n_Aa = 25, n_aa = 50,
         SAFE = TRUE,
         SAFE_boots = 1e6,
         verbose = T,
         effect_type = "lnHWE_A")
#' [That looks correct]

#       Delta         SAFE_BC
# Point -0.3302       -0.3279
# SE     0.2275       0.2824
# Var    0.0517       0.0797


# 8. Just to demonstrate, this works with vectors of effect sizes as well: ----------------------------------------

#' [Let's see how much time it adds to do SAFE calc for 4 effect sizes:]

# No SAFE calculation
s <- Sys.time()
eff_size(x1 = c(15, 8, 14, 30),
         x2 = c(2, 5, 15, 3),
         sd1 = c(1, 2, 0.5, 1.2),
         sd2 = c(0.7, 1, .75, 1),
         n1 = c(15, 10, 10, 13),
         n2 = c(14, 11, 9, 12),
         effect_type = "lnRoM",
         derivative_order = "first",
         SAFE = FALSE,
         SAFE_boots = 1e6)
e <- Sys.time()
e - s # 0.008823872 Secs

# with SAFE bootstrapping:
s <- Sys.time()
eff_size(x1 = c(15, 8, 14, 30),
         x2 = c(2, 5, 15, 3),
         sd1 = c(1, 2, 0.5, 1.2),
         sd2 = c(0.7, 1, .75, 1),
         n1 = c(15, 10, 10, 13),
         n2 = c(14, 11, 9, 12),
         effect_type = "lnRoM",
         derivative_order = "first",
         SAFE = TRUE,
         SAFE_boots = 1e6)
e <- Sys.time()
e - s # 1.038215 secs.

#' [Should probably include a parallel option...]
