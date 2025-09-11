
#' Diagnose raw data for lnRoM, for suggestions as to ideal parent family and thus effect size
#'
#' Blah blah blah
#'
#' @import data.table
#' @import crayon
#' @param effect_family Probably unnecessary unless we expand beyond lnRoM
#' @param x1 Numeric vector of means of group 1
#' @param x2 Numeric vector of means of group 2
#' @param sd1 Numeric vector of SDs of group 1
#' @param sd2 Numeric vector of SDs of group 2
#' @param n1 Numeric vector of sample sizes for group 1
#' @param n2 Numeric vector of sample sizes for group 2
#' @param data_type Vector of data types, if known. Accepted types include: "proportion", "counts", "continuous" or NA. If NA, data type is treated as continuous and unbounded.
#' @return A vector of recommended effect sizes within the lnRoM family
#' @export
diagnose_effects <- function(effect_family = "lnRoM",
                             x1, x2, sd1, sd2, n1, n2,
                             t_low = 0.20,
                             t_high = 0.80,
                             delta = 0.05,
                             data_type = NULL){

  # TODO: what do we do if x1/sd1 are Poisson but x2/sd2 are not?
  if(unique(lengths(list(x1, x2, sd1, sd2, n1, n2))) > 1) return(cat(crayon::red("Vectors of unequal lengths")))

  # Create a data.table from vectors to enable fast vector operations instead of loop/lapply
  dat <- data.table(x1, x2, sd1, sd2, n1, n2)

  if(is.null(data_type)){
    dat[, data_type :=  rep("continuous", length(x1))]
  }else{
    dat <- data.table(dat, data_type)
    dat[is.na(data_type), data_type := "continuous"]
  }

  # ---- Initial checks ------------------------------------
  if(any(x1 <= 0) | any(x2 <= 0)) return(cat(crayon::red("Values less than or equal to 0. lnRoM cannot accept 0 or negative values.
                                                         Please add small constant to x1 and x2 for the affected observations and rerun.")))

  # ---- Proportion ------------------------------------
  dat[data_type == "proportion", recommended_eff := "lnRoM_proportion"]


  # ---- Counts ------------------------------------
  # Poisson:
  tolerance <- 10e-3
  # TODO: the notation looks like tolerance is multiplied by max(1, x1). Is that true?
  # Or is it used like the delta statistic in the CV2 section?
  dat[data_type == "counts" &
        abs(sd1^2 - x1) <= 10e-3 * max(1, x1),
      x1_rec := "lnRoM_Poisson"]

  dat[data_type == "counts" &
        abs(sd2^2 - x2) <= 10e-3 * max(1, x2),
      x2_rec := "lnRoM_Poisson"]

  # sd^2 > x:
  # TODO: Tolerance for this comparison too?
  # If so: sd1_x_delta := sd1^2 - x1; sd1_x_delta > 0 +- tolerance
  dat[data_type == "counts" &
        sd1^2 > x1 &
        is.na(x1_rec),
      rho_NB1 := (sd1^2 - x1)/x1]
  dat[data_type == "counts" &
        sd2^2 > x2 &
        is.na(x1_rec),
      rho_NB2 := (sd2^2 - x1)/x2]

  dat[data_type == "counts" &
        sd1^2 > x1 &
        is.na(x1_rec),
     x1_rec := fcase(rho_NB1 <= 0.2, "lnRoM_Poisson",
                     rho_NB1 >= 0.2 & rho_NB1 <= 1, "lnRoM_negative_binomial",
                     rho_NB1 >= 1, "lnRoM_Gamma")]

  dat[data_type == "counts" &
        sd2^2 > x2 &
        is.na(x1_rec),
      x2_rec := fcase(rho_NB2 <= 0.2, "lnRoM_Poisson",
                      rho_NB2 >= 0.2 & rho_NB2 <= 1, "lnRoM_negative_binomial",
                      rho_NB2 >= 1, "lnRoM_Gamma")]

  # ------ Continuous, non-bounded, or unknown data types: -------------------------------
  dat[data_type %in% c("continuous"), CV1_2 := (sd1^2 / x1^2)]
  dat[data_type %in% c("continuous"), CV2_2 := (sd2^2 / x2^2)]

  # TODO. What's teh equation for normal vs lognormal SEs?
  # equivalently, compute both Normal–
  # Delta and Log-Normal SEs and use the larger one (conservative). To O(n−1) the biases
  # coincide at s2/(2nX2).

  dat[data_type %in% c("continuous"), CV1_2_tlow_min := CV1_2 - t_low - delta]
  dat[data_type %in% c("continuous"), CV1_2_tlow_max := CV1_2 - t_low + delta]
  dat[data_type %in% c("continuous"), CV1_2_thigh_min := CV1_2 - t_high - delta]
  dat[data_type %in% c("continuous"), CV1_2_thigh_max := CV1_2 - t_high + delta]

  # TODO Need to double check high versus low thresholds. Much easier to code without thinking about that....
  dat[data_type %in% c("continuous") &
        n1 >= 30,
      x1_rec := fcase(# Low dispersion:
        CV1_2_tlow_min < 0 & CV1_2_tlow_max < 0, "lnRoM",
        CV1_2_tlow_min < 0 & CV1_2_tlow_max > 0, "lnRoM delta ambiguous",
        # Moderate dispersion:
        CV1_2_tlow_min > 0 & CV1_2_tlow_max > 0 & CV1_2_thigh_min < 0 & CV1_2_thigh_max < 0, "lnRoM_lognormal",
        (CV1_2_tlow_min > 0 & CV1_2_tlow_max < 0) | (CV1_2_thigh_min > 0 & CV1_2_thigh_max < 0), "lnRoM_lognormal delta ambiguous",
        # High dispersion:
        CV1_2_thigh_min > 0 & CV1_2_thigh_max > 0, "lnRoM_Gamma",
        CV1_2_thigh_min > 0 & CV1_2_thigh_max < 0, "lnRoM_Gamma delta ambiguous")]

  #
  dat[data_type %in% c("continuous") &
        n2 >= 30,
      x2_rec := fcase(# Low dispersion:
        CV2_2_tlow_min < 0 & CV2_2_tlow_max < 0, "lnRoM",
        CV2_2_tlow_min < 0 & CV2_2_tlow_max > 0, "lnRoM delta ambiguous",
        # Moderate dispersion:
        CV2_2_tlow_min > 0 & CV2_2_tlow_max > 0 & CV2_2_thigh_min < 0 & CV2_2_thigh_max < 0, "lnRoM_lognormal",
        (CV1_2_tlow_min > 0 & CV1_2_tlow_max < 0) | (CV1_2_thigh_min > 0 & CV1_2_thigh_max < 0), "lnRoM_lognormal delta ambiguous",
        # High dispersion:
        CV2_2_thigh_min > 0 & CV2_2_thigh_max > 0, "lnRoM_Gamma",
        CV2_2_thigh_min > 0 & CV2_2_thigh_max < 0, "lnRoM_Gamma delta ambiguous")]

  # Synthesize recommendations -------------------------------------------------
  # INCORPORATE OTHER TYPES OF AMBIGUITY HERE:
  # TODO
  dat[, recommended_eff := fcase(x1_rec == x2_rec, x1_rec,
                                 x1_rec != x2_rec, paste("x1=", x1_rec, "; x2=", x2_rec))]



}
#

# Create sample data:
