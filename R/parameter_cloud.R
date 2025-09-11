

# So many while loops. Lots of duplicated code.
# I think I should encapsulate the draws. so they can be in 1 while loop

# Is it even possible to make the SAFE function generalizable?
# I guess if the custom formulas fit into the templates...x1, x2, sd1, sd2, a, b, c, d...



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -------------------------------------
# Compile means and sigma matrices. If sigma_matrix is specified then only the means are compiled...
# Could maybe be tidier.Thinking about to more broadly generalize the SAFE call...
.determine_sigma_matrix <- function(family,
                                    sigma_matrix = NULL,
                                    input){
  if(family == "1_normal"){
    if(is.null(sigma_matrix)){
      sigma_matrix <- input$sd / sqrt(input$n)
    }
    means <- c(x = input$x)

  }else if(family %in% c("2_multivariate_normal")){
    if(is.null(sigma_matrix)){

      sigma_matrix <- matrix(data = c((input$sd1^2 / input$n1),                    (input$r*input$sd1*input$sd2)/input$n1, #  / n1 add this to sd1^2
                                      (input$r*input$sd1*input$sd2)/input$n1,      (input$sd2^2 / input$n2)), #  / n2 add this to sd2^2
                             nrow = 2, ncol = 2)

    }
    means <- c(x1 = input$x1, x2 = input$x2)

  }else if(family == "4_multivariate_normal_wishart"){
    # So this is complicated...The wishart call is for the other 2 and has a sigma_matrix / n....confusing as hell
    if(is.null(sigma_matrix)){

      sigma_matrix <- matrix(c(input$sd1^2, input$r*input$sd1*input$sd2,
                               input$r*input$sd1*input$sd2, input$sd2^2),
                             2, 2)

    }
    means <- c(x1 = input$x1, x2 = input$x2)

    # means <- c(x1 = input$x1, x2 = input$x2)
  }else if(family == "4_multivariate_normal"){
    if(is.null(sigma_matrix)){

      sigma_matrix <- matrix(data = c(input$sd1^2/input$n1,                  (input$r*input$sd1*input$sd2)/input$n1, 0,                                                  0,
                                      (input$r*input$sd1*input$sd2)/input$n1, input$sd2^2/input$n2,                   0,                                                  0,
                                      0,                                      0,                                      (2*input$sd1^4)/(input$n1-1),                       ((2*input$r^2*input$sd1^2*input$sd2^2)/(input$n1-1)),
                                      0,                                      0,                                      (2*input$r^2*input$sd1^2*input$sd2^2)/(input$n1-1), (2*input$sd2^4)/(input$n2-1)),
                             nrow = 4,
                             ncol = 4)

    }
    means <- c(x1 = input$x1, x2 = input$x2,
               v1 = input$sd1^2, v2 = input$sd2^2)
  }
  # # Other sigma_matrix options:
  # if(family %in% c("2_multivariate_normal",
  #                  "2_multivariate_lognormal",
  #                  "2_multivariate_Gamma",
  #                  "2_multivariate_negative_binomial",
  #                  "2_multivariate_Poisson",
  #                  "2_multivariate_Beta")){
  #   means <- c(input$x1, input$x2)
  #   sigma_matrix <- c(input$sd1, input$sd2)
  #
  # }
  return(list(means, sigma_matrix))

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -------------------------------------
# This will do the actual draws based on family. The parameter_cloud function will coordinate while loops...
.draw_hyperparameters <- function(means_and_sigma = NULL,
                                  input = NULL,
                                  family){

  if(exists("means_and_sigma")){
    means <- means_and_sigma[[1]]
    sigma_matrix <- means_and_sigma[[2]]
  }

  # NORMAL
  if(family %in% c("1_normal")){
    out <- data.table::data.table(x = rnorm(n=SAFE_boots,
                                            mean = means,
                                            sd = sigma_matrix))
    return(out)

  }
  if(family %in% c("4_multivariate_normal",
                   "2_multivariate_normal")){
    out <- MASS::mvrnorm(n = SAFE_boots,
                         mu = means,
                         Sigma = sigma_matrix) |>
      as.data.frame()

    return(data.table::setDT(out))
  }
  if(family %in% c("4_multivariate_normal_wishart")){
    out <- MASS::mvrnorm(n = SAFE_boots,
                         mu = means,
                         Sigma = sigma_matrix (sigma_matrix / c(input$n1, sqrt(input$n1*input$n2), sqrt(input$n1*input$n2), input$n2))) |>
      as.data.frame()
    setDT(out)

    wishart.out <- stats::rWishart(SAFE_boots,
                                   df = (input$n1-1),
                                   Sigma = sigma_matrix)
    out[, sd1 := sqrt(wishart.out[1, 1, ] / (input$n1 - 1))]
    out[, sd2 := sqrt(wishart.out[2, 2, ] / (input$n2 - 1))]
    out

    return(out)
  }
  # Count data --------------------------------------------------------------
  if(sim_family == "2_binomial"){ # lnRR

    out <- data.table(a = rbinom(SAFE_boots, input$n1, input$a / input$n1) |> as.double(),
                      c = rbinom(SAFE_boots, input$n2, input$c / input$n2) |> as.double())
    out[, n1 := input$n1]
    out[, n2 := input$n2]

    # for lnRR, add 0.5 to affected group and 1 to the appropriate n
    out[a == 0, `:=` (a = a + 0.5,
                      n1 = n1 + 1) ]
    out[c == 0, `:=` (c = c + 0.5,
                      n2 = n2 + 1) ]
    return(out)

  }else if(sim_family == "3_multinomial"){
    N <- (input$n_AA + input$n_Aa + input$n_aa)
    out <- stats::rmultinom(n = SAFE_boots,
                            size = N,
                            prob = c(n_AA = input$n_AA/N,
                                     n_Aa = input$n_Aa/N,
                                     n_aa = input$n_aa/N)) |>
      t() |> # For some reason these are returned WIDE, with 3 rows and 1e6 columns. Weird. Was freezing computer
      as.data.frame()

    data.table::setDT(out)
    out[, `:=` (n_AA = as.double(n_AA),
                n_Aa = as.double(n_Aa),
                n_aa = as.double(n_aa))]

    # Add 0.5 to all groups
    if(nrow(out[(n_AA == 0 | n_Aa == 0 | n_aa == 0), ]) > 0){
      out[, `:=` (n_AA = n_AA + 0.5,
                  n_Aa = n_Aa + 0.5,
                  n_aa = n_aa + 0.5)]
    }

    return(out)
  }else if(sim_family == "4_binomial"){ # this is lnOR

    if(!all(c("n1", "n2") %in% names(input))){
      input$n1 <- input$a + input$b
      input$n2 <- input$c + input$d
    }
    out <- data.table(a = rbinom(SAFE_boots, input$n1, input$a / input$n1) |> as.double(),
                      b = rbinom(SAFE_boots, input$n1, input$b / input$n1) |> as.double(),
                      c = rbinom(SAFE_boots, input$n2, input$c / input$n2) |> as.double(),
                      d = rbinom(SAFE_boots, input$n2, input$d / input$n2) |> as.double())

    # Add 0.5 to rows with ANY zero
    if(nrow(out[(a == 0 | b == 0 | c == 0 | d == 0), ]) > 0){
      out[, `:=` (a = a + 0.5,
                  b = b + 0.5,
                  c = c + 0.5,
                  d = d + 0.5)]
    }
    return(out)
  }

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -------------------------------------
#' Calculate effect sizes
#'
#' Used internally
#'
#' @import data.table
#' @import crayon
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rbinom rWishart rmultinom
#' @param family blah
#' @param input blah
#' @param cloud_filtering_rule blah
#' @param sigma_matrix blah
#' @param SAFE_boots blah
#' @return A data.table with effect sizes and sample variances
#' @export
parameter_cloud <- function(family,
                            input,
                            cloud_filtering_rule = NA,
                            sigma_matrix = NULL,
                            SAFE_boots = 1e6){
  # Determine sigma matrix and means:
  if(family %in% c("4_multivariate_normal",
                   "4_multivariate_normal_wishart",
                   "2_multivariate_normal",
                   "1_normal",
                   "2_multivariate_normal"#,
                   # "2_multivariate_lognormal",
                   # "2_multivariate_Gamma",
                   # "2_multivariate_negative_binomial",
                   # "2_multivariate_Poisson",
                   # "2_multivariate_Beta"
                   )){

    mean_and_sigma <- .determine_sigma_matrix(family = family,
                                              sigma_matrix,
                                              input)
  }
    # Now, draw samples:
    n_boots <- 0
    i <- 1
    # For a timeoutbackup exist option:
    s <- Sys.time()
    cloud <- list()

    # Do all in a while loop although most won't do more than 1 iteration (no cloud filtering rules)
    while(n_boots < SAFE_boots &
          as.numeric(Sys.time() - s) < 60){

      cloud[[i]] <- .draw_hyperparameters(mean_and_sigma = ifelse(exists("mean_and_sigma"),
                                                         mean_and_sigma, NULL),
                                          input = input,
                                          family = family)

      if(!is.na(cloud_filtering_rule)){
        cloud[[i]] <- cloud[[i]][eval(parse(text = cloud_filtering_rule)), ]
      }

      n_boots <- n_boots + nrow(cloud[[i]])
      i <- i + 1
    }
    cloud <- data.table::rbindlist(cloud)
    cloud <- cloud[1:SAFE_boots, ] # Trim to SAFE_boots
    return(cloud)
}

