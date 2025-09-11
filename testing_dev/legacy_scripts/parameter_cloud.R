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

  # Construct Gaussian sigma matrices ------------------------------------------------
  if(family == "1_normal"){
    if(is.null(sigma_matrix)){
      sigma_matrix <- input$sd / sqrt(input$n)
    }
    means <- c(x = input$x)

  }else if(family %in% c("2_multivariate_normal",
                         "2_multivariate_lognormal",
                         "2_multivariate_Gamma",
                         "2_multivariate_negative_binomial",
                         "2_multivariate_Poisson",
                         "2_multivariate_Beta")){
    if(is.null(sigma_matrix)){

      sigma_matrix <- matrix(data = c((input$sd1^2 / input$n1),                    (input$r*input$sd1*input$sd2)/input$n1, #  / n1 add this to sd1^2
                                      (input$r*input$sd1*input$sd2)/input$n1,      (input$sd2^2 / input$n2)), #  / n2 add this to sd2^2
                             nrow = 2, ncol = 2)

    }
    means <- c(x1 = input$x1, x2 = input$x2)

  }else if(family == "4_multivariate_normal_wishart"){
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----------------------------------------
  # Create Gaussian clouds ------------------------------------------------------------
  if(family == "1_normal"){

    out <- data.table::data.table(x = rnorm(n=SAFE_boots,
                                mean = means,
                                sd = sigma_matrix))
    return(out)

  }else if(family %in% c("4_multivariate_normal",
                             "2_multivariate_normal")){

    out <- list()

    n_boots <- 0
    i <- 1
    s <- Sys.time()
    while(n_boots < SAFE_boots &
          as.numeric(Sys.time() - s) < 60){

      out[[i]] <- MASS::mvrnorm(n = SAFE_boots,
                          mu = means,
                          Sigma = sigma_matrix) |>
        as.data.frame()

      data.table::setDT(out[[i]])

      if(!is.na(cloud_filtering_rule)){
        out[[i]] <- out[[i]][eval(parse(text = cloud_filtering_rule)), ]
      }

      n_boots <- n_boots + nrow(out[[i]])
      i <- i + 1
    } # *END WHILE LOOP*
    out <- rbindlist(out)
    out <- out[1:SAFE_boots, ] # Trim to SAFE_boots

    if(family %in% c("4_multivariate_normal")){
      out[, `:=` (sd1 = sqrt(v1), sd2 = sqrt(v2))]
      out[, `:=` (v1 = NULL, v2 = NULL)]
    }

    return(out)
  }else if(family %in% c("4_multivariate_normal_wishart")){

    out <- MASS::mvrnorm(n = SAFE_boots,
                   mu = means,
                   Sigma = (sigma_matrix / c(input$n1, sqrt(input$n1*input$n2), sqrt(input$n1*input$n2), input$n2))) |>
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
  if(family == "2_binomial"){
    if(!all(c("n1", "n2") %in% names(input))){
      input$n1 <- input$a + input$b
      input$n2 <- input$c + input$d
    }
    out <- data.table::data.table(a = stats::rbinom(SAFE_boots, input$n1, input$a / input$n1) |>
                                    as.double(),
                                  c = stats::rbinom(SAFE_boots, input$n2, input$c / input$n2) |>
                                    as.double())
    out[,  `:=` (b = input$n1 - a,
                 d = input$n2 - c)]
    out

    # Add 0.5 to rows with ANY zero
    out[(a == 0 | b == 0 | c == 0 | d == 0),]
    out[(a == 0 | b == 0 | c == 0 | d == 0),
        `:=` (a = a + 0.5,
              b = b + 0.5,
              c = c + 0.5,
              d = d + 0.5)]
    return(out)

  }else if(family == "2_multinomial"){
    # Also need to remove 0s from this...
    N <- (input$n_AA + input$n_Aa + input$n_aa)
    out <- list()

    n_boots <- 0
    i <- 1
    s <- Sys.time()
    while(n_boots < SAFE_boots &
          as.numeric(Sys.time() - s) < 60){

      out[[i]] <- stats::rmultinom(n = SAFE_boots,
                            size = N,
                            prob = c(n_AA = input$n_AA/N,
                                     n_Aa = input$n_Aa/N,
                                     n_aa = input$n_aa/N)) |>
        t() |> # For some reason these are returned WIDE, with 3 rows and 1e6 columns. Weird. Was freezing computer
        as.data.frame()

      data.table::setDT(out[[i]])

      if(!is.na(cloud_filtering_rule)){
        out[[i]] <- out[[i]][eval(parse(text = cloud_filtering_rule)), ]
      }

      n_boots <- n_boots + nrow(out[[i]])
      i <- i + 1
    } # *END WHILE LOOP*
    out <- data.table::rbindlist(out)
    out <- out[1:SAFE_boots, ] # Trim to SAFE_boots

    return(out)
  }

  # Other distributions -----------------------------------------------------
  # "2_multivariate_lognormal",
  # "2_multivariate_Gamma",
  # "2_multivariate_negative_binomial",
  # "2_multivariate_Poisson",
  # "2_multivariate_Beta"
  if(family == "2_multivariate_lognormal"){


  }else if(family == "2_multivariate_Gamma"){


  }else if(family == "2_multivariate_negative_binomial"){

  }else if(family == "2_multivariate_Poisson"){

  }else if(family == "2_multivariate_Beta"){

  }




}
