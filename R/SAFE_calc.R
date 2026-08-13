#' SAFE calculation
#'
#' Used internally to calculate SAFE effect sizes. This function calls parameter_cloud and then transforms hyperparameters
#'
#' @import data.table
#' @importFrom stats sd
#' @param formulas blah
#' @param input blah
#' @param plugin_effect blah
#' @param custom_sigma blah
#' @param SAFE_boots blah
#' @param SAFE_max_secs blah
#' @return A data.table with effect sizes and sample variances calculated by SAFE.
.SAFE_calc <- function(formulas,
                       input,
                       plugin_effect,
                       custom_sigma,
                       SAFE_boots = 1e6,
                       SAFE_max_secs){

  # If paired, add n1 and n2 to input...
  if(formulas$paired_design == TRUE){
    input$n1 <- input$n
    input$n2 <- input$n
  }

  # Construct sigma matrices ------------------------------------------------
  if(formulas$SAFE_family %in% "1_normal"){
    if(is.null(custom_sigma)){
      sigma_matrix <- input$sd / sqrt(input$n)
    }
    means <- c(x = input$x)

  }else if(formulas$SAFE_family %in% "2_multivariate_normal"){
    if(is.null(custom_sigma)){

      sigma_matrix <- matrix(data = c((input$sd1^2 / input$n1),                    (input$r*input$sd1*input$sd2)/input$n1, #  / n1 add this to sd1^2
                                      (input$r*input$sd1*input$sd2)/input$n1,      (input$sd2^2 / input$n2)), #  / n2 add this to sd2^2
                             nrow = 2, ncol = 2)

    }
    means <- c(x1 = input$x1, x2 = input$x2)

  }else if(formulas$SAFE_family == "4_multivariate_normal_wishart"){
    if(is.null(custom_sigma)){

      sigma_matrix <- matrix(c(input$sd1^2, input$r*input$sd1*input$sd2,
                               input$r*input$sd1*input$sd2, input$sd2^2),
                             2, 2)

    }
    means <- c(x1 = input$x1, x2 = input$x2)

    # means <- c(x1 = input$x1, x2 = input$x2)
  }else if(formulas$SAFE_family == "4_multivariate_normal"){
    if(is.null(custom_sigma)){

      sigma_matrix <- matrix(data = c(input$sd1^2/input$n1,                  (input$r*input$sd1*input$sd2)/input$n1, 0,                                                  0,
                                      (input$r*input$sd1*input$sd2)/input$n1, input$sd2^2/input$n2,                   0,                                                  0,
                                      0,                                      0,                                      (2*input$sd1^4)/(input$n1-1),                       ((2*input$r^2*input$sd1^2*input$sd2^2)/(input$n1-1)),
                                      0,                                      0,                                      (2*input$r^2*input$sd1^2*input$sd2^2)/(input$n1-1), (2*input$sd2^4)/(input$n2-1)),
                             nrow = 4,
                             ncol = 4)

    }
    means <- c(x1 = input$x1, x2 = input$x2,
               v1 = input$sd1^2, v2 = input$sd2^2)

  }else if(formulas$SAFE_family %in% c("2_multinomial_as_normal")){

    if(is.null(custom_sigma)){

      if(!"n1" %in% names(input)){
        input$n1 <- input$a + input$b
        input$n2 <- input$c + input$d
      }
      input$p1 <- input$a / input$n1
      input$p2 <- input$c / input$n2

      # This is variance, which is what mvrnorm wants:
      input$v1 <- input$p1 * (1 - input$p1) #/ input$n1
      input$v2 <- input$p2 * (1 - input$p2) #/ input$n2
      input$r <- 0

      # The top-left and bottom-right corneres formerly were sd1^2 / n1 & sd2^2 / n2
      sigma_matrix <- matrix(data = c((input$v1 / input$n1),                    (input$r*input$v1*input$v2)/input$n1, #  / n1 add this to sd1^2
                                      (input$r*input$v1*input$v2)/input$n1,      (input$v2 / input$n2)), #  / n2 add this to sd2^2
                             nrow = 2, ncol = 2)

    }
    means <- c(p1 = input$a/input$n1, p2 = input$c/input$n2)
  }

  if(!is.null(custom_sigma)){
    sigma_matrix <- custom_sigma
  }

  # Cloud generating while loop -------------------------------------------------------------

  cloud_list <- list()
  cloud_length <- 0
  boots_remaining <- SAFE_boots
  i <- 1
  start_time <- Sys.time()

  while(cloud_length < SAFE_boots &&
        as.numeric(Sys.time() - start_time) < SAFE_max_secs){

    # Create Gaussian cloud_list[[i]]s ------------------------------------------------------------
    if(unique(formulas$SAFE_family == "1_normal")){

      cloud_list[[i]] <- data.table::data.table(x = rnorm(n=boots_remaining,
                                              mean = var_guide$mean,
                                              sd = sigma_matrix))

    }else if(formulas$SAFE_family %in% c("4_multivariate_normal",
                                         "2_multivariate_normal",
                                         "2_multinomial_as_normal")){

      cloud_list[[i]] <- MASS::mvrnorm(n = boots_remaining,
                                       mu = means,
                                       Sigma = sigma_matrix) |>
        as.data.frame() |>
        data.table::setDT()
      names(cloud_list[[i]]) <- names(means)

      # Back convert the variance hyperparameters to SD
      if(formulas$SAFE_family == "4_multivariate_normal"){
        cloud_list[[i]][, `:=` (sd1 = sqrt(v1), sd2 = sqrt(v2))]
        cloud_list[[i]][, `:=` (v1 = NULL, v2 = NULL)]
      }
      if(formulas$SAFE_family == "2_multinomial_as_normal"){
        cloud_list[[i]][, `:=` (n1 = input$n1, n2 = input$n2)]
        cloud_list[[i]][, `:=` (a = round(p1 * n1),
                                c = round(p2 * n2))]
        cloud_list[[i]][, `:=` (b = n1 - a,
                                d = n2 - c)]
      }
      #' [I really don't like this degree of specificity of effect_type manipulation inside the function]
      if(formulas$effect_size == "lnRR"){
        cloud_list[[i]][a == 0, `:=` (a = a + 0.5,
                                      n1 = n1 + 1) ]
        cloud_list[[i]][c == 0, `:=` (c = c + 0.5,
                                      n2 = n2 + 1) ]
      }
      if(formulas$effect_size == "lnOR"){
        cloud_list[[i]][(a == 0 | b == 0 | c == 0 | d == 0), `:=`
                        (a = a + 0.5,
                          b = b + 0.5,
                          c = c + 0.5,
                          d = d + 0.5)]
      }

    }else if(formulas$SAFE_family %in% c("4_multivariate_normal_wishart")){

      cloud_list[[i]] <- MASS::mvrnorm(n = boots_remaining,
                                       mu = means,
                                       Sigma = (sigma_matrix / c(input$n1, sqrt(input$n1*input$n2), sqrt(input$n1*input$n2), input$n2))) |>
        as.data.frame() |>
        data.table::setDT()
      names(cloud_list[[i]]) <- names(means)

      #
      wishart.cloud <-  stats::rWishart(boots_remaining,
                                        df = (input$n1-1),
                                        Sigma = sigma_matrix)

      cloud_list[[i]][, sd1 := sqrt(wishart.cloud[1, 1, ] / (input$n1 - 1))]
      cloud_list[[i]][, sd2 := sqrt(wishart.cloud[2, 2, ] / (input$n2 - 1))]
      cloud_list[[i]]

    }

    # Create Binomial clouds --------------------------------------------------------------
    if(formulas$SAFE_family == "2_binomial"){ # lnRR

      cloud_list[[i]] <- data.table::data.table(a = rbinom(boots_remaining, input$n1, input$a / input$n1) |> as.double(),
                                                c = rbinom(boots_remaining, input$n2, input$c / input$n2) |> as.double())
      cloud_list[[i]][, n1 := input$n1]
      cloud_list[[i]][, n2 := input$n2]

      cloud_list[[i]][a == 0, `:=` (a = a + 0.5,
                                    n1 = n1 + 1) ]
      cloud_list[[i]][c == 0, `:=` (c = c + 0.5,
                                    n2 = n2 + 1) ]

    }else if(formulas$SAFE_family == "4_binomial"){ # this is lnOR
      if(!all(c("n1", "n2") %in% names(input))){
        input$n1 <- input$a + input$b
        input$n2 <- input$c + input$d
      }

      cloud_list[[i]] <- data.table::data.table(a = rbinom(boots_remaining, input$n1, input$a / input$n1) |> as.double(),
                                                c = rbinom(boots_remaining, input$n2, input$c / input$n2) |> as.double())

      cloud_list[[i]][, `:=` (b = input$n1 - a,
                              d = input$n2 - c)]

      cloud_list[[i]][(a == 0 | b == 0 | c == 0 | d == 0), `:=`
                      (a = a + 0.5,
                        b = b + 0.5,
                        c = c + 0.5,
                        d = d + 0.5)]

    }else if(formulas$SAFE_family == "3_multinomial"){
      N <- (input$n_AA + input$n_Aa + input$n_aa)
      cloud_list[[i]] <- stats::rmultinom(n = boots_remaining,
                                          size = N,
                                          prob = c(n_AA = input$n_AA/N,
                                                   n_Aa = input$n_Aa/N,
                                                   n_aa = input$n_aa/N)) |>
        t() |>
        as.data.frame()

      data.table::setDT(cloud)
      cloud_list[[i]][, `:=` (n_AA = as.double(n_AA),
                              n_Aa = as.double(n_Aa),
                              n_aa = as.double(n_aa))]

      cloud_list[[i]][(n_AA == 0 | n_Aa == 0 | n_aa == 0),
                      `:=` (n_AA = n_AA + 0.5,
                            n_Aa = n_Aa + 0.5,
                            n_aa = n_aa + 0.5)]
    }

    # Filter  -----------------------------------------------------------------
    cloud_list[[i]] <- cloud_list[[i]][eval(parse(text = formulas$cloud_filtering_rules)), ]

    # Add missing inputs (e.g., n)
    cloud_list[[i]] <- data.table::data.table(cloud_list[[i]],
                                  input[!names(input) %in% names(cloud_list[[i]])] |> unlist() |> t() |> data.table())
    # Calculate estimates -------------------------------------------------

    # Convert cloud
    cloud_list[[i]]$yi_first <- suppressWarnings(.calc_effect(formulas = formulas[calc_type == "point_estimate" &
                                                                                    derivative == "first", ],
                                                              input = cloud_list[[i]])$yi_first)

    # Filter out NAs:
    cloud_list[[i]] <- cloud_list[[i]][!is.na(yi_first), ]

    cloud_length <- sapply(cloud_list, nrow) |> sum()
    boots_remaining <- SAFE_boots - cloud_length
    i <- i+1

  } # End while loop
  elapsed_time <- Sys.time() - start_time

  cloud <- rbindlist(cloud_list)

  if(nrow(cloud) > 0){

    cloud <- cloud[1:SAFE_boots, ] # Filter extra

    # bias corrected estimate of sampling variance and SE:
    SE_safe <- sd(cloud$yi_first)
    vi_safe <- SE_safe^2

    if(formulas$SAFE_bias_correction == "yes"){
      bias_SAFE <- mean(cloud$yi_first) - plugin_effect

      yi_safe <- plugin_effect - bias_SAFE
    }else{ #lnM doesn't do bias-correction:
      yi_safe <- mean(cloud$yi_first)
    }

    return(data.table::data.table(yi_safe = yi_safe,
                      vi_safe = vi_safe,
                      number_SAFE_iterations = i,
                      number_SAFE_bootstraps = nrow(cloud)),
                      SAFE_elapsed_secs = elapsed_time)
  }else{
    return(data.table::data.table(yi_safe = NA,
                      vi_safe = NA,
                      number_SAFE_iterations = i,
                      number_SAFE_bootstraps = 0,
                      SAFE_elapsed_secs = elapsed_time))
  }
}

