#' Neural network SIC
#'
#' @param X Input data
#' @param y Response vector
#' @param q Number of hidden nodes
#' @param method One of c("single", "group", "twostep")
#' @param grouping Grouping for "group" method
#' @param eps_start Starting epsilon value for telescoping
#' @param t Number of epsilon values
#' @param r Rate of decay for epsilon
#' @param n_iter Number of iterations for nlm after first step
#' @param n_init Number of random initialisations for first step
#' @param tau Value to swap from group to single method in "twostep"
#'
#' @return Function value
#' @export
nn_sic <- function(X, y, q, method = c("single", "group", "twostep"),
                   grouping = c("input", "hidden"),
                   eps_start = 10, t = 100, r = 0.87, n_iter = 1,
                   n_init = 1, tau = 0) {

  # check if X is design matrix
  if(all(X[, 1] == 1)) {
    p <- ncol(X) - 1
  } else {
    p <- ncol(X)
    X <- cbind(1, X)
  }

  if (method[1] == "single") {
    nn_logl <- nn_logl_sic
  } else if (method[1] == "group" & grouping == "input") {
    nn_logl <- nn_logl_gsic_input
  } else if (method[1] == "group" & grouping == "hidden") {
    nn_logl <- nn_logl_gsic_hidden
  } else if (method[1] == "twostep" & grouping == "input") {
    nn_logl <- nn_logl_gsic_input
    nn_logl2 <- nn_logl_sic
  } else if (method[1] == "twostep" & grouping == "hidden") {
    nn_logl <- nn_logl_gsic_hidden
    nn_logl2 <- nn_logl_sic
  }

  if (tau != 0 & method[1] != "twostep") {
    warning("Warning: Change point not applicable when method is not twostep.")
    tau <- 0
  }

  n <- nrow(X)

  eps <- sapply(1:t, \(x) eps_start * r ^ (x - 1))

  theta_start <- matrix(rnorm(n_init * ((p + 2) * q + 2), sd = 1/3),
                        nrow = n_init, byrow = TRUE)

  weight_matrix <- matrix(NA, nrow = t, ncol = (p + 2) * q + 2)

  for (i in 1:(length(eps) - tau)) {

    if (i == 1) {

      min_temp <- Inf

      for (j in 1:n_init) {
        nlm_nn_sic <- nlm(nn_logl, theta_start[j, ], y, X, q, eps = eps[i],
                          iterlim = 2000)

        if (nlm_nn_sic$minimum < min_temp) {
          weight_matrix[i, ] <- nlm_nn_sic$estimate
          min_temp <- nlm_nn_sic$minimum
        }

      }

    } else {

      nlm_nn_sic <- nlm(nn_logl, weight_matrix[i - 1, ], y, X, q, eps = eps[i],
                        iterlim = n_iter, check.analyticals = FALSE)
      weight_matrix[i, ] <- nlm_nn_sic$estimate

    }
  }

  if(method[1] == "twostep" & tau != 0) {
    for (i in (tau + 1):length(eps)) {

      nlm_nn_sic <- nlm(nn_logl2, weight_matrix[i - 1, ], y, X, q, eps = eps[i],
                        iterlim = n_iter)

      weight_matrix[i, ] <- nlm_nn_sic$estimate

    }
  }



  return(list("estimate" = nlm_nn_sic$estimate,
              "weight_matrix" = weight_matrix))

}
