#' SIC log-likelihood
#'
#' @param theta_s A vector of the neural network parameters and log(sigma)
#' @param y Response vector
#' @param X Input matrix
#' @param q Number of hidden nodes
#' @param eps Epsilon value
#'
#' @return neagtive penalised log-likelihood value
#' @export
nn_logl_sic <- function(theta_s, y, X, q, eps){

  n <- nrow(X)
  p <- dim(X)[2] - 1

  g <- theta_s[(((p + 1) * q + 1):((p + 2) * q + 1))]
  wt <- matrix(theta_s[1:((p + 1) * q)], nrow = p + 1)

  H <- cbind(1, 1 / (1 + exp(-X%*%wt)))

  sigma2 <- exp(theta_s[length(theta_s)]) ^ 2

  sse <- sum((y - H%*%g)^2)

  logl <- - (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - sse/(2 * sigma2)

  intercepts <- c((p + 1) * (1:q - 1) + 1, (p + 1) * q + 1)
  weights <- c(1:((p + 2) * q + 1))[-intercepts]

  pen <- (log(n) / 2) * (sum(phi(theta_s[weights], eps)) + length(intercepts))

  mlogl_p <- - logl + pen

  # for logl derivative
  d1 <- (- 1 / sigma2) * (y - H %*% g)

  d2 <- d1 %*% t(g[-1]) * sigmoid(X%*%wt) * (1 - sigmoid(X%*%wt))

  # for pen derivative
  weights_ind <- rep(0, ((p + 2) * q + 2))

  weights_ind[weights] <- 1

  theta_tilde <- theta_s * weights_ind

  attr(mlogl_p, "gradient") <-
    c(as.vector(t(X) %*% d2), t(d1) %*% H, n - (sse / (sigma2))) +
    (log(n) / 2) * ((2 * theta_tilde * eps ^ 2) / (theta_tilde ^ 2 + eps ^ 2) ^ 2)



  return(mlogl_p)
}

#' Group SIC log-likelihood for hidden nodes
#'
#' @param theta_s A vector of the neural network parameters and log(sigma)
#' @param y Response vector
#' @param X Input matrix
#' @param q Number of hidden nodes
#' @param eps Epsilon value (default = 0)
#'
#' @return neagtive penalised log-likelihood value
#' @export
nn_logl_gsic_hidden <- function(theta_s, y, X, q, eps = 0){

  n <- nrow(X)
  p <- dim(X)[2] - 1
  g <- theta_s[(((p + 1) * q + 1):((p + 2) * q + 1))]
  wt <- matrix(theta_s[1:((p + 1) * q)], nrow = p + 1)
  H <- cbind(1, 1 / (1 + exp(-X%*%wt)))

  sigma2 <- exp(theta_s[length(theta_s)]) ^ 2

  sse <- sum((y - H%*%g)^2)

  logl <- - (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - sse/(2 * sigma2)

  intercepts <- c((p + 1) * (1:q - 1) + 1, (p + 1) * q + 1)
  weights <- c(1:((p + 2) * q + 1))[-intercepts]

  weights_grouped <- sapply(1:p, \(i)
                            sapply(X = 1:q, FUN = function(x) (x - 1) *
                                     (p + 1) + 1 + i))

  weights_grouped <- t(weights_grouped)

  g_ind <- ((p + 1) * q + 2):((p + 1) * q + 1 + q)

  weights_grouped <- rbind(weights_grouped, g_ind)
  rownames(weights_grouped) <- NULL

  group_pen <- sum(apply(sapply(1:q, \(x) theta_s[weights_grouped[, x]]), 2,
                         phi_g, eps = eps))

  pen <- (log(n) / 2) * (group_pen + length(intercepts))

  mlogl_p <- - logl + pen

  # for logl derivative
  d1 <- (- 1 / sigma2) * (y - H %*% g)

  d2 <- d1 %*% t(g[-1]) * sigmoid(X%*%wt) * (1 - sigmoid(X%*%wt))

  # for pen derivative
  weights_ind <- rep(0, ((p + 2) * q + 2))

  weights_ind[weights] <- 1

  theta_tilde <- theta_s * weights_ind

  theta_tilde_denom <- theta_tilde

  theta_tilde_denom[
    weights[!(weights %in% g_ind)]] <- rep(apply(sapply(1:q,
                                                        \(x) theta_s[weights_grouped[, x]]),
                                                 2, \(y) sum(y ^ 2)), each = p)

  theta_tilde_denom[g_ind] <- apply(sapply(1:q,
                                           \(x) theta_s[weights_grouped[, x]]),
                                    2, \(y) sum(y ^ 2))

  p_vec <- c(rep(c(0, rep(p + 1, times = p)), times = q), 0, rep(p + 1, times = q), 0)


  attr(mlogl_p, "gradient") <-
    c(as.vector(t(X) %*% d2), t(d1) %*% H, n - (sse / (sigma2))) +
    (log(n) / 2) * ((2 * theta_tilde * eps ^ 2) / (theta_tilde_denom ^ 2 + eps ^ 2) ^ 2) * p_vec

  return(mlogl_p)
}

#' Group SIC log-likelihood for input nodes
#'
#' @param theta_s A vector of the neural network parameters and log(sigma)
#' @param y Response vector
#' @param X Input matrix
#' @param q Number of hidden nodes
#' @param eps Epsilon value (default = 0)
#'
#' @return neagtive penalised log-likelihood value
#' @export
nn_logl_gsic_input <- function(theta_s, y, X, q, eps = 0){

  n <- nrow(X)
  p <- dim(X)[2] - 1
  g <- theta_s[(((p + 1) * q + 1):((p + 2) * q + 1))]
  wt <- matrix(theta_s[1:((p + 1) * q)], nrow = p + 1)
  H <- cbind(1, 1 / (1 + exp(-X%*%wt)))

  sigma2 <- exp(theta_s[length(theta_s)]) ^ 2

  sse <- sum((y - H%*%g)^2)

  logl <- - (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - sse/(2 * sigma2)

  intercepts <- c((p + 1) * (1:q - 1) + 1, (p + 1) * q + 1)
  weights <- c(1:((p + 2) * q + 1))[-intercepts]

  weights_grouped <- sapply(1:p, \(i)
                            sapply(X = 1:q, FUN = function(x) (x - 1) *
                                     (p + 1) + 1 + i))

  g_ind <- ((p + 1) * q + 2):((p + 1) * q + 1 + q)

  group_pen <- sum(apply(sapply(1:p, \(x) theta_s[weights_grouped[, x]]), 2,
                         phi_g, eps = eps))
  g_pen <- sum(phi(theta_s[g_ind], eps))

  pen <- (log(n) / 2) * (group_pen + g_pen + length(intercepts))

  mlogl_p <- - logl + pen

  # for logl derivative
  d1 <- (- 1 / sigma2) * (y - H %*% g)

  d2 <- d1 %*% t(g[-1]) * sigmoid(X%*%wt) * (1 - sigmoid(X%*%wt))

  # for pen derivative
  weights_ind <- rep(0, ((p + 2) * q + 2))

  weights_ind[weights] <- 1

  theta_tilde <- theta_s * weights_ind

  theta_tilde_denom <- c(rep(c(0,
                               apply(sapply(1:p,
                                            \(x) theta_s[weights_grouped[, x]]),
                                     2, \(y) sum(y ^ 2))),
                             times = q),
                         theta_tilde[c(((p + 1) * q + 1):(((p + 2) * q + 1)))],
                         0)
  q_vec <- c(rep(c(0, rep(q, times = p)), times = q), 0, rep(1, times = q), 0)


  attr(mlogl_p, "gradient") <-
    c(as.vector(t(X) %*% d2), t(d1) %*% H, n - (sse / (sigma2))) +
    (log(n) / 2) * ((2 * theta_tilde * eps ^ 2) / (theta_tilde_denom ^ 2 + eps ^ 2) ^ 2) * q_vec

  return(mlogl_p)
}
