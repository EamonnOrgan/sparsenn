#' Generate all weight symmetries
#'
#' @param W weight vector
#' @param p number of input nodes
#' @param q number of hidden nodes
#'
#' @return matrix containing all weight symmetries for W
#' @export
generate_symmetries <- function(W, p, q) {

  # number of symmetries
  n_symmetries <- 2 ^ q * factorial(q)

  W_symmetries <- matrix(NA, nrow = n_symmetries, ncol = length(W))

  # list of permunations of hidden nodes
  perm_list <- combinat::permn(1:q)

  # row index for sign flipping after performing node flipping
  row_ind <- sapply(1:(2 ^ q), \(x) 1 + factorial(q) * (x - 1))

  # generate list for which nodes need signflipping
  sign_flip <- c(0, unlist(lapply(1:q, combinat::combn, x = 1:q, simplify = FALSE),
                           recursive = FALSE))

  W_symmetries[1, ] <- W

  W_mat <- cbind(
    matrix(W_symmetries[1, 1:(q * (p + 1))], byrow = TRUE, ncol = (p + 1)),
    W_symmetries[1, ((p + 1) * q + 2):((p + 1) * q + q + 1)])

  for (i in 2:length(perm_list)) {

    W_mat_temp <- W_mat[perm_list[[i]], ]

    W_symmetries[i, ] <- c(t(W_mat_temp[, 1:(p + 1)]),
                           W_symmetries[1, (p + 1) * q + 1],
                           t(W_mat_temp[, (p + 2)]))
  }


  for (j in 2:(2 ^ q)) {


    # gets indicies of W that need to be multiplied by a minus
    ind_vec <- sapply(sign_flip[j],
                      \(y) sapply(1:length(y),
                                  \(x) c(1:(p + 1) + (y[x] - 1) * (p + 1),
                                         (p + 1) * q + y[x]
                                         + 1)))

    W_symmetries[row_ind[j], ] <- W_symmetries[1, ]

    W_symmetries[row_ind[j], ind_vec] <- - W_symmetries[row_ind[j], ind_vec]

    W_symmetries[row_ind[j], (p + 1) * q + 1] <-
      W_symmetries[row_ind[j], (p + 1) * q + 1] +
      sum(W_symmetries[1, (p + 1) * q + sign_flip[[j]] + 1])

    W_mat <- cbind(matrix(
      W_symmetries[row_ind[j], 1:(q * (p + 1))], byrow = TRUE, ncol = (p + 1)),
      W_symmetries[row_ind[j], ((p + 1) * q + 2):((p + 1) * q + q + 1)])

    for (i in 2:length(perm_list)) {

      W_mat_temp <- W_mat[perm_list[[i]], ]

      W_symmetries[row_ind[j] + i - 1, ] <- c(t(W_mat_temp[, 1:(p + 1)]),
                                              W_symmetries[row_ind[j], (p + 1) * q + 1],
                                              t(W_mat_temp[, (p + 2)]))
    }
  }

  return(W_symmetries)
}


#' Compare two weight vectors to match up symmetries
#'
#' @param W_pred weight vector to change
#' @param W_true comparison weight vector
#' @param p number of input nodes
#' @param q number of hidden nodes
#' @param input_only only compare using input-to-hidden weights (useful if
#' hidden-to-output weights grow too large due to redundancy)
#'
#' @return matrix containing all weight symmetries for W
#' @export
compare_symmetries <- function(W_pred, W_true, p, q, input_only = FALSE) {

  symmetries <- generate_symmetries(W_pred, p, q)

  # ensure all signs in output layer the same
  filter_ind <- which(
    apply(
      symmetries[, ((p + 1) * q + 2):((p + 2) * q + 1)],
      1,
      \(x) all.equal(sign(x),
                     ifelse(sign(W_true[((p + 1) * q + 2):((p + 2) * q + 1)]) == 0,
                            1, sign(W_true[((p + 1) * q + 2):((p + 2) * q + 1)])))
    ) == "TRUE")

  if (input_only) {
    symmetries <- symmetries[, 1:((p + 1) * q)]
    W_true <- W_true[1:((p + 1) * q)]
  } else {
    symmetries <- symmetries[filter_ind,]
  }

  edist <- rowSums(
    (matrix(rep(W_true, nrow(symmetries)), nrow = nrow(symmetries), byrow = TRUE)
     - symmetries) ^ 2)

  ind <- which.min(edist)

  return(symmetries[ind, ])

}

#' Generate covariance matrix for MVN data
#'
#' @param p number of input nodes
#' @param rho correlation of AR(1) process
#'
#' @return covariance matrix
#' @export
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) -
                    (1:p - 1))
  rho^exponent
}
