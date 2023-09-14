#' Smooth approximation to L0 norm
#'
#' @param x Value
#' @param eps Epsilon value (default = 0)
#'
#' @return Function value
#' @export
phi <- function (x, eps = 0) {
  return(x ^ 2 / (x ^ 2 + eps ^ 2))
}

#' Smooth approximation to L0 norm for groups
#'
#' @param x Value
#' @param eps Epsilon value (default = 0)
#'
#' @return Function value
#' @export
phi_g <- function (x, eps = 0) {
  (sum(x ^ 2) / (sum(x ^ 2) + eps ^ 2)) * length(x)
}
