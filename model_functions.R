#' The Kendall's \eqn{\tau} formulas are list below:
#'
#' The Clayton copula Kendall's \eqn{\tau = \eta/(2+\eta)}.
#'
#' The Gumbel copula Kendall's \eqn{\tau = 1 - 1/\eta}.
#'
#' The Frank copula Kendall's \eqn{\tau = 1+4\{D_1(\eta)-1\}/\eta},
#' in which \eqn{D_1(\eta) = \frac{1}{\eta} \int_{0}^{\eta} \frac{t}{e^t-1}dt}.
#'
#' The AMH copula Kendall's \eqn{\tau =  1-2\{(1-\eta)^2 \log (1-\eta) + \eta\}/(3\eta^2)}.
#'
#' The Joe copula Kendall's \eqn{\tau = 1 - 4 \sum_{k=1}^{\infty} \frac{1}{k(\eta k+2)\{\eta(k-1)+2\}}}.
#'
#' The Two-parameter copula (\code{Copula2}) Kendall's \eqn{\tau = 1-2\alpha\kappa/(2\kappa+1)}. \cr
#'
tau_copula <- function(eta, copula){
  
  if (tolower(copula) == "copula2"){
    alpha <- eta[1]
    kappa <- eta[2]
    output <- 1 - 2 * alpha * kappa/(1 + 2 * kappa)
  }
  
  if (tolower(copula) != "copula2") {
    output <- tau(archmCopula(tolower(copula), param = eta, dim = 2))
  }
  
  return(output)
}
