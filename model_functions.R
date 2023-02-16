# Tranformation function; 
# para is . inside G(.)
# r is the parameter that determines function, as illustrated below.
G_fun <- function(para,r) {
  
  if (r <= 2) { # Box-Cox transformations; if r = 1, then PH model
    output = (1/r)*((1+para)^r-1)
  }
  
  else { # Logarithmic transfomations; if r = 3, then PO model
    output = (1/(r-2))*log(1+para*(r-2)) 
  }
  return(output)
}
# Bernstein polynomial function: j from 0 to m, where m is degree; 
# [l,u] is specified range of time;
# t for observed times
# output is a Bernstein polynomial as function of t
bern <- function(j,m,l,u,t){
  b = (factorial(m)/(factorial(j)*factorial(m-j)))*(((t-l)/(u-l))^j)*((1-(t-l)/(u-l))^(m-j))
  return(b)
}

# n<-10
# k<-4
# choose(n,k)
# factorial(n)/(factorial(k)*factorial(n-k))

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