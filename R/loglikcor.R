#' Compute the likelihood for the correlation-upstream model.
#'
#' @export
#' @param wt Observationnal data, matrix of dimension N*2, where N stand for the number of replicate.
#' The first column correspond to the knocked-out gene.
#' @param ko Interventionnal data, matrix of dimension N*2, where N stand for the number of replicate.
#' The first column correspond to the knocked-out gene.
#' @param par Set of parameter, such that the first correspond to alpha,
#' the second and third to the mean, and the last two ones to the log of the standard deviation
#' @return The value of the log likelihood for the correlation-upstream model.
#' @examples
#' wt <- rnorm(n=10,mean=-2,sd=0.8)
#' wt <- cbind(wt,2*wt+rnorm(n=10,mean=1,sd=1.7))
#' ko <- rnorm(n=10,mean=0,sd=0.01)
#' ko <- cbind(wt,2*wt+rnorm(n=10,mean=1,sd=1.7))
#' par <- c(-2,1,0.3,-4,10)
#' loglik.cor(wt,ko,par)

loglik.cor <- function(wt,ko,par) {
  alpha <- par[1]; mu <- par[2:3]; sigma <- exp(par[4:5])
  res <- sum(dnorm(wt[,1],mu[1],sigma[1],log <- TRUE))
  res <- res+sum(dnorm(wt[,2],alpha*wt[,1]+mu[2],sigma[2],log <- TRUE))
  res <- res+sum(dnorm(ko[,2],alpha*mu[1]+mu[2],
                    sqrt(alpha^2*sigma[1]^2+sigma[2]^2),log <- TRUE))
  return(-res)
}
