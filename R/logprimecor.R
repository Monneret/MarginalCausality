#' Compute the gradient of the log-likelihood for the correlation-upstream model.
#'
#' @export
#' @param wt Observationnal data, matrix of dimension N*2, where N stand for the number of replicate.
#' The first column correspond to the knocked-out gene.
#' @param ko Interventionnal data, matrix of dimension N*2, where N stand for the number of replicate.
#' The first column correspond to the knocked-out gene.
#' @param par Set of parameter, such that the first correspond to alpha,
#' the second and third to the mean, and the last two ones to the log of the standard deviation
#' @return The value of the gradient of the log likelihood for the correlation-upstream model.
#' @examples
#' wt=rnorm(n=10,mean=-2,sd=0.8)
#' wt=cbind(wt,2*wt+rnorm(n=10,mean=1,sd=1.7))
#' ko=rnorm(n=10,mean=0,sd=0.01)
#' ko=cbind(wt,2*wt+rnorm(n=10,mean=1,sd=1.7))
#' par=c(-2,1,0.3,-4,10)
#' logprime.cor(wt,ko,par)

logprime.cor <- function(wt,ko,par) {
  alpha <- par[1]; mu <- par[2:3]; sigma <- exp(par[4:5])
  pmuX <- sum(wt[,2]-mu[2]-alpha*wt[,1])/sigma[2]^2+
    sum(ko[,2]-mu[2]-alpha*mu[1])/(alpha^2*sigma[1]^2+sigma[2]^2)
  pmuG <- sum(wt[,1]-mu[1])/sigma[1]^2+
    alpha*sum(ko[,2]-mu[2]-alpha*mu[1])/(alpha^2*sigma[1]^2+sigma[2]^2)
  psigmaX <- sum(-1/(sigma[2])+(wt[,2]-mu[2]-alpha*wt[,1])^2/(sigma[2]^3))+
    sum(-sigma[2]/(alpha^2*sigma[1]^2+sigma[2]^2)+(sigma[2]*(ko[,2]-mu[2]-alpha*mu[1])^2)/(alpha^2*sigma[1]^2+sigma[2]^2)^2)
  psigmaG <- sum(-1/sigma[1]+(wt[,1]-mu[1])^2/sigma[1]^3)+
    sum(-alpha^2*sigma[1]/(alpha^2*sigma[1]^2+sigma[2]^2)+(alpha^2*sigma[1]*(ko[,2]-mu[2]-alpha*mu[1])^2)/(alpha^2*sigma[1]^2+sigma[2]^2)^2)
  palpha <- sum(wt[,1]*(wt[,2]-mu[2]-alpha*wt[,1]))/sigma[2]^2+
    sum(-alpha*sigma[1]^2/(alpha^2*sigma[1]^2+sigma[2]^2)+
          (ko[,2]-mu[2]-alpha*mu[1])*(alpha*sigma[1]^2*(ko[,2]-mu[2])+mu[1]*sigma[2]^2)/(alpha^2*sigma[1]^2+sigma[2]^2)^2)
  return(c(pmuX = pmuX,pmuG = pmuG,psigmaX = psigmaX,psigmaG = psigmaG,palpha = palpha))
}
