#' Compute bayes factor associated with a marginal linear causal model.
#'
#' @export
#' @param wt Observationnal data, matrix of dimension N*2, where N stand for the number of replicate.
#' The first column correspond to the knocked-out gene.
#' @param ko Interventionnal data, matrix of dimension N*2, where N stand for the number of replicate.
#' The first column correspond to the knocked-out gene.
#' @param par Set of initial parameter for the optimisation, such that the first correspond to alpha,
#' the second and third to the mean, and the last two ones to the log of the standard deviation
#' @param gradient Do we use the gradient? Default set to TRUE.
#' @param prior Do we want to promote one model? Default to 0.5.
#' @return The return of the two numerical optimisation, with the optimal value for each parameters, and the value of the bayes factor.
#' @examples
#' wt <- rnorm(n=10,mean=-2,sd=0.8)
#' wt <- cbind(wt,2*wt+rnorm(n=10,mean=1,sd=1.7))
#' ko <- rnorm(n=10,mean=0,sd=0.01)
#' ko <- cbind(wt,2*wt+rnorm(n=10,mean=1,sd=1.7))
#' par <- c(-2,1,0.3,-4,10)
#' post.causal(wt,ko,par)$B

post.causal <- function(wt,ko,par=c(0,apply(wt,2,mean),log(apply(wt,2,sd))),gradient=FALSE,prior=0.5) {
  if (gradient == TRUE){
    ld <- function(par){
      return(loglik.down(wt,ko,par))
    }
    lc <- function(par){
      return(loglik.cor(wt,ko,par))
    }
    gd <- function(par){
      return(logprime.down(wt,ko,par))
    }
    gc <- function(par){
      return(logprime.cor(wt,ko,par))
    }
    opt1 <- optim(fn = ld, gr  =  gd, par = par)
    opt0 <- optim(fn = lc, gr =  gc, par = par)
  } else {
    ld <- function(par){
      return(loglik.down(wt,ko,par))
    }
    lc <- function(par){
      return(loglik.cor(wt,ko,par))
    }
    opt1 <- optim(fn = ld,par = par)
    opt0 <- optim(fn = lc,par = par)
  }
  B <- exp(-opt0$value+opt1$value)
  return(list(cor=opt0,down=opt1,B=B))
}
