#' Compute parameters for the correlation or the upstream model based on downstream parameters.
#'
#' @export
#' @param par Parameters of the downstream model
#' @param to A character to decide the model of the returning parameters. Accept "cor" for correlation model, and "up" for upstream model.
#' @return Parameters of the "to" model.
#' @examples
#' wtd <- rnorm(n=10,mean=-2,sd=0.8)
#' wtd <- cbind(wtd,2*wtd+rnorm(n=10,mean=1,sd=1.7))
#' par <- c(2,-2,1,log(0.8),log(1.7))
#' parup <- undownstream(par,to="up")
#' wtu <- rnorm(n=10,mean=parup[3],sd=exp(parup[5]))
#' wtu <- cbind(parup[1]*wtu+rnorm(n=10,mean=parup[2],sd=exp(parup[4])),wtu)
#' parcor <- undownstream(par,to="cor")
#' cor(wtd[,1],wtu[,1])
#' cor(wtd[,2],wtu[,2])
#' parcor[1]

undownstream<-function(par,to){
  alpha <- par[1]; mu <- par[2:3]; sigma <- exp(par[4:5])
  if (to=="cor"){
    rho <- alpha*sigma[1]/sqrt(alpha^2*sigma[1]^2+sigma[2]^2)
    m <- c(mu[1],mu[2]+alpha*mu[1])
    s <- c(sigma[1],sqrt(alpha^2*sigma[1]^2+sigma[2]^2))
    return(par=c(rho=rho,m=m,s=s))
  } else if (to=="up"){
    beta <- alpha*sigma[1]^2/(alpha^2*sigma[1]^2+sigma[2]^2)
    mutilde <- c(mu[1]-alpha*sigma[1]^2*(mu[2]+alpha*mu[1])/(alpha^2*sigma[1]^2+sigma[2]^2),mu[2]+alpha*mu[1])
    sigmatilde <- c(sqrt(sigma[1]^2-alpha^2*sigma[1]^4/(alpha^2*sigma[1]^2+sigma[2]^2)),sqrt(alpha^2*sigma[1]^2+sigma[2]^2))
    return(par=c(beta=beta,mu=mutilde,sigma=sigmatilde))
  } else {
    stop("Specify a good model, either \'cor\' or \'up\'.")
  }
}
