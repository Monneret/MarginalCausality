% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/undownstream.R
\name{undownstream}
\alias{undownstream}
\title{Compute parameters for the correlation or the upstream model based on downstream parameters.}
\usage{
undownstream(par, to)
}
\arguments{
\item{par}{Parameters of the downstream model}

\item{to}{A character to decide the model of the returning parameters. Accept "cor" for correlation model, and "up" for upstream model.}
}
\value{
Parameters of the "to" model.
}
\description{
Compute parameters for the correlation or the upstream model based on downstream parameters.
}
\examples{
wtd <- rnorm(n=10,mean=-2,sd=0.8)
wtd <- cbind(wtd,2*wtd+rnorm(n=10,mean=1,sd=1.7))
par <- c(2,-2,1,log(0.8),log(1.7))
parup <- undownstream(par,to="up")
wtu <- rnorm(n=10,mean=parup[3],sd=exp(parup[5]))
wtu <- cbind(parup[1]*wtu+rnorm(n=10,mean=parup[2],sd=exp(parup[4])),wtu)
parcor <- undownstream(par,to="cor")
cor(wtd[,1],wtu[,1])
cor(wtd[,2],wtu[,2])
parcor[1]
}

