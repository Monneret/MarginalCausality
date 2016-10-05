require(mvtnorm)
set.seed(42)
n=500
##down
sim.down=function(n,alpha=0.3,mu=c(2.0,-3.1),sigma=c(0.3,0.5)){
  wt=ko=matrix(NA,n,2)
  wt[,1]=rnorm(n,mu[1],sigma[1])
  ko[,1]=rnorm(n,0,10e-10)
  wt[,2]=rnorm(n,alpha*wt[,1]+mu[2],sigma[2])
  ko[,2]=rnorm(n,alpha*ko[,1]+mu[2],sigma[2])
  return(list(n=n,theta=list(alpha=alpha,mu=mu,lsigma=log(sigma)),wt=wt,ko=ko))
}

nrep=100
est=matrix(NA,nrep,5)
for (rep in 1:nrep) {
  sim=sim.down(n=n,alpha=-0.6)
  wt=sim$wt; ko=sim$ko
  ld <- function(par){
    return(loglik.down(wt,ko,par))
  }
  opt1=optim(fn=ld,par=c(0,apply(wt,2,mean),log(apply(wt,2,sd))))
  est[rep,]=opt1$par
}
res=rbind(unlist(sim$theta),apply(est,2,mean),apply(est,2,sd))
rownames(res)=c("ref","mean","sd")

test_that("Test the value of down parameters asymptotically", {
  expect_true(abs(res[2,1]-res[1,1])<0.01)
  expect_true(abs(res[2,2]-res[1,2])<0.01)
  expect_true(abs(res[2,3]-res[1,3])<0.01)
  expect_true(abs(res[2,4]-res[1,4])<0.01)
  expect_true(abs(res[2,5]-res[1,5])<0.01)
})

#cor
sim.cor=function(n,alpha=0.3,mu=c(2.0,-3.1),sigma=c(0.3,0.5)){
  m=c(mu[1],alpha*mu[1]+mu[2])
  s=c(sigma[1],sqrt(alpha^2*sigma[1]^2+sigma[2]^2))
  rho=alpha*s[1]/s[2]
  S=tcrossprod(s)*matrix(c(1,rho,rho,1),2,2)
  wt=rmvnorm(n,m,S)
  ko=rmvnorm(n,m,S)
  ko[,1]=rnorm(n,0,10e-10)
  return(list(n=n,theta=list(alpha=alpha,mu=mu,lsigma=log(sigma)),wt=wt,ko=ko))
}

nrep=100
est=matrix(NA,nrep,5)
for (rep in 1:nrep) {
  sim=sim.cor(n=n,alpha=-0.6)
  wt=sim$wt; ko=sim$ko
  lc <- function(par){
    return(loglik.cor(wt,ko,par))
  }
  opt0=optim(fn=lc,par=c(0,apply(wt,2,mean),log(apply(wt,2,sd))))
  est[rep,]=opt0$par
}
res=rbind(unlist(sim$theta),
          apply(est,2,mean),
          apply(est,2,sd))
rownames(res)=c("ref","mean","sd")

test_that("Test the value of up-cor parameters asymptotically", {
  expect_true(abs(res[2,1]-res[1,1])<0.01)
  expect_true(abs(res[2,2]-res[1,2])<0.01)
  expect_true(abs(res[2,3]-res[1,3])<0.01)
  expect_true(abs(res[2,4]-res[1,4])<0.01)
  expect_true(abs(res[2,5]-res[1,5])<0.01)
})

detach(name="package:mvtnorm",unload=TRUE)
