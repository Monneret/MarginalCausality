set.seed(27)
n=50

wt <- rnorm(n=n,mean=-2,sd=0.02)
wt <- cbind(wt,2*wt+rnorm(n=n,mean=1,sd=0.07))
ko <- rnorm(n=n,mean=0,sd=10e-5)
ko <- cbind(ko,2*ko+rnorm(n=n,mean=1,sd=0.07))
par <- c(-2,1,0.3,-4,10)

Bdown=log10(post.causal(wt,ko,par=c(0,apply(wt,2,mean),log(apply(wt,2,sd))),gradient=FALSE,prior=0.5)$B)

wt <- rnorm(n=n,mean=-2,sd=0.02)
wt <- cbind(2*wt+rnorm(n=n,mean=1,sd=0.07),wt)
ko <- cbind(rnorm(n=n,mean=0,sd=10e-5),rnorm(n=n,mean=-2,sd=0.02))
Bup=log10(post.causal(wt,ko,par=c(0,apply(wt,2,mean),log(apply(wt,2,sd))),gradient=FALSE,prior=0.5)$B)

test_that("Value of a particular bayes factor, for a clearly down model and a clearly cor model", {
  expect_true(Bdown< -2)
  expect_true(Bup>2)
})
