set.seed(42)
n=50
nrep=100

wt <- rnorm(n=n,mean=-2,sd=0.8)
wt <- cbind(wt,2*wt+rnorm(n=n,mean=1,sd=1.7))
ko <- rnorm(n=n,mean=0,sd=0.01)
ko <- cbind(wt,2*wt+rnorm(n=n,mean=1,sd=1.7))
par <- c(-2,1,0.3,-4,10)

dalpha=Galpha=dmuG=GmuG=dmuX=GmuX=dsigmaG=GsigmaG=dsigmaX=GsigmaX=vector(mode="double",length=nrep)

for (i in 1:nrep){
  x=par+c(i/nrep,0,0,0,0)
  dalpha[i]=loglik.down(wt,ko,x)
  Galpha[i]=logprime.down(wt,ko,x)["palpha"]
  x=par+c(0,i/nrep,0,0,0)
  dmuG[i]=loglik.down(wt,ko,x)
  GmuG[i]=logprime.down(wt,ko,x)["pmuG"]
  x=par+c(0,0,i/nrep,0,0)
  dmuX[i]=loglik.down(wt,ko,x)
  GmuX[i]=logprime.down(wt,ko,x)["pmuX"]
  x=par+c(0,0,0,i/nrep,0)
  dsigmaG[i]=loglik.down(wt,ko,x)
  GsigmaG[i]=logprime.down(wt,ko,x)["psigmaG"]
  x=par+c(0,0,0,0,i/nrep)
  dsigmaX[i]=loglik.down(wt,ko,x)
  GsigmaX[i]=logprime.down(wt,ko,x)["psigmaX"]
}

test_that("Test the value of up-cor parameters asymptotically", {
  expect_equal(signif(Galpha[1],digits=2),signif((dalpha[1]-dalpha[2])*nrep,digits=2))
  expect_equal(signif(GmuG[1],digits = 2),signif((dmuG[1]-dmuG[2])*nrep,digits=2))
  expect_equal(signif(GmuX[1],digits=2),signif((dmuX[1]-dmuX[2])*nrep,digits = 2))
})

dalpha=Galpha=dmuG=GmuG=dmuX=GmuX=dsigmaG=GsigmaG=dsigmaX=GsigmaX=vector(mode="double",length=nrep)

for (i in 1:nrep){
  x=par+c(i/nrep,0,0,0,0)
  dalpha[i]=loglik.cor(wt,ko,x)
  Galpha[i]=logprime.cor(wt,ko,x)["palpha"]
  x=par+c(0,i/nrep,0,0,0)
  dmuG[i]=loglik.cor(wt,ko,x)
  GmuG[i]=logprime.cor(wt,ko,x)["pmuG"]
  x=par+c(0,0,i/nrep,0,0)
  dmuX[i]=loglik.cor(wt,ko,x)
  GmuX[i]=logprime.cor(wt,ko,x)["pmuX"]
  x=par+c(0,0,0,i/nrep,0)
  dsigmaG[i]=loglik.cor(wt,ko,x)
  GsigmaG[i]=logprime.cor(wt,ko,x)["psigmaG"]
  x=par+c(0,0,0,0,i/nrep)
  dsigmaX[i]=loglik.cor(wt,ko,x)
  GsigmaX[i]=logprime.cor(wt,ko,x)["psigmaX"]
}

test_that("Test the value of up-cor parameters asymptotically", {
  expect_equal(signif(Galpha[1],digits=2),signif((dalpha[1]-dalpha[2])*nrep,digits=2))
  expect_equal(signif(GmuG[1],digits=2),signif((dmuG[1]-dmuG[2])*nrep,digits=2))
  expect_equal(signif(GmuX[1],digits=2),signif((dmuX[1]-dmuX[2])*nrep,digits = 2))
})

