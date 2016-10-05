n=100
wtd <- rnorm(n=n,mean=-2,sd=0.8)
wtd <- cbind(wtd,2*wtd+rnorm(n=n,mean=1,sd=1.7))

par <- c(2,-2,1,log(0.8),log(1.7))
parup <- undownstream(par,to="up")

wtu <- rnorm(n=n,mean=parup[3],sd=exp(parup[5]))
wtu <- cbind(parup[1]*wtu+rnorm(n=n,mean=parup[2],sd=exp(parup[4])),wtu)

parcor <- undownstream(par,to="cor")

cor(wtd[,1],wtd[,2])
cor(wtu[,1],wtu[,2])


test_that("Equivalence of parameters", {
  expect_equivalent(parcor[1]*parcor[4]/parcor[5],parup[1])
  expect_equivalent(parcor[1]*parcor[5]/parcor[4],par[1])
  expect_equivalent(sqrt(parcor[4]^2-parup[1]^2*parcor[5]^2),parup[4])
  expect_equivalent(parup[1]*parup[3]+parup[2],parcor[2])
})
