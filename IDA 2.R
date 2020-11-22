library(maxLik)

load("dataex2.Rdata")
dataex2
#optim
ll <- function(x, r, mu){
  sigma = 1.5
  -sum(r * dnorm(x, mu, sigma, log = TRUE) + (1-r) * pnorm(x, mu, sigma, log = TRUE))
}
opt <- optim(8, ll, x = dataex2$X, r = dataex2$R, method = "BFGS")

#maxLik
ll <- function(mu){
  sigma = 1.5
  sum(r * dnorm(x, mu, sigma, log = TRUE) + (1-r) * pnorm(x, mu, sigma, log = TRUE))
}
require(maxLik)
x = dataex2$X; r = dataex2$R
mle <- maxLik(logLik = ll, start = c(5))

mugrid <- seq(0,11,len = 200)
n <- length(mugrid)
res <- numeric(n)
for (i in 1:n){
  res[i] <- ll(mu = mugrid[i])
}
plot(mugrid, res, type = "l", xlab = expression(mu), ylab = "log likelihood")
abline(v = mle$estimate, col = "red")



load("dataex4.Rdata")
dataex4

mis <- which(is.na(dataex4$Y))
xm = dataex4$X[-mis]
xn = dataex4$X[mis]
ym = dataex4$Y[-mis]
yn = dataex4$Y[mis]

p <- function(x,beta){
  beta0 = beta[1]
  beta1 = beta[2]
  exp(beta0 + x*beta1)/(1 + exp(beta0 + x*beta1))
}
Q <- function(beta){
  beta0 = beta[1]
  beta1 = beta[2]
  s1 = sum(ym*(beta0 + xm*beta1))
  s2 = -sum(log(1 + exp(beta0 + dataex4$X*beta1)))
  s3 = sum((beta0 + xn * beta1)*p(xn, beta.old))
  -(s1 + s2 + s3)
}

diff = 1
eps = 0.00001
i = 1
beta.old = c(1,2)
while(diff > eps){
  opt <- optim(c(1,2),Q)
  diff = sqrt(sum((beta.old - opt$par)^2))
  beta.old = opt$par
  i = i + 1
}
opt



load('dataex5.Rdata')
dataex5
em <- function(y, theta0, eps){
  n = length(y)
  theta <- theta0
  p <- theta[1]
  mu <- theta[2]
  sigma <- theta[3]
  lambda <- theta[4]
  diff <- 1
  while (diff > eps){
    theta.old <- theta
    
    #E-step
    ptilde1 <- p*dlnorm(y, mu, sigma)
    ptilde2 <- (1-p)*dexp(y, lambda)
    ptilde <- ptilde1/(ptilde1 + ptilde2)
    
    #M-step
    p <- mean(ptilde)
    mu <- sum(log(y)*ptilde)/sum(ptilde)
    sigma <- sqrt(sum((log(y)-mu)^2*ptilde)/sum(ptilde))
    lambda <- sum(1-ptilde)/sum((1-ptilde)*y)
    theta <- c(p, mu, sigma, lambda)
    diff <- sum(abs(theta-theta.old))
  }
  return(theta)
}

res <- em(y = dataex5, theta0 = c(0.1, 1, 0.5, 2), eps = 0.00001)
pest <- res[1]
muest <- res[2]
sigma2est <- res[3]^2
lambdaest <- res[4]
pest; muest; sigma2est;lambdaest


hist(dataex5,freq = FALSE, xlab = "Sample",
     ylab = "Density", col = "yellow",cex.main = 1.5, cex.lab = 1.5,ylim = c(0,0.2), cex.axis = 1.4, main = "Histogram of data with estimated density superimposed.")
curve(pest*dlnorm(x, muest, sigma2est) + (1-pest)*dexp(x, lambdaest), add = TRUE, lwd = 2, col = "red")
axis(side = 4)



