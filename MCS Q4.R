#Generates poisson variable from exponential variable.
library(moments)
library(tseries)
generate <- function(lambda, n) {
  x <- c()
  for (i in 1:n) {
    X <- 0
    sum <- 0
    while (TRUE) {
      exponential <- rexp(1, rate = 1)
      sum <- sum + exponential
      if (sum < lambda) {
        X <- X + 1
      }
      else {
        x[i] <- X 
        break
      }
    }
  }
  return (x)
}

#Summing random variables - random walk:
n <- 10000
Ns <- c(10, 100, 1000, 5000)
lambda = 1
for (N in Ns) {

ZN <- rep(NA,n)

for (i in 1:n) {
  ZN[i] <- sum(generate(lambda = lambda, N))
} 

x <- seq(0, max(ZN), length = max(ZN) + 1) 
FZN <- dpois(x, lambda = lambda*N)
hist(ZN, breaks = max(ZN) - min(ZN) , probability = TRUE, main = paste("Histogram with N value: ", N, sep = " "))
lines(FZN)


CDF <- rep(NA,n)
for (i in 1:(n)) {
  CDF[i] <- sum(FZN[1:i])
}
ZNcdf = ecdf(ZN)
plot(ZNcdf, col = "red", lty = 1, main = paste("Theoretical vs Empirical CDFs with N value: ", N, sep = " "))
lines(CDF, col = "blue", lty = 2)
legend(x = "right", legend = c("Fit", "Observed"), col = c("red", "blue"), lty = 1:2)



expectation <- lambda
variance <- lambda 

UN <- (ZN - expectation*N)/sqrt(variance*N)
normSeq <- seq(-10, 10, 0.1)
y = dnorm(normSeq)
hist(UN, breaks = "FD", probability = TRUE, main = paste("Histogram of UN with N value:", N, sep = " "))
lines(normSeq, y)

print(jarque.bera.test(UN))
lr <- diff(UN)
mlr <- mean(lr)
stdlr <- sd(lr)
JBstat <- rep(NA,length(UN))
number <- length(lr)
for (i in 1:n) {
  lrth <- rnorm(number,mean=mlr,sd=stdlr)
  JBstat[i] <- number*(skewness(lrth)^2+(kurtosis(lrth)-3)^2/4)/6
}
hist(JBstat, main = paste("Histogram of JB statistic with N value: ", N, sep = " "))
}

