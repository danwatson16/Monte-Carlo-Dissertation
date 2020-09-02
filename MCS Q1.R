### 3 algorithms to produce poisson random variables ### 

#Using exponential random variables to produce poisson random variables.

generate <- function(lambda, n) {
  x <- c()
  for (i in 1:n) {
  X <- 0
  sum <- 0
  while (TRUE) {
    exponential <- rexp(1, rate = 1) #Generating an exponential r.v
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

lambda_values = c(10, 50) #Vector of lambda values - the rates
sample_sizes <- c(1e+03, 1e+04, 1e+05, 1e+06) #Vector of sample sizes
chi_sq_values <- c()
for (lambda in lambda_values) {
  for (i in sample_sizes) {

result <- generate(lambda = lambda, i) #Generating r.vs with specificied lambda and sample size

k <- seq(0,max(result), length = max(result)+1)  #Theoretical poisson k values
y <- exp(-lambda)*((lambda**k)/factorial(k)) #Theoretical poisson pmf
hist(result, breaks = (max(result)-min(result)), probability = TRUE, main = paste("Histogram with sample size", i,"and rate", lambda, sep = " "))
lines(k,y) #Histogram of the random variables aswell as a theoretical curve

observed_vals <- c()
for (j in 0:max(result)) {
  observed_vals[j+1] <- length(which(result==j)) #Frequencies of each bin
}
y <- y*i #Theoretical frequencies

chi_sq <- sum((((observed_vals - y)**2))/y) #Chi-squared value

append(chi_sq_values, chi_sq)

df <- length(observed_vals) - 2 #number of bins - number of parameters estimating

chisq_table <- qchisq(p = 0.10, df = df)
Accept_H0 <- chi_sq <= chisq_table
print(Accept_H0)
#Null hypothesis that distriubitons are the same is 
#rejected with this sample size at significance 0.1.
  }
}

#Using unform rvs
generate1 <- function(lambda, n) {
  x <- c()
  for (i in 1:n) {
    X <- 0
    product <- 1
    while (TRUE) {
      uniform <- runif(1, min = 0, max = 1)
      product <- product*uniform
      if (product > exp(-lambda)) {
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
for (i in sample_sizes) {
result1 <- generate1(lambda = lambda, i)
k <- seq(0,max(result1), length = max(result1)+1) 
y <- exp(-lambda)*((lambda**k)/factorial(k))
hist(result1, breaks = max(result1), probability = TRUE)
lines(k,y)

observed_vals <- c()
for (j in 0:max(result1)) {
  observed_vals[j+1] <- length(which(result1==j)) #Frequencies of each bin
}
y <- y*i #Theoretical frequencies
chi_sq <- sum((((observed_vals - y)**2))/y) #Chi-squared value

df <- length(observed_vals) - 1 #number of bins - number of parameters estimating

chisq_table <- qchisq(p = 0.10, df = df)
Accept_H0 <- chi_sq <= chisq_table
print(Accept_H0)

print(chi_sq) #Null hypothesis that distriubitons are the same is usually rejected.
print(ks.test(result1, y))
print(chisq.test(observed_vals, y))

}
#Second method seems to produce more accurate distributions.