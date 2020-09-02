# 2 Marbles, 2 pots. 
# 1 step includes:
# - 1 marble randomly selected (prob of 1/2)
# - Then decide if the selected marble moves (1/3 prob it does move)
# - If it does move, randomly select where it moves to (1/3 prob for each of the 3 new pots)

# transition probability matrix.
library(expm)

probs <- matrix(rep(0,9),nrow=3,ncol=3)
probs[1,1] <- 2/3
probs[1,2] <- 1/3
probs[2,1] <- 1/6
probs[2,2] <- 2/3
probs[2,3] <- 1/6
probs[3,2] <- 1/3
probs[3,3] <- 2/3
probs

stationary <- matrix(rep(0,3),nrow=3,ncol=1)
stationary[1] <- 1/4 
stationary[2] <- 1/2
stationary[3] <- 1/4
stationary

probs %^% 20
probs %^% 2

X <- rep(1,2) # individual descriptions: initial state with all the marbles on the left section
Niter <- 100000 # number of iterations
K <- rep(0,Niter) # vector of states
Z <- 0 # auxiliary variable

for (i in 1:Niter) {
  k <- sample(2, size = 1) # selects a marble
  if (runif(1)<1/3){X[k] <- -X[k] + 1} else {X[k] <- X[k]} 
    # flips with probability 1/3
  Z = sum(X) # number of marbles on the left box
  K[i] <- Z # vector of states as a function of time
  Z <- 0 
  }



freq <- rep(0,3)
for (j in 1:3) {
  freq[j] <- length(which(K==j-1))/Niter #Calculating frequencies of states
}
freq

# Empirical vs invariant distribution

state = c(0:2) #vector of states

plot(state,freq,xlab="State",ylab="Probability of state", col="blue", main="Empirical and theoretical distributions")
lines(state,stationary,type="h")
legend("topright", c("Empirical distribution","Theoretical invariant distribution."), fill=c("blue", "black"), cex = 0.55)


######### Running Mean (Converges to theoretical mean as number of steps increase)

n <- 10 #Number of steps

Mean = seq(0,Niter/n-1) #Initialising vector of means

for (i in 1:Niter/n) {Mean[i]=mean(K[1:i])} #Running mean

m <- 1 # expected value of k (0.25 * 0 + 0.5 * 1 + 0.25 * 2) = 1

TMean = rep(m,Niter/n -1) #Vector of expected mean

plot (1:(Niter/n), Mean, xlab="Number of Steps", ylab="Mean State", main="Running mean state")
lines(TMean, col="red")
legend("topright", c("Mean States","Expected Value"), fill=c("black", "red"))

Mean
# Running variance (sampled every n Monte Carlo steps)

Var = seq(0,Niter/n-1) #Vector of variances

for (i in 1:Niter/n) {Var[i]=var(K[1:i])} #Running variance

v <- 0.5 # Theoretical variance of k 

TVar = rep(v,Niter/n-1)

plot(Var, xlab="Number of Steps", ylab="Variance", main = "Running variance")
lines(TVar, col = "red")
legend("bottomright", c("Mean States","Expected Value"), fill=c("black", "red"))

