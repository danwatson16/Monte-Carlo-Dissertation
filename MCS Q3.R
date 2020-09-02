
kTJ = 1.9 # temperature
betaJ = 1/kTJ # inverse temperature
Niter = 100000 # number of Monte Carlo steps

s <- rep(0,5) #Vector of states

for (i in 1:5) {s[i]= 1} # initial condition: all spins down
H = -betaJ*(s[1]*s[2]+s[2]*s[3]+s[3]*s[4]+s[4]*s[1] + s[1]*s[5] + s[2]*s[5] + s[3]*s[5] + s[4]*s[5]) # initial energy

E <-rep(0,Niter+1) #Vector of energy
M <-rep(0,Niter+1) #Vector of magnetisation
E[1] <- H #initial energy
M[1] <- sum(s) #initial magnetization

for (i in 1:Niter) {
  k <- sample(5, size = 1) # selects a spin
  Hs = -(s[1]*s[2]+s[2]*s[3]+s[3]*s[4]+s[4]*s[1] + s[1]*s[5] + s[2]*s[5] + s[3]*s[5] + s[4]*s[5]) #old hamiltonian
  s[k]<- -s[k] # the spin is flipped
  Hy = -(s[1]*s[2]+s[2]*s[3]+s[3]*s[4]+s[4]*s[1] + s[1]*s[5] + s[2]*s[5] + s[3]*s[5] + s[4]*s[5]) #new hamiltonian
  DeltaH <- -betaJ*(Hy - Hs) 
  alpha <- min(1,exp(DeltaH)) #Generate a random uniform number in [0,1]. If less then alpha, accept otherwise reject.
  U <- runif(1)
  if (U<=alpha) #the move is accepted, energy and magnetization are updated
  {E[i+1]<-E[i]+DeltaH
  M[i+1]<-sum(s)}
  else #The move is and the new energies and magnetisations are the same as the old.
  {s[k] <- -s[k]
  E[i+1]<- E[i]
  M[i+1]<- M[i]}
}


# Average absolute magnetization per spin

memp <- mean(abs(M[5000:10000]))/5

# Running average (sampled every n Monte Carlo steps)
n = 10

Mean = seq(0,Niter/n-1)

for (i in 1:Niter/n) {Mean[i]=mean(abs(M[1:i]))/5} #Calculating the mean

m <- 0.8889 # expected value of m

TMean = rep(m,Niter/n -1) #Theoretical mean

plot (Mean,xlab="MC steps per spin",ylab="|m| (kT/J = 1.9)", main = "Running average vs expected value")
lines(TMean, col = "blue")
legend("topright", c("Expected value","Running mean."), fill=c("blue", "black"))

# Plot of results
kTJ <- c(0.1,0.3,0.5,0.7,0.9,1,1.1,1.3,1.5,1.7,1.9)
memp <- c(1,1,0.999,0.999,0.997,0.995,0.990,0.976,0.954,0.924,0.889)
m <- (2*exp(8/kTJ) + (6/5) + (8/5)*exp(-2/kTJ) + (8/5) + (24/5)*exp(2/kTJ) + (4/5)*exp(-4/kTJ)) / ( 16*cosh(2/kTJ) + 10 + 2*exp(8/kTJ) + 4*exp(-4/kTJ) )
plot(kTJ,memp,xlab="kT/J",ylab="|m|",col="red", main = "The effect of increase in temperature on magnetisation")
lines(kTJ,m)
legend("topright", c("Expected value","Empirical mean"), fill=c("black", "red"))

#Energy
energy <- (32*betaJ*sinh(2/kTJ) + 16*betaJ*exp(8/kTJ) - 16*betaJ*exp(-4/kTJ))/(16*cosh(2/kTJ) + 10 + 2*exp(8/kTJ) + 4*exp(-4/kTJ))/5
plot(energy, xlab = "kT/J", ylab = "Energy", main = "Energy change as temperature increases")

MeanEnergy <- seq(0, Niter/n-1)
for (i in 1:Niter/n) {MeanEnergy[i]=mean(E[1:i]/5)} #Running energy mean
plot(MeanEnergy, main = "Running energy", xlab = "Energy", ylab = "MC steps per spin")
