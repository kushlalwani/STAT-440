### To Estimate the integral of u^2 as u goes from 0 to 1:
U <- runif(1e6)
mean(U^2)

### To Estimate the integral of exp(-x) as x goes from 0 to 1:
## Method A
X <- runif(1e5)
mean(exp(-X))
var(exp(-X))/1e5

## Method B
Y <- rexp(1e5)
mean(tmp <- Y<1)
var(Y<1) / 1e5

# WLLN demonstration
x <- rexp(1e5, rate=3)
plot(cumsum(x)/(1:1e5), type="l")
abline(h=1/3, lty=2, col=2)

### Tail Probability of Standard Normal (from Week 03b)
set.seed(999)

# naive simulation
X = rnorm(100000)
mean(X>100)




# Try importance sampling with a normal proposal "by hand" using z=10:
set.seed(999)
Y <- 10 + rnorm(1e4)
mean((Y>10) * dnorm(Y) / dnorm(Y-10))

# Try importance sampling with an exponential proposal "by hand" using z=10:
set.seed(999)
Y <- rexp(1e4) + 10 
mean ((Y>10) * dnorm(Y) / dexp(Y-10))
mean (exp(-Y^2/2)/sqrt(2*pi) / exp(-(Y-10)))

## Quickly work through Homework #2 Exercise 4(a):
X1 <- runif(1e4) * 3
mean(Y1 <- 3*X^2*exp(-X))

## Quickly work through Homework #2 Exercise 4(b):
X2 <- rgamma(1e4, shape=3)
mean(Y2 <- 2*(X<3))

## Simulation study comparing simple MC estimator to antithetic MC estimator:
set.seed(440)
cHat <- -0.16237
h <- function(x) dnorm(x)
n <- 1000
reps <- 750
simple <- 1:reps
antithetic <- 1:reps
control <- 1:reps
for (i in 1:reps) {
  U <- runif(2*n)
  simple[i] <- mean(h(U))
  # Next implement antithetic variate method
  V <- U[1:n] # only use the first n, so the comparison is valid
  antithetic[i] <- mean((h(V) + h(1-V))/2)
  control[i] <- mean(h(U) - cHat*(U^2 - 1/3))
}

U <- runif(1000)
plot(dnorm(U), U^2)

###########
# Week 05b stuff

x <- seq(0, 5, len=200)
plot(x, dchisq(x, 1), type="l")

## Test of Abbot again Roche, Example 3 from Week 05b
test_data = data.frame(
  testname=c("Roche", "Abbott", "MEDsan", "Siemens"),
  pos_tests=c(42, 38, 37, 40),
  neg_tests=c(42, 46, 47, 32)
)

sampls = rhyper(100000, 80, 88, 84)
mean(sampls <= 42)
phyper(42, 80, 88, 84)

fisher.test(as.matrix(test_data[1:2,c("pos_tests", "neg_tests")]),
            alternative="less")

chisq.test(as.matrix(test_data[1:2,c("pos_tests", "neg_tests")]))

####
#Week 06a: PSU vs Pitt

df = read.csv("https://tinyurl.com/PSUVsPittCSV")
MyHFATestStatistic <-
mean(df[df[,1]=="StateCollege",3]) - mean(df[df[,1]=="Pittsburgh",3])
MyHFATestStatistic

same_dist_perm_test = function( n_perms, xs, ys, test_statistic ) {
  #' @param n_perms number of permutations to generate
  #' @param xs vector of samples from distribution X
  #' @param yx vector of samples from distribution y
  #' @param test_statistic function that calculates the test statistic
  # calculate the number of samples in X and Y
  n = length(xs)
  m = length(ys)
  # define labels (1 = X samples, 0 = y samples)
  labels = c(rep(1, n), rep(0, m))
  all_data = c(xs, ys)
  # for every permutation replication
  replicate(
    n_perms, {
      # permute label orders
      permuted_labels = sample(labels)
      # generate new test statistic under permutation
      test_statistic(all_data[permuted_labels == 1],
                     all_data[permuted_labels == 0])
    }
  )
}


score_diff_perms = same_dist_perm_test(
  10000,
  df[df$Location == "StateCollege", "Diff"],
  df[df$Location == "Pittsburgh", "Diff"],
  function(a, b) { mean(a) - mean(b) }
)
hist(score_diff_perms)
abline(v=7.377)

###########
## Week 06b: Confidence intervals for the population median from 
## an exponential distribution
set.seed(440)
n <- 51
X <- rexp(n)
M <- median(X)
#####
# Nonparametric Bootstrapping approximates the distribution 
# of the sample median
# of a sample of size 51 from the UNKNOWN distribution by replacing
# that unknown distribution by the empirical distribution.
AllBootMedians <- 1:1e4
for (i in 1:1e4) {
  Xboot <- sample(X, 51, replace=TRUE)
  AllBootMedians[i] <- median(Xboot)
}
# 95% CI based on normality of sample median:
M + c(-1.96, 1.96) * sd(AllBootMedians)
# 95% CI based on quantiles of bootstrap sample:
quantile(AllBootMedians, c(0.025, 0.975))

# Parametric bootstrap based on the assumption that
# F is exponential(mean=mu):
AllParBootMedians <- 1:1e4
for (i in 1:1e4) {
  XPboot <- rexp(n, rate=1/mean(X))
  AllParBootMedians[i] <- median(XPboot)
}
# Conf int based on sample SD from parametric bootstrap values:
M + c(-1.96, 1.96) * sd(AllParBootMedians)
# Conf int based on sample quantiles from parametric bootstrap values:
quantile(AllParBootMedians, c(0.025, 0.975))

## Side question: What's the true SD of this median? Answer via simulation:
AllMedians <- 1:1e4
for (i in 1:1e4) {
  X <- rexp(n)
  AllMedians[i] <- median(X)
}
hist(AllMedians)
sd(AllMedians)

####### AN example of eigen-decomposition of a cov matrix:
set.seed(440)
X <- matrix(runif(20), 10, 2)
Sigma <- cov(X)
e <- eigen(Sigma)
M <- e$vec %*% diag(sqrt(e$val))

############### 
## Bayes problem: Prior on theta is Gamma(4,1)
## X is Poisson(theta)
## Posterior: Gamma(X + 4 , 2)
theta <- seq(0, 50, len=200) 
plot(theta, dgamma(theta, shape=4, rate=1), type="l")
lines(theta, dgamma(theta, shape=49, rate = 2), col=2, lty=2)

################
## Naive Bayes example (relevant for HW #5, Exercise 1)
path <- "https://raw.githubusercontent.com/data-8/textbook/main/assets/data"
wine <- read.csv(paste(path, "wine.csv", sep=""))
names(wine)
table(wine$Class)
pairs(wine[,2:5])
pairs(wine[,2:5], col=wine$Class)
pairs(wine[,c(2, 3, 5)], col=wine$Class)

#### Brute force approach:
colMeans(as.matrix(wine[c(2,3,5)])[wine$Class==1, ])
colMeans(as.matrix(wine[c(2,3,5)])[wine$Class==2, ])
colMeans(as.matrix(wine[c(2,3,5)])[wine$Class==3, ])

apply(as.matrix(wine[c(2,3,5)])[wine$Class==1, ] , MARGIN=2, FUN=var) 
apply(as.matrix(wine[c(2,3,5)])[wine$Class==2, ] , MARGIN=2, FUN=var) 
apply(as.matrix(wine[c(2,3,5)])[wine$Class==3, ] , MARGIN=2, FUN=var) 

### Pick a random row from the wine dataset:
wine[sample(nrow(wine), 1), c(2, 3, 5)]

### What about a prior distribution?
table(wine$Class) / nrow(wine)

### Explore the case where X ~ Bin(n,p) and the prior on p is beta:

# generate 1000 Bernoulli variables with (true) p=0.7
set.seed(440)
x_samps <- rbinom(1000, 1, .7)
phat <- mean(x_samps)

# Create a function to plot the posterior distribution
plot_posterior = function(Xdata, prior_a, prior_b) {
  n <- length(Xdata) 
  X <- sum(Xdata)
  xs <- seq(.001, .999, .001)
  plot(xs, dbeta(xs, X + prior_a, n - X + prior_b), type="l",
       main = "Dotted line is posterior mean")
  # The posterior Beta parameters are (X + prior_a, n - X + prior_b)
  # The posterior mean is therefore X + prior_a / (n + prior_a + prior(b))
  post_mean <- (X + prior_a) / (n + prior_a + prior_b)
  abline(v = post_mean, lty=2, col=2)
  axis(1, at = post_mean, lab=round(post_mean,2), 
       col.axis =2, col.ticks = 2)
}

plot_posterior(x_samps, 2, 2)

#####################
# Attempt at MCMC using a symmetric proposal distribution

f <- function (mu, sigsq) {# Take into account x=1.2
  dnorm(mu) * 1/sigsq^2 * exp(-1/sigsq) * (sigsq>0) *
  exp(-(1.2-mu)^2 / (2 * sigsq))
}

# Step 1: Initialization (remember, C mean current and P means proposed)
muC <- 1.2 ; sigsqC <- 1
set.seed(440)

AllTheMuValues <- rep(0, 1e4)
AllTheSigSqValues <- rep(0, 1e4)
for( i in 1:1e4) {
# Step 2: Generate proposed values
muP <- muC + rnorm(1)/10 ; sigsqP <- sigsqC + rnorm(1)/10
# Step 3: Calculate the acceptance ratio
A <- f(muP, sigsqP) / f(muC, sigsqC)
# Steps 4 and 5: Generate U~unif(0,1) and compare with A
if ( (U<-runif(1)) < A ) { muC <- muP ; sigsqC <- sigsqP}
# Step 6: Save the values
AllTheMuValues[i] <- muC ; AllTheSigSqValues[i] <- sigsqC
}

################
## Naive Bayes example (relevant for HW #5, Exercise 1)
path <- "https://raw.githubusercontent.com/data-8/textbook/main/assets/data"
wine <- read.csv(paste(path, "wine.csv", sep=""))
names(wine)
table(wine$Class)
pairs(wine[,c(2, 3, 5)], col=wine$Class)

# Three-fold cross-validation: break dataset into three groups:
set.seed(440)
GroupA <- sample(1:178, 59)
GroupB <- sample((1:178)[-GroupA], 59)
GroupC <- (1:178)[-c(GroupA, GroupB)]

#### Brute force approach for "training" (finding means and variances):
#Means for Groups A and B, Class 1 through 3
colMeans(as.matrix(wine[c(GroupA, GroupC), c(2, 3, 5)])[wine[c(GroupA, GroupC),]$Class==1, ])
colMeans(as.matrix(wine[c(GroupA, GroupC), c(2, 3, 5)])[wine[c(GroupA, GroupC),]$Class==2, ])
colMeans(as.matrix(wine[c(GroupA, GroupC), c(2, 3, 5)])[wine[c(GroupA, GroupC),]$Class==3, ])

### Doing this a bit more systematically...
### Part (b) asks for nine means and nine variances for the three variables
### In my case, this is columns 2, 3, and 5; yours will be different
Variables <- as.matrix(wine[, c(2, 3, 5)]) 

### As an alternative to colMeans, we can use "apply" as follows:
MuHat <- NULL
SigmaHat <- NULL
for (category in 1:3) {
  subset <- wine$Class == category
  MuHat <- rbind(MuHat, apply (Variables[subset,], MARGIN = 2, FUN = mean))
  SigmaHat <- rbind(SigmaHat, apply (Variables[subset,], MARGIN = 2, FUN = sd))
}
rownames(MuHat) <- rownames(SigmaHat) <- paste("Class", 1:3)



### The remainder of the code assumes that we are 
### using GroupA and GroupC as the training rows and GroupB as the testing rows
### as in one section of part (e)
MuHat <- NULL
SigmaHat <- NULL
for (category in 1:3) {
  # Choose the subset of rows that are in the right category AND in Group A or B
  subset <- (wine$Class == category) & ((1:178) %in% c(GroupA, GroupC))
  MuHat <- rbind(MuHat, apply (Variables[subset,], MARGIN = 2, FUN = mean))
  SigmaHat <- rbind(SigmaHat, apply (Variables[subset,], MARGIN = 2, FUN = var))
}
rownames(MuHat) <- rownames(SigmaHat) <- paste("Class", 1:3)

### For implementing the formulas in part (c), here's a function that takes as inputs:
###   -- a 3-vector of X values
###   -- a 3x3 matrix of mu values, one mean per every variable and category
###   -- a 3x3 matrix of sig^2 values, one variance per every variable and category
### and outputs three values, namely:
###   -- P(X1, X2, X3 | C=j) for j=1, 2, 3
JointProbs <- function(x, mu, sigma) {
  ans <- rep(0,3) 
  for (j in 1:3) {
    ans[j] <- prod(dnorm(x, mean=mu[j,], sd=sigma[j,]))
  } 
  return(ans)
}

### For part (d), you'll need to figure out how to use the output from JointProbs
### to calculate the conditional probabilities. 
Probs <- matrix(0, length(GroupB), 3)
colnames(Probs) <- paste("Class", 1:3)
for (i in 1:length(GroupB)) {
  JP <- JointProbs(Variables[GroupB,][i,], MuHat, SigmaHat)
  Probs[i, ] <- JP / sum(JP)
}

### If you only need to know which category is the most likely, as in the
### second part of (d), then you can just find the maximum value
### returned by that function, is in the code below:
Tests <- matrix(0, length(GroupB), 2)
colnames(Tests) <- c("True", "Predicted")
for (i in 1:length(GroupB)) {
  Tests[i, 1] <- wine$Class[GroupB][i]
  Tests[i, 2] <- which.max(JointProbs(Variables[GroupB,][i,], MuHat, SigmaHat))
}

### There are several ways to summarize the test results. 
### One is to simply report in a table form which are the correct and which are
### the predicted categories:
table(Tests[,1], Tests[,2]) # Rows = True; Columns = Predicted



## Implement Newton-Raphson to find a zero of f(x) = x^2 - 15
xNow <- 20

f <- function(x) x^2-15
fPrime <- function (x) 2*x
xNext <- xNow - f(xNow) / fPrime(xNow)

print(sqrt(15) - (xNext <- xNext - f(xNext) / fPrime(xNext)))



### Finding an MLE for theta if X ~ Binomial(5, theta)
### This method will use Newton's method

X <- 3 # Let's assume that this is the value we observed (so the MLE equals 3/5)

LogLik <- function(theta) {
  log(choose(5,X)) + X*log(theta) + (5-X)*log(1-theta) 
}

# Plot the loglik for theta in (0,1):
thetaSeq <- seq(0,1,len=200)
plot(thetaSeq, LogLik(thetaSeq), type="l", ylim=c(-12,0))

# We also need first and second derivatives to search for solution to lPrime=0
lPrime <- function(theta) {
  X/theta - (5-X)/(1-theta)
}
lDblPrime <- function(theta) {
  -X/theta^2 - (5-X)/(1-theta)^2
}

# Initialize value of thetaNow:
thetaNow <- 0.3

# Mark position of thetaNow, then update it using one Newton step:
abline(v = thetaNow, lty=2, col=2)
thetaNow <- thetaNow - lPrime(thetaNow) / lDblPrime(thetaNow)


# Use optim function to do the same maximization problem:
optim(0.3, LogLik, control = list(fnscale = -1))


### Use Newton-Raphson to find an MLE for a sample from
### Beta(theta1, theta2)

# Read in the dataset (copied from HW #7 Exercise 2)
mlbbat2010 <- read.csv("https://www.openintro.org/data/csv/mlbbat10.csv")
ba <- mlbbat2010$bat_avg
ab <- mlbbat2010$at_bat
x <- ba[ab >= 50]

# Loglikelihood function (copied from HW #7 Exercise 1)
ell <- function(theta, data = x) {
a <- theta[1]
b <- theta[2]
n <- length(data)
return( n * lgamma(a+b) - n * lgamma(a) - n * lgamma(b) +
          (a-1) * sum(log(data)) + (b-1)* sum(log(1 - data)))
}

# Test loglikelihood function on theta=(20,60)
theta <- c(20, 60)
ell(theta)

# Define gradient-vector function using formulas derived in class:
gradient <- function(theta, data = x) {
  a <- theta[1]
  b <- theta[2]
  n <- length(data)
  return( c(n * digamma(a+b) - n * digamma(a) + sum(log(data)),  
    n * digamma(a+b) - n * digamma(b) + sum(log(1-data)))
  )
}

gradient(theta) # Test to make sure this gives a 2-vector

# Define Hessian-matrix function using formulas derived in class:
Hessian <- function(theta, data = x) {
  a <- theta[1]
  b <- theta[2]
  n <- length(data)
  return( n * trigamma(a+b) - n*diag(c(trigamma(a), trigamma(b))))
}

Hessian(theta) # Test to make sure this gives a 2x2 matrix

# Now we can implement Newton-Raphson:
print(theta <- theta - solve(Hessian(theta)) %*% gradient(theta))

### Sample from a mixture of two gamma densities (Week 13a notes)
n <- 10000 # number of samples
x <- rep(0,n) # to store the samples
shape <- c(2,2) # shapes of the two components
scale <- c(0.5,1) # scales of the two components
for(i in 1:n){
  if(runif(1)<0.7)
    z=1 else z=2
    x[i] = rgamma(1,scale=scale[z],shape=shape[z])
}

hist(x, nclass=20, prob=TRUE)
xseq <- seq(0,10,len=200)
lines(xseq, dgamma(xseq, 2, scale=0.5), col=2, lwd=3)
lines(xseq, dgamma(xseq, 2, scale=1), col=3, lwd=3)
lines(xseq, 0.7* dgamma(xseq, 2, scale=0.5) + 0.3*dgamma(xseq, 2, scale=1), col=4,lwd=3)

### Build the log-likelihood:
mix_loglik = function(pi1, data = x){
  sum(log(pi1 * dgamma(data ,scale=0.5, shape=2) +
            (1-pi1)*dgamma(data, scale=1, shape=2)))
}

pi1seq <- seq(0, 1, len=200)
LogLikPi1 <- rep(0,200) 
for(i in 1:200) {
  LogLikPi1[i] <- mix_loglik(pi1seq[i])
}
plot(pi1seq, LogLikPi1, type="l")
abline(v=0.7, col=2, lty=2)

# "Eyeball" method of finding the MLE:
pi1seq <- seq(0.7096, 0.7104, len=200)
LogLikPi1 <- rep(0,200) 
for(i in 1:200) {
  LogLikPi1[i] <- mix_loglik(pi1seq[i])
}
plot(pi1seq, LogLikPi1, type="l")
abline(v=0.7, col=2, lty=2)

piNow <- 0.6
# Now repeat the following lines over and over until convergence:
abline(v=piNow, lty=2, col=4)
w <- dgamma(x, scale=0.5, shape=2) * piNow /
   (dgamma(x, scale=0.5, shape=2) * piNow + dgamma(x, scale=1, shape=2) * (1-piNow))
piNow <- mean(w)



#### Game of Thrones Logisitic Regression Example from Week 13b:
GoT = read.csv("https://sites.psu.edu/drh20/files/2023/11/got_characters.csv")
head(GoT)


## Shortcut for fitting logistic regression model using R:
got_mod = glm(
  isAlive ~ isMale + isNoble + popularity + house,
  data=GoT, family=binomial(link="logit")
)

### Set up Newton-Raphson
Y <- GoT$isAlive # Response (True/False)

## To Be Continued...

## Need to create a lot of indicator covariates fromm column 2
X <- matrix(GoT[,c(2:5)])


#### Game of Thrones Logistic Regression Example from Week 13b:
GoT = read.csv("https://sites.psu.edu/drh20/files/2023/11/got_characters.csv")
head(GoT)

## Shortcut for fitting logistic regression model using R:
got_mod = glm(
  isAlive ~ isMale + isNoble + popularity + house,
  data=GoT, family=binomial(link="logit")
)

### Set up Newton-Raphson in several steps. First, need  
### X (nxp matrix of covariates) and Y (nx1) vector of responses
Y <- GoT$isAlive # Response (True/False)

## For X, create a lot of indicator covariates from the categorical column 'house'.
## Use the function directly from the Week 13b notes for this:
house_xs = sapply(
  # names of all houses except the first
  # one (which is "other")
  names(table(GoT$house))[-1],
  # for each house name
  function(house_name) {
    # create a new column based on whether
    # the characters are in that house
    as.integer(GoT$house == house_name)
  }
)
## Now we can create the X matrix:
X <- cbind(Intercept=1, Male = GoT$isMale, Noble = GoT$isNoble, 
           Pop = GoT$popularity, house_xs) 
## Starting beta value of all zeros. My using the 1st row of X we can
## give beta the same set of variable labels:
beta <- 0*X[1,] 
names(beta)

## Now for the Newton-Raphson. Calculate the gradient and Hessian
## using the Week 13b, page 2 formulas:
theta <- exp(X %*% beta) / (1 + exp( X %*% beta))
gradient <- t(X) %*% (Y - theta)
## For the Hessian, it's wasteful to calculate a HUGE diagonal matrix of mostly
## zeros, so we can accomplish the same thing using the sweep function in R:
Hessian <- -t(X) %*% sweep(X, MARGIN = 1, STATS = theta*(1-theta), FUN = '*')

## Finally, the Newton update to beta:
beta <- print(beta - solve(Hessian) %*% gradient)

## Repeat the lines above several times. 
## Check the difference between the R-fitted model and the Newton-Raphson solution:
coef(got_mod) - beta



### PCA dataset (from Week 15a notes)
df <- read.csv("https://sites.psu.edu/drh20/files/2023/11/video_game_lengths.csv")
head(df)
names(df)
dim(df)
X <- as.matrix(df[,4:9])
dim(X)
t(X) %*% X
cm <- colMeans(X)
X <- sweep(X, MARGIN=2, STATS=cm, FUN='-')
# TopHat question: How are X^TX and cov(X) related?
(t(X) %*% X) / cov(X)
S <- cov(X)
# In parallel with Exercise 2(a) on HW#8:
e <- eigen(S)
sv <- svd(X)
# Check that the eigenvectors are orthogonal:
round(e$vectors %*% t(e$vectors),10)
# Check that the V matrix of the SVD is orthogonal:
round(sv$v %*% t(sv$v),10)
# Compare the V matrix of the SVD with the eigenvectors:
# First column:
cbind(e$vec[,1], sv$v[,1])
# Second column:
cbind(e$vec[,2], sv$v[,2])
# Third column:
cbind(e$vec[,3], sv$v[,3])
# Fourth column:
cbind(e$vec[,4], sv$v[,4])
### ...and so on.
## Conclusion: sv$v and e$vec are the same except in some cases their
## signs (of particular columns) are opposite
# Compare the eigenvalues to the singular values:
rbind(e$val, (sv$d)^2/1211)
# Conclusion: Squaring the singular values of X and dividing by n-1
# gives the eigenvalues of t(X) %*% X / (n-1)

### Let's do PCA using the eigendecomposition:
## The first principal component is the eigenvector with the largest eigenvalue:
e$vectors[,1]

# The proportion of variability explained by this PC is its e-val divided by
# the sum of the e-values:
e$values[1] / sum(e$values)
# We can write "The first PC explains 92.3% of the variation in those 6 variables."
# Based on the values of that first column of V (i.e., the first e-vector),
# The first PC is roughly the average of the three "Complete" variables.

# What about the 2nd PC? Look at the second column of the e-vector matrix:
e$vectors[,2]
# It explains this fraction of the variability:
e$values[2] / sum(e$values)
# Based on the loadings (i.e., the coefficients in the linear combo expressed
# by the 2nd PC), it looks like that 2nd PC is a contrast between var 4 and
# variables 5&6.

# So far, we can explain over 97% of the variability and we haven't needed those
# first three variables. What's up?
# Let's take a look at the covariances:
round(S, 0)

# Due to the large difference in the scale of variation between "MainStory" and "Complete" 
# variables, let's try a PCA on the correlation matrix:
S2 <- cor(X)
e2 <- eigen(S2)
# The first PC for the correlation structure:
e2$vectors[,1]
# How much does it explain?
e2$values / sum(e2$values)
# Clearly, the first PC is the sample mean of all six variables. It explains 76%
# of the total variation in these (standardized) data.

