## -----------------------------------------------------------------------------
lmR <- function(y,X){
  if(length(y)!=nrow(X)){
    print("The matrix is not inversable！")
    return()
  }
  else return(solve(t(X)%*%X)%*%t(X)%*%y)
}

## -----------------------------------------------------------------------------
x <- 1:5
y <- c(4,6,9,2,5)
plot(x,y)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1000
u <- runif(n)
x <- ifelse(u<=1/2, log(2*u), -log(2*(1-u)))
hist(x, prob = TRUE, xlim = c(-6,6), ylim = c(0,0.5))
#lines(density(x))
y <- seq(-6, 6, 0.01)
lines(y, 1/2*exp(-abs(y)))

## -----------------------------------------------------------------------------
mybeta <- function(a, b, n){
  j<-k<-0
  y <- numeric(n)
  while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1) #random variate from g(.)
    if (x^(a-1) * (1-x)^(b-1) > u) {
      #we accept x
      k <- k + 1
      y[k] <- x
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
a <- 3
b <- 2
n <- 1000
set.seed(123)
x <- mybeta(a, b, n)
hist(x, prob = TRUE, xlim = c(0,1), ylim = c(0,2))
y <- seq(0, 1, 0.01)
lines(y, 1/beta(a, b)*y^(a-1)*(1-y)^(b-1))

## -----------------------------------------------------------------------------
myEpanechnikov <- function(n){
  x <- numeric(n)
  for(i in 1:n){
      u <- runif(3, min = -1, max = 1)
    if(abs(u[3])>=abs(u[2]) & abs(u[3])>=abs(u[1])) x[i] <- u[2]
    else x[i] <- u[3]
  }
return(x)
}

## -----------------------------------------------------------------------------
n <- 1000
#set.seed(123)
x <- myEpanechnikov(n)
hist(x, prob = TRUE, xlim = c(-1,1), ylim = c(0,1))
lines(density(x))

## ----eval=FALSE, include=FALSE------------------------------------------------
#  ## The acceptance-rejection method can not be applied to this problem, mainly because we generate a continuous random variable X=x from g and then we check whether f(x)/cg(x)>U. If so, we accept this x. Here since the objective random variable is discrete, we can not map a continuous one to a discrete one one by one. Thus even we get X, we can not get X from f.
#  
#  mysample <- function (vec, size, prob = NULL) {
#    l <- length(vec)
#    if (is.null(prob)) prob <- rep(1/l, l)
#    sort.vec <- sort(vec)
#    sort.prob <- prob[order(vec)]
#    j <- k <- 0
#    y <- numeric(size)
#    while (k < size) {
#      u <- runif(1)
#      j <- j + 1
#      x <- runif(1) #random variate from g(.)
#      temp <- which(sort.vec >= x)
#      if (length(temp) != 0) {
#        temp.prob <- sort.prob[length(temp)]
#        #find the probability of the smallest value no less than x
#        if (temp.prob > u) {
#          #we accept x
#          k <- k + 1
#          y[k] <- sort.vec[length(temp)]
#        }
#      }
#      # if length(temp) == 0, we drop this and do another, because we can not map x to a value from f. The map is on the left of each point value.
#    }
#    return(y)
#  }

## -----------------------------------------------------------------------------
mysample <- function (x, size, prob = rep(1/length(x), length(x))) {
  sort.x <- sort(x)
  sort.prob <- prob[order(x)]
  cumprob <- cumsum(sort.prob)
  u <- runif(size)
  r <- findInterval(u, cumprob) + 1 # the index of the interval is from 0 to length(x), each interval is on the left of each point in x and concluding x.
  return(sort.x[r])
}

mysample(1:10, 3)


## -----------------------------------------------------------------------------
set.seed(12345)
l <- c(0.1,0.5,1)
d <- 1
n <- 1e6
K <- 100
pi_hat <- numeric(0)
for(k in 1:K) {
  X <- runif(n,0,d/2)
  Y <- runif(n,0,pi/2)
  pihat <- numeric(3)
  for(i in 1:3) {
    pihat[i] <- 2*l[i]/d/mean(l[i]/2*sin(Y)>X)
  }
  pi_hat <- rbind(pi_hat, pihat)
}

# the value of pi_hat
apply(pi_hat, 2, mean)

# asymptotic variance of pi_hat
apply(pi_hat, 2, var)/(K-1)

## -----------------------------------------------------------------------------
# the simple Monte Carlo method
m <- 1e6
set.seed(123)
u <- runif(m)
mc_s <- exp(u)
thetahat_s <- mean(mc_s)
v_s <- var(mc_s)/(m-1)

# the antithetic variate approach
m <- 1e6
set.seed(123)
u <- runif(m/2)
mc_a <- (exp(u)+exp(1-u))/2
thetahat_a <- mean(mc_a)
v_a <- var(mc_a)/(m/2-1)

# an empirical estimate of the percent reduction in variance using the antithetic variate
(v_s-v_a)/v_s

## -----------------------------------------------------------------------------
set.seed(123)
m <- 1000

# phi_1 inverse transformation
f1 <- function(x){
  sqrt(2*pi*exp(1))/x
}
IF1 <- function(u){
  sqrt(1-2*log(1-u))
}
u <- runif(m)
x <- IF1(u)
estimate1 <- mean(1/f1(x))
var1 <- var(1/f1(x))

# phi_2 acceptance-rejection algorithm
f2 <- function(x){
  1/x^2/(1-pnorm(1))
}
j <- k <- 0
y <- numeric(m)
while(k<m){
  u <- runif(1)
  j <- j + 1
  x <- rnorm(1) #random variate from g(.)
  if (ifelse(x>1,1,0)>u) {
    #accept x
    k <- k + 1
    y[k] <- x
  }
}
estimate2 <- mean(1/f2(y))
var2 <- var(1/f2(y))

# theoretical integral
g <- function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)
}
value <- integrate(g,1,Inf)$value
cat(value, estimate1, estimate2)
cat(var1, var2)

## -----------------------------------------------------------------------------
# stratified importance sampling estimate
M <- 20
k <- 5
m <- M/k
set.seed(123)
f <- list()
IF <- list()
gf <- list()
a <- seq(0,1,by=1/k)
g <- function(x){
  exp(-x)/(1+x^2)
}
estimate <- numeric(k)
var <- numeric(k)
for(j in 1:5){
  f[[j]] <- function(x){
    exp(-x)/(exp(-a[j])-exp(-a[j+1]))
  }
  IF[[j]] <- function(u){
    -log(exp(-a[j])-(exp(-a[j])-exp(-a[j+1]))*u)
  }
  gf[[j]] <- function(x){
    (exp(-a[j])-exp(-a[j+1]))/(1+x^2)
  }
  u <- runif(m)
  estimate[j] <- mean(gf[[j]](IF[[j]](u)))
  var[j] <- var(gf[[j]](IF[[j]](u)))
}
cat(integrate(g,0,1)$value, sum(estimate), sum(var))
cat((0.09658794^2-8.817495e-05)/0.09658794^2)

## -----------------------------------------------------------------------------
set.seed(123)
m <- 1000
n <- 20
alpha <- 0.05
sample <- replicate(m, expr = {
x <- rchisq(n,2)
} )
c1 <- apply(sample,2,mean)+qt(alpha/2,n-1)*apply(sample,2,sd)/sqrt(n)
c2 <- apply(sample,2,mean)-qt(alpha/2,n-1)*apply(sample,2,sd)/sqrt(n)
mean(ifelse(c1<2 & c2>2,1,0))              

## -----------------------------------------------------------------------------
set.seed(123)
m <- 1000
n <- 20
alpha <- 0.05
sample1 <- replicate(m, expr = {
  x <- rchisq(n,1)
} )
sample2 <- replicate(m, expr = {
  x <- runif(n,0,2)
} )
sample3 <- replicate(m, expr = {
  x <- rexp(n)
} )

test <- function(sample,mu){
  c1 <- mu+qt(alpha/2,n-1)*apply(sample,2,sd)/sqrt(n)
  c2 <- mu+qt(1-alpha/2,n-1)*apply(sample,2,sd)/sqrt(n)
  xbar <- apply(sample,2,mean)
  mean(ifelse(xbar<c1|xbar>c2,1,0))
}
test(sample1,1)
test(sample2,1)
test(sample3,1)

## -----------------------------------------------------------------------------
library(foreach)
hypothesis <- function(method){
set.seed(10)
estimate <- foreach(i = 1:1000, .combine = "rbind") %do% {
  p0 <- runif(950)
  p1 <- rbeta(50,0.1,1)
  p <- c(p0,p1)
  p0.adjust <- p.adjust(p, method = method)[1:950]
  p1.adjust <- p.adjust(p, method = method)[951:1000]
  V <- sum(p0.adjust<0.1)
  S <- sum(p1.adjust<0.1)
  R <- V+S
  result <- c(V>=1, V/R, S/50)
  names(result) <- c('FWER','FDR','TPR')
  result
}
apply(estimate,2,mean)
}
hypothesis('bonferroni')
hypothesis('BH')

## -----------------------------------------------------------------------------
result <- function(n){
set.seed(123)
m <- 100
B <- 100
simulation <- foreach(i=1:m, .combine = 'rbind') %do% {
  x <- rexp(n,rate=2)
  lambda.hat <- 1/mean(x)
  bootstrap <- foreach(i=1:B, .combine = 'c') %do% {
    b <- sample(1:n, size = n, replace = TRUE)
    1/mean(x[b])
  }
  bias.bootstrap <- mean(bootstrap)-lambda.hat
  standard.error.bootstrap <- sd(bootstrap)
  temp <- c(bias.bootstrap,standard.error.bootstrap)
  names(temp) <- c('bias.bootstrap','standard.error.bootstrap')
  temp
}
cat('theoretical results:',2/(n-1),2*n/(n-1)/sqrt(n-2),'\n')
apply(simulation, 2, mean)
}

result(5)
result(10)
result(20)


## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(123)
B <- 200
n <- nrow(law)
alpha <- 0.05
R.hat <- cor(law$LSAT,law$GPA)
bootstrap <- foreach(i = 1:B, .combine = 'c') %do% {
  b <- sample(1:n, size = n, replace = TRUE)
  LSAT <- law$LSAT[b]
  GPA <- law$GPA[b]
  R.hat.b <- cor(LSAT, GPA)
  R.hat.b.k <- foreach(j = 1:B, .combine = 'c') %do% {
    k <- sample(1:n, size = n, replace = TRUE)
    LSAT <- law$LSAT[k]
    GPA <- law$GPA[k]
    cor(LSAT, GPA)
  }  
  t.b <- (R.hat.b-R.hat)/sd(R.hat.b.k)
}
Qt <- quantile(bootstrap, c(alpha/2, 1-alpha/2), type = 1)
R.hat+Qt*sd(bootstrap)
R.hat

## ----warning=FALSE------------------------------------------------------------
library(boot)
stat.1 <- function(data, index) { # index is the sample after resampling
  mean(data[index])
}
stat.2 <- function(data, index) { # index is the sample after resampling
  1/mean(data[index])
}
set.seed(123)
boot.1 <- boot(data = aircondit[,1], statistic = stat.1, R = 1e3)
set.seed(123)
boot.2 <- boot(data = aircondit[,1], statistic = stat.2, R = 1e3)
ci.1 <- boot.ci(boot.1, type = c("norm","basic","perc","bca"))
ci.2 <- boot.ci(boot.2, type = c("norm","basic","perc","bca"))
print(ci.1) # 1/lambda
print(ci.2) # lambda

compare <- cbind(
  rbind(ci.1$normal[2:3],ci.1$basic[4:5],ci.1$percent[4:5],ci.1$bca[4:5]),
  rbind(1/ci.2$normal[2:3],1/ci.2$basic[4:5],1/ci.2$percent[5:4],1/ci.2$bca[5:4]))
rownames(compare) <- c("norm","basic","perc","bca")
colnames(compare) <- c('ci.1.left','ci.1.right','1/ci.2.left','1/ci.2.right')
print(compare)

## ----warning=FALSE------------------------------------------------------------
library(bootstrap)
data <- as.matrix(scor)
n <- nrow(data)

stat <- function(data) {
  sigma.hat <-  var(data)
  lambda <- eigen(sigma.hat)$values
  theta.hat <- sort(lambda)[1]/sum(lambda)
  return(theta.hat)
}

theta.hat <- stat(data)
theta.jack <- numeric(n)
for (i in 1:n)
theta.jack[i] <- stat(data[-i,])
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2))
cat('bias: ',bias,'\n','se: ',se)

## ----warning=FALSE------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-two-out samples
for (index in 1:choose(n,2)) {
k <- combn(n,2)[,index]
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
yhat4 <- exp(logyhat4)
e4[k] <- magnetic[k] - yhat4
}
print(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)))

## -----------------------------------------------------------------------------
library(twosamples)
pcvm <- function(x,y){
R <- 999 #number of replicates
n <- length(x)
m <- length(y)
z <- c(x, y) #pooled sample
K <- 1:(n+m)
D <- numeric(R) #storage for replicates
D0 <- cvm_test(x, y)[1]
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
D[i] <- cvm_test(x1, y1)[1]
}
p <- mean(c(D0, D) >= D0)
return(p)
}

# Example 8.1 & 8.2
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
pcvm(x,y)
detach(chickwts)

## -----------------------------------------------------------------------------
set.seed(10)
max.num.extreme <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}

R <- 999 #number of replicates
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
z <- c(x, y) #pooled sample
K <- 1:(n1+n2)
D <- numeric(R) #storage for replicates
D0 <- max.num.extreme(x, y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n1, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
D[i] <- max.num.extreme(x1,y1)
}
p <- mean(c(D0, D) >= D0)
print(p)

## -----------------------------------------------------------------------------
set.seed(123)
LogisticIntercept <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,lambda = 1)
  x2 <- rexp(N, rate = 1)
  x3 <- sample(0:1,N,replace=TRUE)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3); p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-100,100))
  #round(unlist(solution),5)[1:3]
  alpha <- solution$root
  return(alpha)
}
N <- 1e6
b1 <- 0
b2 <- 1
b3 <- -1
f0 <- c(0.1,0.01,0.001,0.0001)

alpha <- numeric(4)
for(i in 1:4)
alpha[i] <- LogisticIntercept(N,b1,b2,b3,f0[i])
names(alpha) <- c('f0=0.1','f0=0.01','f0=0.001','f0=0.0001')
print(alpha)

plot(-log(f0),alpha)

## -----------------------------------------------------------------------------
f <- function(x){
  exp(-abs(x))/2
}

rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
  y <- rnorm(1, x[i-1], sigma)
  if (u[i] <= (f(y) / f(x[i-1])))
    x[i] <- y
  else {
    x[i] <- x[i-1]
    k <- k + 1
  }
}
return(list(x=x, k=k))
}

N <- 2000
x0 <- 10
sigma <- c(.05, .5, 2, 16)
res <- list()
# par(mfrow = c(2, 2))
acceptance.rate <- numeric(4)
for(i in 1:4){
  res[[i]] <- rw.Metropolis(sigma[i],x0,N)
  # plot(1:N,res[[i]]$x,xlab = '',ylab = 'points',type = 'l')
  acceptance.rate[i] <- (N-res[[i]]$k)/N
}
names(acceptance.rate) <- c('sigma=0.05','sigma=0.5','sigma=2','sigma=16')
print(acceptance.rate)

## -----------------------------------------------------------------------------
#initialize constants and parameters
set.seed(123)
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- 0.9 #correlation
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2
###### generate the chain #####
X[1, ] <- c(mu1, mu2) #initialize
for (i in 2:N) {
x2 <- X[i-1, 2]
m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
X[i, 1] <- rnorm(1, m1, s1)
x1 <- X[i, 1]
m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]

plot(x[,1],x[,2])
fit <- lm(x[,2]~1+x[,1])
summary(fit)

## -----------------------------------------------------------------------------
library(coda)
set.seed(123)

f <- function(x, sigma) {
if (any(x < 0)) return (0)
stopifnot(sigma > 0)
return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

Rayleigh.chain <- function(m,sigma){
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rchisq(1, df = xt)
num <- f(y, sigma) * dchisq(xt, df = y)
den <- f(xt, sigma) * dchisq(y, df = xt)
if (u[i] <= num/den) x[i] <- y 
else {
x[i] <- xt
k <- k+1 #y is rejected
}
}
return(x)
}

Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

#generate the chains
m <- 15000 #length of chains
sigma <- 4
k <- 1
mcmc <- mcmc.list()
mcmc[[k]] <- as.mcmc(Rayleigh.chain(m,sigma))
X <- Rayleigh.chain(m,sigma)
rhat <- 10
while(isTRUE(rhat[k] >= 1.2)){
  X <- rbind(X,Rayleigh.chain(m,sigma))
  psi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi))
    psi[i,] <- psi[i,] / (1:ncol(psi))
  rhat <- c(rhat,Gelman.Rubin(psi))
  k <- k+1
  mcmc[[k]] <- as.mcmc(Rayleigh.chain(m,sigma))
}
print(k)
gelman.diag(mcmc)
gelman.plot(mcmc)

## -----------------------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- u + 1
n <- length(u)

f <- function(lambda){
  sum((u*exp(-lambda*u)-v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v)))
}

# MLE
uniroot(f,c(0,1))$root

# EM
lambda.hat <- 1/mean(u)
for(k in 1:1e4){
  g <- function(lambda){
    temp <- f(lambda.hat[k])+1/lambda.hat[k]
    n/lambda-temp+f(lambda)
  }
  lambda.hat <- c(lambda.hat, uniroot(g(lambda.hat[k]),c(0,1))$root)
  if(abs(lambda.hat[k+1]-lambda.hat[k])<1e-7) break
}
lambda.hat[k]

## -----------------------------------------------------------------------------
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}

#enter the payoff matrix
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
s.A <- solve.game(A)
round(cbind(s.A$x, s.A$y), 7)
B <- A+2
s.B <- solve.game(B)
round(cbind(s.B$x, s.B$y), 7)
value.A <- t(s.A$x)%*%(A)%*%(s.A$y)
value.B <- t(s.B$x)%*%(A+2)%*%(s.B$y)
round(c(value.A,value.B),3)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
library(dplyr)
df <- data.frame(x=1:3,y=c('a','b','c'),z=c(4,6,10))
lapply(df %>% select_if(is.numeric),scale01)

## -----------------------------------------------------------------------------
df <- data.frame(x=1:3,y=c(4,6,10))
vapply(df,sd,numeric(1))
df <- data.frame(x=1:3,y=c('a','b','c'),z=c(4,6,10))
vapply(df %>% select_if(is.numeric),sd,numeric(1))

