## A Test file for the stabledist stuff in R
source('dist-stableMode.R')
source('utils.R')
source('dpqr-stable.R')

## we seek to test all paths of the code, there are
## two that we are concerned with (holding pm=1 as a
## constant). 
## The first path is with alpha=1 and beta=0

test  <- 1      # enable test mode

alpha <- 1
beta  <- 0
n     <- 1
gamma <- 10
delta <- 0
pm    <- 1

val1 <- rstable(n, alpha, beta, gamma, delta, pm, test)
v1 = c(alpha,beta,n,gamma,delta,pm,val1)        # the vector we will write out for cross-check

## The second path is when alpha!=1 or beta!=0
## we have a couple of tests
alpha <- 2
beta  <- 0.5
n     <- 1
gamma <- 2
delta <- 1
pm    <- 1
val2 <- rstable(n, alpha, beta, gamma, delta, pm, test)
v2 = c(alpha,beta,n,gamma,delta,pm,val2)        # the vector we will write out for cross-check

alpha <- 2
beta  <- 0.5
n     <- 1
gamma <- 4
delta <- 1
pm    <- 1
val3 <- rstable(n, alpha, beta, gamma, delta, pm, test)
v3 = c(alpha,beta,n,gamma,delta,pm,val3)        # the vector we will write out for cross-check

matout = rbind(v1,v2,v3)

write.table(matout, file="stabledist_test.txt", sep=",", append=FALSE, row.names=FALSE, col.names=FALSE)