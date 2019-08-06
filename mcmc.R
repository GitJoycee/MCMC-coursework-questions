######### b)
library(MASS)
data <- c(0.17,0.22,0.37,0.38,0.5,0.58,0.72,0.85,0.9,0.93,1,1,1,1,1,1,1,1,1,1) # initial data given

dist <- function(data,alpha,beta){    # Conditional distribution of alpha to be used in MH step
  (1/(gamma(alpha))^20)*(beta^((20*alpha)))*(prod(data^(alpha-1)))*alpha*exp(-alpha)
}

gibbs<- function(start.alpha,start.beta,n.sims){
  n <- length(data)
  y.mis <- integer(0)
  
  res <- matrix(NA,nrow=n.sims,ncol=2)  #Matrix of alpha,beta samples
  res[1,] <- c(start.alpha,start.beta)
  
  for (i in 2:n.sims){
    for(j in 1:10){
      while(TRUE){       #Simulates gamma random variable conditioned on it's value being greater than 1
        y.mis[j]<-rgamma(1,res[i-1,1],res[i-1,2])  
        if(y.mis[j]>1) break
      }
    }
    data <- replace(data,11:20,y.mis)      # augments the data with missing values
    x <- sum(data)
    
    res[i,2]<- rgamma(1,shape=n*res[i-1,1]+2,rate= x+1)   #Gibbs step for beta using conditional distribution of beta
    
    can <- rnorm(1,res[i-1,1],1)    #Generates a proposal from normal distribution for MH step
    acc <- min(1,dist(data,can,res[i,2])/dist(data,res[i-1],res[i,2]))   #Acceptance probability using earlier dist funciton
    
    u<-runif(1)
    if (u<acc){
      res[i,1]<-can
    }else{
      res[i,1]<-res[i-1,1]
    }
    print(i)    #Just to check progress - slow algorithm (can be improved?)
    Sys.sleep(0.01)
    flush.console()
  }
  ans <- list(sample=res,mis=y.mis)
}
samp <- gibbs(2,1,10000)$sample


plot(samp[,1],samp[,2])

mu1 <- mean(samp[,1]) # mean of alpha
mu2 <- mean(samp[,2])   # mean of beta
sigma <- var(samp)    # variance- covariance matrix of alpha and beta

######### f)

post <- function(data,alpha,beta){  # posterior we are trying to sample from
  (1/(gamma(alpha))^20)*(beta^((20*alpha)+1))*(prod(data^(alpha-1)))*alpha*exp(-beta*(sum(data)+1)-alpha)
}
con <- matrix(NA,10000,2) # Setting up matrix of alpha and beta values
con[1,] <- c(mu1,mu2) # Taking initial step to be estimates of alpha and beta

for(i in 2:10000){    #independence sampler
  for(j in 1:10){
    while(TRUE){       #Simulates gamma random variable conditioned on it's value being greater than 1
      y.mis[j] <- rgamma(1,con[i-1,1],con[i-1,2])  
      if(y.mis[j]>1) break
    }
  }
  data <- replace(data,11:20,y.mis)
  
  cand <- mvrnorm(1,c(con[i-1,1],con[i-1,2]),sigma)  #proposal from bivariate normal with mean vector and variance matrix taken from previous estimate
  accu <- min(1,post(data,cand[1],cand[2])/post(data,con[i-1,1],con[i-1,2])) # acceptance probability
  if(is.nan(accu)==TRUE){ # adjusting for possible NaN's produced when value is near 0
    accu=0
  }
  u <- runif(1)
  if(u<accu){
    con[i,] <- cand
  }else{
    con[i,] <- con[i-1,]
  }
  print(i)    #Just to check progress - slow algorithm (can be improved?)
  Sys.sleep(0.01)
  flush.console()
}

plot(con[,1],con[,2])

  
     