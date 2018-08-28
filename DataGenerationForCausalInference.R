### This generates the Simulation Data as a dataframe
# by Daniel Jacob (daniel.jacob@hu-berlin.de) 

# Arguments to specify are: 

# N = Number of Observations (real number)
# k = Number of Covariates (real number)
# random_d = treatment assignment: (Either T for random assignment or F for confounding on X)
# theta = treatment effect: (Either real number for only one theta, or "binary" {0.1,0.3} or "con" for continuous values (0.1,0.3))
# var = Size of the Variance (Noise-level)

#Required Packages
library(clusterGeneration)
library(mvtnorm)


datagen <- function(N,k,random_d,theta,var) {
  
  N=N 
  k=k 
  b=1/(1:k)
  # = Generate covariance matrix of z = #
  sigma=genPositiveDefMat(k,"unifcorrmat")$Sigma
  sigma=cov2cor(sigma)
  
  
  z=rmvnorm(N,sigma=sigma) # = Generate z = #
  
  
  ### Options for D (m_(X))
  if(random_d==T) {d <- rep(c(0,1),length.out=N)} else {d_prop <- pnorm(z%*%b) # D is dependent on Z 
  d <- as.numeric(rbinom(N,prob=d_prop, size=1))}
  
  
  ### Options for theta
  if (theta=="con") {theta_s <- as.vector(sin(z%*%b)^2) 
  theta <- (theta_s-min(theta_s))*(0.3-0.1)/(max(theta_s)-min(theta_s))+0.1} else if(theta=="binary") {theta_low <- rbinom(N,pnorm((z[,6]*(z[,1]%*%t(z[,5]))*z[,2])^2),size=1)
  theta <- ifelse(theta_low==1,0.3,0.1)} else if(theta=="big"){theta_big <- rbinom(N,pnorm((z[,6]*(z[,1]%*%t(z[,5]))*z[,2])^2),size=1)
  theta <- ifelse(theta_big==1,1,0.4)}  else {theta==theta}
  

  g=as.vector(cos(z%*%b)^2)
  
  y=theta*d+g+rnorm(N,var)
  
  data <- as.data.frame(y)
  data <- cbind(data,theta,d,z)
  colnames(data) <- c("y","theta","d",c(paste0("V", 1:k)))
  
  return(data)
} 


### Example 
test <- datagen(N=5000,k=20,random_d=F,theta="binary",var=0.25)
summary(test)

