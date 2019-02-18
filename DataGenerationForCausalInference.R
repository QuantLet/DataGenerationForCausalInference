### This generates the Simulation Data as a dataframe
# by Daniel Jacob (daniel.jacob@hu-berlin.de) 

# Arguments to specify are: 
# Y = Outcome (dependend variable). Either continuous or binary. 
# N = Number of observations (real number)
# k = Number of covariates (real number). At least 10 
# random_d = treatment assignment: (Either T for random assignment or F for confounding on X)
# theta = treatment effect: (Either real number for only one theta, or "binary" {0.1,0.3}, "con" for continuous values (0.1,0.3) or "big" for {1,0.4})
# var = Size of the variance (Noise-level)

#Required Packages
if(!require("clusterGeneration")) install.packages("clusterGeneration"); library("clusterGeneration")
if(!require("mvtnorm")) install.packages("mvtnorm"); library("mvtnorm")



datagen <- function(N,y,k,random_d,theta,var) {
  
  N = N
  k = k
  b = 1 / (1:k)
  # = Generate covariance matrix of z = #
  sigma <- genPositiveDefMat(k, "unifcorrmat")$Sigma
  sigma <- cov2cor(sigma)
  
  
  z <- rmvnorm(N, sigma = sigma) # = Generate z = #
  
  
  ### Options for D (m_(X))
  if (random_d == T) {
    d <- rep(c(0, 1), length.out = N)
  } else {
    d_prop <- pnorm(z %*% b) # D is dependent on Za
    d <- as.numeric(rbinom(N, prob = d_prop, size = 1))
  }
  
  
  ### Options for theta
  if (theta == "con") {
    theta_s <- as.vector(sin(z %*% b) ^ 2)
    theta <-
      (theta_s - min(theta_s)) * (0.3 - 0.1) / (max(theta_s) - min(theta_s)) +
      0.1
  } else if (theta == "binary") {
    theta_low <-
      rbinom(N, pnorm((z[, 6] * (z[, 1] %*% t(
        z[, 5]
      )) * z[, 2]) ^ 2), size = 1)
    theta <-
      ifelse(theta_low == 1, 0.3, 0.1)
  } else if (theta == "big") {
    theta_big <-
      rbinom(N, pnorm((z[, 6] * (z[, 1] %*% t(
        z[, 5]
      )) * z[, 2]) ^ 2), size = 1)
    theta <- ifelse(theta_big == 1, 1, 0.4)
  }  else {
    theta == theta
  }
  
  
  g <- as.vector(cos(z %*% b) ^ 2)
  
  if(y=="binary") {
    y1 <- theta * d + g 
    y1.1 <- rbinom(N,prob=pnorm(scale(y1)),size=1)
    y <- y1.1
  } else {y <- theta * d + g + rnorm(N,0,var)}
  
  data <- as.data.frame(y)
  data <- cbind(data, theta, d, z)
  colnames(data) <- c("y", "theta", "d", c(paste0("V", 1:k)))
  
  return(data)
}


### Example
dataset <- datagen(y="binary",N = 2000, k = 20, random_d = F, theta = "binary", var = 1)
summary(dataset)
str(dataset)
