[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **DataGenerationForCausalInference** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml
Name of Quantlet: DataGenerationForCausalInference

Published in: Masterthesis 'Causal Inference using Machine Learning'

Description: Generates synthetic data in form of a partial linear model to apply simulations for causal inference estimation. The parameter of interest is the treatment or uplift effect for a binary treatment assignment. 

Keywords: synthetic data, causal inference, simulation, data generation, partial linear model, treatment effect, uplift, high-dimensional

Author: Daniel Jacob


Submitted: 2018/08/24

Output: 
- Partial linear Model
- Output variable (continuous)
- Treatment paramter (different options)
- Treatment assignment (binary)
- Covariates 

```

# Data Generation for Causal Inference Simulations
Generates synthetic data to apply simulations for causal inference.

The basic model used in this function is a partially linear regression model with extensions: 

![img](http://latex.codecogs.com/svg.latex?Y%3D%5Ctheta_%7B0%7DD%2Bg_%7B0%7D%28X%29%2BU%2C%5C%5C%0D%0AD%3Dm_%7B0%7D%28X%29%2BV%2C%5C%5C%0D%0A%5Ctheta_%7B0%7D%3Dt_%7B0%7D%28Z%29%2BW%0D%0A)

The data generating process creates data of the following form.
Note that all variables are randomly generated which is why the distribution might slightly change every time.

![Data Distribution](https://github.com/QuantLet/Data_Generation/blob/master/DataGen_Distribution_Plot_different_theta.png)



### R Code
```r
### This generates the Simulation Data as a dataframe
# by Daniel Jacob (daniel.jacob@hu-berlin.de) 

# Arguments to specify are: 

# N = Number of observations (real number)
# k = Number of covariates (real number)
# random_d = treatment assignment: (Either T for random assignment or F for confounding on X)
# theta = treatment effect: (Either real number for only one theta, or "binary" {0.1,0.3} or "con" for continuous values (0.1,0.3))
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
    y1 <- theta * d + g + pnorm(rnorm(N,0,var))
    y1.1 <- rbinom(N,prob=pnorm(scale(y1)),size=1)
    #y1.1 <- (y1 - min(y1)) * (1) / (max(y1) - min(y1)) + 0
    #y <-  rbinom(N,prob=y1.1,size=1)
    y <- y1.1
  } else {y <- y1}
  
  data <- as.data.frame(y)
  data <- cbind(data, theta, d, z)
  colnames(data) <- c("y", "theta", "d", c(paste0("V", 1:k)))
  
  return(data)
}


### Example
dataset <- datagen(y="binary",N = 2000, k = 20, random_d = F, theta = "binary", var = 1)
summary(dataset)
str(dataset)


```
