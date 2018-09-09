# Time Varying Quadratic Discriminant Analysis

# EEG Data Application

# install.packages("glasso", dependencies = T)
Sys.getenv("R_LIBS_USER")
library(glasso)
library(MASS)


# Kernel Function, using Gaussian Kernel
GetKernel.Gaussian <- function(x, y, h = 0.2){
  t <- abs(x - y) / h
  if(abs(t) <= 1){Kernelvalue = 0.5}
  else {Kernelvalue <- 0}
  return(Kernelvalue)
}

# Epanechnikov kernel
GetKernel.Epan <- function(x, y, h){
  t <- abs(x - y) / h
  if(abs(t) <= 1){
    Kernelvalue <- 0.75 * (1 - t^2)
  }
  else {Kernelvalue <- 0}

  return(Kernelvalue)
}


# Tri-Cube kernel
GetKernel.TriCu <- function(x, y, h){
  t <- abs(x - y) / h
  if(abs(t) <= 1){
    Kernelvalue <- (1 - abs(t)^3)^3
  }
  else {Kernelvalue <- 0}

  return(Kernelvalue)
}


# Function for estimating functional graphical model
# X.Data should contain trail, timepoint and new time
# Need to load more data, otherwise there may not observations at same time point and resulting in
# no covariance matrix
Est.Kernel <- function(X.Data, Time.new, h, flag.kernel){
  n.trail <- length(unique(X.Data$Trail))
  Time.seq <- sort(unique(X.Data$Time))
  
  KernelSum <- 0
  KernelSum.Vec <- rep(0, 64)
  KernelSum.Cov <- matrix(0, 64, 64)
  for (j in 1 : length(Time.seq)){
    SubDataset <- X.Data[which(X.Data$Time == Time.seq[j]), c(3 : 66)]
    X.bar <- colMeans(SubDataset)
    X.Cov <- cov(SubDataset)
    if(flag.kernel == 1){KernelValue = GetKernel.Gaussian(Time.seq[j], Time.new, h)}
    else if(flag.kernel == 2) {KernelValue = GetKernel.Epan(Time.seq[j], Time.new, h)}
    else {KernelValue = GetKernel.TriCu(Time.seq[j], Time.new, h)}
    KernelSum <- KernelSum + KernelValue 
    KernelSum.Vec <- KernelSum.Vec + (KernelValue * X.bar)
    KernelSum.Cov <- KernelSum.Cov + (KernelValue * X.Cov)
  }
  
  FGG.Mean <- KernelSum.Vec / KernelSum
  FGG.Cov <- KernelSum.Cov / KernelSum
  return(list(Est.Mean = FGG.Mean, Est.Cov = FGG.Cov))
}


# Estimating the precision Matrix with Glasso, use function 'glasso'


# This is the main function of the procedure
# Input: New observation, tunning parameter, at specific time point 
# Output: Posterior probability of classification, New observation is class 1
Time.QDA <- function(Testing, Training.X, Training.Y, h = 0.2, Time.seq.new, priorX = 0.5, flag.kernel){
  
  priorY <- 1 - priorX
  d <- length(Time.seq.new)
  Vote <- rep(0, d)
  
  Test.set <- Testing[, c(3 : 66)]
  
  for (i in 1 : length(Time.seq.new)){
    
    Est.class1 <- Est.Kernel(Training.X, Time.new = Time.seq.new[i], h = 0.2, flag.kernel = flag.kernel)
    Est.class2 <- Est.Kernel(Training.Y, Time.new = Time.seq.new[i], h = 0.2, flag.kernel = flag.kernel)
  
  
  # Glasso Procedure
  Prec.Matrix1 <- glasso(Est.class1$Est.Cov, rho = 0.5)$wi
  Prec.Matrix2 <- glasso(Est.class2$Est.Cov, rho = 0.5)$wi
  
  Vec1.temp <- as.vector((Test.set[i, ] - Est.class1$Est.Mean), mode = 'numeric')
  Vec2.temp <- as.vector((Test.set[i, ] - Est.class2$Est.Mean), mode = 'numeric')
  dis1 <- 0.5 * log(det(Prec.Matrix1)) - 0.5 * t(Vec1.temp) %*% Prec.Matrix1 %*% Vec1.temp + log(priorX)
  dis2 <- 0.5 * log(det(Prec.Matrix2)) - 0.5 * t(Vec2.temp) %*% Prec.Matrix2 %*% Vec2.temp + log(priorY)
    
    if (dis1 > dis2){Vote[i] = 1}
  }
  Post.Prob <- mean(Vote)
  return(Post.Prob)
}



