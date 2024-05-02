## metadata
rm(list = ls())
install.packages("devtools")
install.packages("coda")
library(devtools)
library(coda)


## Standalone function
## Monte Carlo sampling event handler function
MC <- function(y, a, b, r = NULL, m = NULL, smc, prior, sampling, predictive = FALSE, summary = FALSE) {
  if (prior=="gamma" & sampling=="poisson" & predictive == FALSE & summary == FALSE){
    return (gammaPoisPosterior(y, a, b, smc))
  } else if (prior=="gamma" & sampling=="poisson" & predictive == FALSE & summary == TRUE){
    return(simulationSummary(gammaPoisPosterior(y, a, b, smc), type="MC"))
  } else if (prior=="gamma" & sampling=="poisson" & predictive == TRUE & summary == FALSE){
    return (gammaPoisPredictive(y, a, b, smc))
  } else if (prior=="gamma" & sampling=="poisson" & predictive == TRUE & summary == TRUE){
    return (simulationSummary(gammaPoisPredictive(y, a, b, smc), type="MC"))

  } else if (prior=="beta" & sampling=="bernoulli" & predictive == FALSE & summary == FALSE){
    return (betaBernoulliPosterior(y, a, b, smc))
  } else if (prior=="beta" & sampling=="bernoulli" & predictive == FALSE & summary == TRUE){
    return (simulationSummary(betaBernoulliPosterior(y, a, b, smc), type="MC"))
  } else if (prior=="beta" & sampling=="bernoulli" & predictive == TRUE & summary == FALSE){
    return (betaBernoulliPredictive(y, a, b, smc))
  } else if (prior=="beta" & sampling=="bernoulli" & predictive == TRUE & summary == TRUE){
    return (simulationSummary(betaBernoulliPredictive(y, a, b, smc), type="MC"))

  } else if (prior=="beta" & sampling=="binomial" & predictive == FALSE & summary == FALSE){
    return (betaBinomialPosterior(y, a, b, m, smc))
  } else if (prior=="beta" & sampling=="binomial" & predictive == FALSE & summary == TRUE){
    return (simulationSummary(betaBinomialPosterior(y, a, b, m, smc), type="MC"))
  } else if (prior=="beta" & sampling=="binomial" & predictive == TRUE & summary == FALSE){
    return (betaBinomialPredictive(y, a, b, m, smc))
  } else if (prior=="beta" & sampling=="binomial" & predictive == TRUE & summary == TRUE){
    return (simulationSummary(betaBinomialPredictive(y, a, b, m, smc), type="MC"))


  } else if (prior=="beta" & sampling=="negativebinomial" & predictive == FALSE & summary == FALSE){
    return (betaNegativeBinomialPosterior(y, a, b, r, smc))
  } else if (prior=="beta" & sampling=="negativebinomial" & predictive == FALSE & summary == TRUE){
    return (simulationSummary(betaNegativeBinomialPosterior(y, a, b, r, smc), type="MC"))
  } else if (prior=="beta" & sampling=="negativebinomial" & predictive == TRUE & summary == FALSE){
    return (betaNegativeBinomialPredictive(y, a, b, r, smc))
  } else if (prior=="beta" & sampling=="negativebinomial" & predictive == TRUE & summary == TRUE){
    return (simulationSummary(betaNegativeBinomialPredictive(y, a, b, r, smc), type="MC"))

  } else if (prior=="beta" & sampling=="geometric"  & predictive == FALSE & summary == FALSE){
    return(betaGeometricPosterior(y, a, b, smc))
  } else if (prior=="beta" & sampling=="geometric"  & predictive == FALSE & summary == TRUE){
    return(simulationSummary(betaGeometricPosterior(y, a, b, smc), type="MC"))

  } else if (prior=="gamma" & sampling=="exponential" & predictive==FALSE & summary == FALSE){
    return(gammaExponentialPosterior(y, a, b, smc))
  } else if (prior=="gamma" & sampling=="exponential" & predictive==FALSE & summary == TRUE){
    return(simulationSummary(gammaExponentialPosterior(y, a, b, smc), type="MC"))
  } else{
    return("ERROR: chosen distribution combinations are not supported!")
  }
}


simulationSummary <- function(values, type){
  if (type=="MC"){
    print(paste("Approximate expected value is: ", mean(values), sep=""))
    print(paste("Approximate variance value is: ", var(values), sep=""))
    print(paste("Monte Carlo standard error is: ", sd(values), sep=""))
    print(hist(values, main = "Approxmate Distribution", xlab = "Value", ylab = "Frequency"))
  } else if (type=="Gibbs") {
    print("there")
    print(paste("Approximate expected value is: ", mean(values), sep=""))
    print(paste("Approximate variance value is: ", var(values), sep=""))
    print(hist(values, main = "Approxmate Distribution", xlab = "Value", ylab = "Frequency"))
    print(acf(values, plot = TRUE))
    print(paste("Effective sample size is: ", coda::effectiveSize(values), sep=""))
  } else if (type=="GibbsMvn"){
    print(paste("Approximate expected value is: ", apply(values, 2, mean), sep=""))
    print(paste("Approximate variance value is: ", apply(values, 2, var), sep=""))
    print(apply(values, 2, acf))
    print(paste("Effective sample size is: ", apply(values, 2, coda::effectiveSize), sep=""))
  }
}



## Normal prior on mean & Normal sampling distribution with variance known, mean posterior outputs
univariateNormMC_varKnown <- function(y, mu0, tau0, sigma, smc, output, summary = FALSE){
  varianceValue = 1 / (1/tau0^2 + length(y) / sigma^2)
  meanValue = (mu0 /tau0^2 + sum(y) / sigma^2) * varianceValue
  thetaValue <- rnorm(smc, meanValue, sqrt(varianceValue))
  predValue <- rnorm(smc, meanValue, sqrt(varianceValue + sigma^2))
  if (output=="mean" & summary == FALSE){
    return(thetaValue)
  } else if (output=="predictive" & summary == FALSE){
    return(predValue)
  } else if (output=="mean" & summary == TRUE){
    return(simulationSummary(thetaValue, "MC"))
  } else if (output=="predictive" & summary == TRUE){
    return(simulationSummary(predValue, "MC"))
  }
}


## Standalone function
## Normal prior on mean, Normal sampling distribution & inverse gamma prior on variance
## posterior outputs for mean & variance
univariateNormMC_varUnknown <- function(y, v0, sigma0, mu0, kappa0, smc, output, summary = FALSE){
  n <- length(y)
  vn <- v0 + n
  sigman_squared <- 0.5 * ((length(y) - 1)*var(y) + v0*sigma0^2 + (n*kappa0)/(n + kappa0)*(mean(y)-mu0)^2)
  sigmaValue <- 1 / rgamma(smc, vn, sigman_squared)
  thetaValue <- rnorm(smc, (kappa0*mu0 + n*mean(y))/(kappa0 + n), sqrt(sigmaValue/(kappa0 + n)))
  if (output=="mean" & summary == FALSE){
    return(thetaValue)
  } else if (output=="variance" & summary == FALSE){
    return(sigmaValue)
  } else if (output=="mean" & summary == TRUE){
    return(simulationSummary(thetaValue, "MC"))
  } else if (output=="variance" & summary == TRUE){
    return(simulationSummary(sigmaValue, "MC"))
  }
}

## Standalone Function
## Gibbs sampler for Mean and Variance of Univariate Normal Distribution
univariateNorm_Gibbs <- function(y, v0, sigma0, mu0, tau0, smc, output, summary = FALSE){
  vn <- v0 + length(y)
  theta_0 <- mean(y)
  var_0 <-  var(y)
  PHI <- matrix(nrow = smc, ncol=2)
  PHI[1,] <- phi <- c(theta_0, var_0)
  PRED <- rnorm(1, theta_0, sqrt(var_0))
  for (i in 2:smc){
    sigma_n <- (v0*sigma0^2 + (length(y) - 1) * var_0 + length(y) * (theta_0-phi[1])^2)
    phi[2] <- 1 / rgamma(1, vn/2, sigma_n/2)
    tau_n <- 1 / (1/tau0^2 + length(y) / phi[2])
    mu_n <- (mu0/tau0^2 + sum(y)/phi[2]) * tau_n
    phi[1] <- rnorm(1, mu_n, sqrt(tau_n))
    PHI[i,] <- phi
    PRED <- rbind(PRED, rnorm(1, phi[1], sqrt(phi[2])))
  }
  if (output=="mean" & summary ==FALSE){
    return (PHI[,1])
  } else if (output=="variance" & summary==FALSE){
    return (PHI[,2])
  } else if (output=="predictive" & summary==FALSE){
    return (PRED)
  } else if (output=="mean" & summary ==TRUE){
    return (simulationSummary(PHI[,1], "Gibbs"))
  } else if (output=="variance" & summary==TRUE){
    return (simulationSummary(PHI[,2], "Gibbs"))
  } else if (output=="predictive" & summary==TRUE){
    return (simulationSummary(PRED, "Gibbs"))
  }
}




## Standalone Function
## Note: this function uses logic from Peter Hoff's textbook: "A First Course in Bayesian
## Statistical Methods"
mvn_Gibbs <- function(y, mu0, Lambda0, v0, S0, smc, output, summary = FALSE) {
  n <- dim(y)[1]
  ybar <- apply(y, 2, mean)
  Sigma <- cov(y)
  THETA <- SIGMA <- PRED <- NULL
  for (i in 1:smc){
    Lambdan <- solve(solve(Lambda0) + n * solve(Sigma))
    mun <- (solve(Lambda0)%*%mu0 + n*solve(Sigma)%*%ybar)
    theta <- rmvnorm(1, Lambdan%*%mun, Lambdan)
    S_theta <- (t(y)-c(theta))%*%t(t(y)-c(theta))
    Sigma <- solve(rwish(1, v0+n, solve(S0 + S_theta)))
    pred <- rmvnorm(1, theta, Sigma)
    PRED <- rbind(PRED, pred)
    THETA <- rbind(THETA, theta)
    SIGMA <- rbind(SIGMA, c(Sigma))
  }
  if (output=="mean" & summary ==FALSE){
    return (THETA)
  } else if (output=="covariance" & summary==FALSE){
    return (SIGMA)
  } else if (output=="predictive" & summary==FALSE){
    return (PRED)
  } else if (output=="mean" & summary ==TRUE){
    return (simulationSummary(THETA, "GibbsMvn"))
  } else if (output=="covariance" & summary==TRUE){
    return ("ERROR: summary statistics not supported for covariance matrix")
  } else if (output=="predictive" & summary==TRUE){
    return (simulationSummary(PRED, "GibbsMvn"))
  }
}


## The below functions are helper functions ------------------------------------

## Gamma prior & Pois sampling distribution with posterior outputs
gammaPoisPosterior <- function(y, a, b, smc){
  ret <- rgamma(smc, a + sum(y), b + length(y))
  return (ret)
}

## Gamma prior & Pois sampling distribution with predictive outputs
gammaPoisPredictive <- function(y, a, b, smc){
  theta <- rgamma(smc, a + sum(y), b + length(y))
  ret <- rpois(smc, theta)
  return (ret)
}

## Beta prior & Bernoulli sampling distribution with posterior outputs
betaBernoulliPosterior <- function(y, a, b, smc) {
  ret <- rbeta(smc, a + sum(y), b + length(y) - sum(y))
  return (ret)
}

## Beta prior & Bernoulli sampling distribution with predictive outputs
betaBernoulliPredictive <- function(y, a, b, smc) {
  theta <- rbeta(smc, a + sum(y), b + length(y) - sum(y))
  ret <- rbinom(smc, 1, theta)
  return (ret)
}

## Beta prior & Binomial sampling distribution with posterior outputs
betaBinomialPosterior <- function(y, a, b, m, smc) {
  ret <- rbeta(smc, a + sum(y), b + m*length(y) - sum(y))
  return (ret)
}

## Beta prior & Binomial sampling distribution with predictive outputs
betaBinomialPredictive <- function(y, a, b, m, smc) {
  theta <- rbeta(smc, a + sum(y), b + m*length(y) - sum(y))
  ret <- rbinom(smc, m, theta)
  return (ret)
}


## Beta prior & Negative Binomial sampling distribution with posterior outputs
betaNegativeBinomialPosterior <- function(y, a, b, r, smc){
  ret <- rbeta(smc, a + r*length(y), b + sum(y))
  return (ret)
}

## Beta prior & Negative Binomial sampling distribution with predictive outputs
betaNegativeBinomialPredictive <- function(y, a, b, r, smc){
  theta <- rbeta(smc, a + r*length(y), b + sum(y))
  ret <- rnbinom(smc, length(y), theta)
  return (ret)
}


## Beta prior & Geometric sampling distribution with posterior outputs
betaGeometricPosterior <- function(y, a, b, smc){
  ret <- rbeta(smc, a + length(y), b + sum(y))
  return (ret)
}

## Gamma prior & Exponential sampling distribution with posterior outputs
gammaExponentialPosterior <- function(y, a, b, smc){
  ret <- rgamma(smc, a + length(y), b + sum(y))
  return (ret)
}

## Note: this helper function is taken from Peter Hoff's textbook: "A First Course in Bayesian
## Statistical Methods"
rwish<-function(n,nu0,S0){
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n){
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

## Note: this helper function is taken from Peter Hoff's textbook: "A First Course in Bayesian
## Statistical Methods"
rmvnorm<-
  function(n,mu,Sigma) {
    p<-length(mu)
    res<-matrix(0,nrow=n,ncol=p)
    if( n>0 & p>0 ) {
      E<-matrix(rnorm(n*p),n,p)
      res<-t(  t(E%*%chol(Sigma)) +c(mu))
    }
    res
  }

