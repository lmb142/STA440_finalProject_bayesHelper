\name{univariateNorm_Gibbs}
\alias{univariateNorm_Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Markov Chain Monte Carlo Estimates for a Univariate Normal Model with Unknown Mean and Variance
}
\description{
This function can output either posterior estimates for the mean parameter, variance parameter or posterior predictive estimates. This function assumes that the prior on the sampling variance is an Inverse Gamma, the prior on the mean parameter is a Normal prior, and the sampling distribution is a univariate Normal distribution. Note, here we sample from the full conditional posterior distributions, and the prior distributions on the mean and variance are independent.
}
\usage{
univariateNorm_Gibbs(y, v0, sigma0, mu0, tau0, smc, output, summary = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data numerical vector}
  \item{v0}{a component of the shape parameter for an inverse gamma distribution (a numerical value)}
  \item{sigma0}{a component of the scale parameter for an inverse gamma distribution (a numerical value)}
 \item{mu0}{the mean parameter for the prior distribution on the mean of the sampling distribution (a numerical value)}
 \item{tau0}{the standard deviation parameter for the prior distribution on the mean parameter of the sampling distribution (a numerical value)}
  \item{smc}{Markov Chain Monte Carlo simulation size (an integer)}
  \item{output}{a parameter to specify whether the posterior mean, posterior variance or posterior predictive estimates will be the output (a string)}
  \item{summary}{an optional parameter indicating if summary output is desired (a boolean)}

}

\details{
The output variable takes a value of "mean" for the posterior mean estimate, "variance" for the posterior variance estimate and "predictive" for the posterior predictive estimate.
}
\value{
A numerical vector of the generated Markov Chain Monte Carlo values, or a series of summary statistics (numerical quantities embedded in strings), Markov Chain Monte Carlo diagnostics, and a histogram based off the documentation of the simulationSummary() function.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lucas Bulczynski
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
univariateNorm_Gibbs(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 500, output="mean", summary=FALSE)
univariateNorm_Gibbs(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 500, output="variance", summary=FALSE)
univariateNorm_Gibbs(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 500, output="mean", summary=TRUE)
univariateNorm_Gibbs(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 500, output="predictive", summary=TRUE)

## The function is currently defined as
function (y, v0, sigma0, mu0, tau0, smc, output, summary = FALSE){
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
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
