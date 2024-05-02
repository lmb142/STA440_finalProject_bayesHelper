\name{univariateNormMC_varUnknown}
\alias{univariateNormMC_varUnknown}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Estimates for a Univariate Normal Model with Unknown Mean and Variance
}
\description{
This function can output either posterior estimates for the mean parameter or variance parameter. This function assumes that the prior on the sampling variance is an Inverse Gamma, the prior on the mean parameter is a Normal prior conditioned on the variance parameter, and the sampling distribution is a univariate Normal distribution.
}
\usage{
univariateNormMC_varUnknown(y, v0, sigma0, mu0, kappa0, smc, output, summary = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{y}{data numerical vector}
  \item{v0}{the shape parameter for an inverse gamma distribution (a numerical value)
%%     ~~Describe \code{v0} here~~
}
  \item{sigma0}{a component of the scale parameter for an inverse gamma distribution (a numerical value)
%%     ~~Describe \code{sigma0} here~~
}
  \item{mu0}{the mean parameter for the prior distribution on the mean of the sampling distribution (a numerical value)}
  \item{kappa0}{a component of the variance parameter for the sampling distribution (representing the number of prior measurements--minimum value 1), due to the parameterization of choice in bayesHelper (a numerical value)
%%     ~~Describe \code{kappa0} here~~
}
  \item{smc}{Monte Carlo simulation size (an integer)}
  \item{output}{a parameter to specify whether the posterior mean or posterior variance estimates will be the output (a string)
%%     ~~Describe \code{output} here~~
}
  \item{summary}{an optional parameter indicating if summary output is desired (a boolean)}

%%     ~~Describe \code{summary} here~~
}
\details{
The output variable takes a value of "mean" for the posterior mean estimate, or "variance" for the posterior variance estimate.
}
\value{
A numerical vector of the generated Monte Carlo values, or a series of summary statistics (numerical quantities embedded in strings) and a histogram based off the documentation of the simulationSummary() function.
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
univariateNormMC_varUnknown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 1e5, output="mean", summary=FALSE)
univariateNormMC_varUnknown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 1e5, output="variance", summary=FALSE)
univariateNormMC_varUnknown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 1e5, output="mean", summary=TRUE)
univariateNormMC_varUnknown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 1, 2, 4, 1, 1e5, output="variance", summary=TRUE)

## The function is currently defined as
function (y, v0, sigma0, mu0, kappa0, smc, output, summary = FALSE){
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
    return(simulationSummary(predValue, "MC"))
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
