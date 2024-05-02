\name{univariateNormMC_varKnown}
\alias{univariateNormMC_varKnown}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Estimates for a Univariate Normal Model with Unknown Mean and Known Variance
}
\description{
This function can output either posterior estimates for the mean parameter, or posterior predictive estimates for the data. This function assumes that the sampling variance is known, the prior on the mean parameter is a Normal prior conditioned on the sampling distribution's variance, and the sampling distribution is a univariate Normal distribution.
}
\usage{
univariateNormMC_varKnown(y, mu0, tau0, sigma, smc, output, summary = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data numerical vector}
  \item{mu0}{the mean parameter for the prior distribution on the mean of the sampling distribution (a numerical value)}
  \item{tau0}{the standard deviation parameter for the prior distribution on the mean parameter of the sampling distribution (a numerical value)
%%     ~~Describe \code{tau} here~~
}
  \item{sigma}{the standard deviation of the sampling distribution (a numerical value)
%%     ~~Describe \code{sigma} here~~
}
  \item{smc}{Monte Carlo simulation size (an integer)}
%%     ~~Describe \code{smc} here~~
  \item{output}{a parameter to specify whether the posterior mean or posterior predictive estimates will be the output (a string)
%%     ~~Describe \code{output} here~~
}
  \item{summary}{an optional parameter indicating if summary output is desired (a boolean)}

%%     ~~Describe \code{summary} here~~
}
\details{
The output variable takes a value of "mean" for the posterior mean estimate, or "predictive" for the posterior predictive estimate.
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
univariateNormMC_varKnown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 3, 2, 5, 1e5, output =
"mean", summary = FALSE)
univariateNormMC_varKnown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 3, 2, 5, 1e5, output =
"predictive", summary = FALSE)
univariateNormMC_varKnown(c(-1,-3,5,8,2,12,-5,7,6,9,3,2), 3, 2, 5, 1e5, output =
"mean", summary = TRUE)

## The function is currently defined as
function (y, mu0, tau0, sigma, smc, output, summary = FALSE){
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
