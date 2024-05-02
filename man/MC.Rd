\name{MC}
\alias{MC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Simulation Generator
}
\description{
This function generates a Monte Carlo simulation for the posterior or posterior predictive distribution with the specified prior and sampling distribution, parameters and sample size.
}
\usage{
MC(y, a, b, r = NULL, m = NULL, smc, prior, sampling, predictive = FALSE, summary = FALSE)
}

\arguments{
  \item{y}{data numerical vector}
  \item{a}{a numerical parameter for certain prior distributions (alpha)}
  \item{b}{a numerical parameter for certain prior distributions (beta)}
  \item{r}{stopping-time numerical parameter for a negative binomial distribution}
  \item{m}{the number of trials}
  \item{smc}{Monte Carlo simulation size (an integer)}
  \item{prior}{prior distribution (a string)}
  \item{sampling}{sampling distribution (a string)}
  \item{predictive}{an optional parameter indicating if the posterior predictive output is desired (a boolean). By default, the value is set to false}
  \item{summary}{an optional parameter indicating if summary output is desired for a given simulation (a boolean). By default, the value is set to false}
}
\details{
This function calls a series of helper functions, based on which prior and sampling distributions are passed to the function, among other parameters. See in the function definition below what these helper functions are named, and access the associated help files to learn more about each one.

The supported prior and sampling distribution pairs for posterior output are: Gamma prior and Poisson sampling, Beta prior and Bernoulli sampling, Beta prior and Binomial sampling, Beta prior and Negative Binomial sampling, Beta prior and Geometric sampling and Gamma prior and Exponential sampling.

The supported prior and sampling distribution pairs for posterior predictive output are: Gamma prior and Poisson sampling, Beta prior and Bernoulli sampling, Beta prior and Binomial sampling and Beta prior and Negative Binomial sampling.
}
\value{
A numerical vector of the generated Monte Carlo values, or a series of summary statistics (numerical quantities embedded in strings) and a histogram generated in the simulationSummary() helper function.
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
MC(c(1,5,3,8,5,9,3,8,4), a = 1, b = 2, smc = 1e5, prior = "gamma", sampling = "poisson")
MC(c(0.1,0.3,0.5, 0.7, 0.4, 0.9), a = 1, b = 2, smc = 1e5, prior = "beta", sampling = "bernoulli", summary=TRUE, predictive=FALSE)
MC(c(0.1,0.3,0.5, 0.7, 0.4, 0.9), a = 2, b = 6, r = 3, smc = 1e5, prior = "beta", sampling = "negativebinomial", summary=FALSE, predictive=FALSE)
MC(c(0.1,0.3,0.5, 0.7, 0.4, 0.9), a = 2, b = 6, smc = 1e5, prior = "gamma", sampling = "exponential", summary=TRUE, predictive=FALSE)

## The function is currently defined as
function (y, a, b, r = NULL, m = NULL, smc, prior, sampling, predictive = FALSE, summary = FALSE) {
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
