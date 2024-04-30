\name{mvn_Gibbs}
\alias{mvn_Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Markov Chain Monte Carlo Estimates for a Multivariate Normal Model with Unknown Mean and Variance
}
\description{
This function can output either posterior estimates for the mean vector parameter (1 x p) or covariance matrix (p x p), or posterior predictive estimates for the data vector (1 x p). This function assumes that the prior on the sampling covariance matrix is a Wishart distribution, the prior on the mean vector is a multivariate Normal prior, and the sampling distribution is a multivariate Normal distribution. Note, here we sample from the full conditional posterior distributions, and the prior distributions on the mean and variance are independent.
}
\usage{
mvn_Gibbs(y, mu0, Lambda0, v0, S0, smc, output = NULL, summary = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the data, represented as a n x p numerical matrix (where n refers to the number of observations)}
  \item{mu0}{the mean vector parameter for the prior distribution on the mean of the sampling distribution (a numerical 1 x p vector)}
  \item{Lambda0}{the covariance matrix parameter for the prior distribution on the mean parameter of the sampling distribution (a numerical p x p matrix)}
  \item{v0}{the degrees of freedom parameter for the Wishart distribution (a numerical value)
%%     ~~Describe \code{v0} here~~
}
  \item{S0}{the scale matrix of the Wishart Distribution (a numerical p x p matrix)
%%     ~~Describe \code{S0} here~~
}
  \item{smc}{Markov Chain Monte Carlo simulation size (an integer)}
  \item{output}{an optional parameter to specify whether the posterior mean vector, posterior covariance matrix, or predictive posterior vector estimates will be the output. By default both are printed in a list, with the posterior mean vector estimates being printed in first (a string)}
  \item{summary}{an optional parameter indicating if summary output is desired (a boolean)}
}
\details{
The output variable takes a value of "mean" for the posterior mean vector estimate, "covariance" for the posterior covariance matrix estimate, or "predictive" for the posterior predictive vector estimates.
}
\value{
A numerical vector of the generated Markov Chain Monte Carlo vector estimates, a numerical matrix, or a series of summary statistics (numerical quantities embedded in strings) and Markov Chain Monte Carlo diagnostics.
}
\references{
Hoff, P. D. (2009). A first course in Bayesian statistical methods. Springer. https://doi.org/10.1007/978-0-387-92407-6.
}
\author{
Lucas Bulczynski
}
\note{
This function uses logic from Peter Hoff's textbook cited below. Furthermore, the example uses publicly available data from the textbook, as devising a set of parameters that fit matrix invertibility requirements proved difficult.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
mvn_Gibbs(matrix(c(59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49, 72, 34,
70, 34, 50, 41, 52, 60, 34, 28, 35, 77, 39, 46, 26, 38, 43, 68,
86, 77, 60, 50, 59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33), ncol=2, nrow=22), c(50, 50), matrix(c(625, 312.5, 312,5, 625), nrow=2, ncol=2), 4, matrix(c(625, 312.5, 312,5, 625), nrow=2, ncol=2)
, 500, output= "mean", summary=FALSE)

mvn_Gibbs(matrix(c(59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49, 72, 34,
70, 34, 50, 41, 52, 60, 34, 28, 35, 77, 39, 46, 26, 38, 43, 68,
86, 77, 60, 50, 59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33), ncol=2, nrow=22), c(50, 50), matrix(c(625, 312.5, 312,5, 625), nrow=2, ncol=2), 4, matrix(c(625, 312.5, 312,5, 625), nrow=2, ncol=2)
, 500, output= "covariance", summary=FALSE)

mvn_Gibbs(matrix(c(59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49, 72, 34,
70, 34, 50, 41, 52, 60, 34, 28, 35, 77, 39, 46, 26, 38, 43, 68,
86, 77, 60, 50, 59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33), ncol=2, nrow=22), c(50, 50), matrix(c(625, 312.5, 312,5, 625), nrow=2, ncol=2), 4, matrix(c(625, 312.5, 312,5, 625), nrow=2, ncol=2)
, 500, output= "mean", summary=TRUE)


## The function is currently defined as
function (y, mu0, Lambda0, v0, S0, smc, output=NULL, summary = FALSE) {
  n <- dim(y)[1]
  ybar <- apply(y, 2, mean)
  Sigma <- cov(y)
  THETA <- SIGMA <- PRED <- NULL
  for (i in 1:smc){
    Lambdan <- solve(solve(Lambda0) + n * solve(Sigma))
    mun <- (solve(Lambda0)\%*\%mu0 + n*solve(Sigma)\%*\%ybar)
    theta <- rmvnorm(1, Lambdan\%*\%mun, Lambdan)
    S_theta <- (t(y)-c(theta))\%*\%t(t(y)-c(theta))
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
  } else{
    return (THETA)
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
