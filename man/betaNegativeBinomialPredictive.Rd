\name{betaNegativeBinomialPredictive}
\alias{betaNegativeBinomialPredictive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Posterior Predictive Estimates for a Beta Prior and Negative Binomial Sampling Distribution
}
\description{
This function generates posterior predictive estimates for the data, when the prior distribution is Beta and the sampling distribution is Negative Binomial.
}
\usage{
betaNegativeBinomialPredictive(y, a, b, r, smc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data numerical vector}
  \item{a}{a numerical parameter for the beta prior distribution (alpha parameter)}
  \item{b}{a numerical parameter for the beta prior distribution (beta parameter)}
  \item{r}{stopping-time numerical parameter for a negative binomial distribution}
  \item{smc}{Monte Carlo simulation size (an integer)}
}
\details{
%%  ~~further notes~~
}
\value{
A numerical vector of the generated Monte Carlo values.
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
betaNegativeBinomialPredictive(c(0.1,0.4,0.6,0.3,0.8), 1, 5, 3, 1e5)

## The function is currently defined as
function (y, a, b, r, smc) {
  theta <- rbeta(smc, a + r*length(y), b + sum(y))
  ret <- rnbinom(smc, length(y), theta)
  return(ret)
  }
}
