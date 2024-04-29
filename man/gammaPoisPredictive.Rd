\name{gammaPoisPredictive}
\alias{gammaPoisPredictive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Posterior Predictive Estimates for a Gamma Prior and Poisson Sampling Distribution
}
\description{
This function generates posterior predictive estimates for the data, when the prior distribution is Gamma and the sampling distribution is Poisson.
}
\usage{
gammaPoisPredictive(y, a, b, smc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data numerical vector}
  \item{a}{a numerical parameter for the gamma prior distribution (alpha parameter)}
  \item{b}{a numerical parameter for the gamma prior distribution (beta parameter)}
  \item{smc}{Monte Carlo simulation (an integer)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
gammaPoisPredictive(c(1,3,5,2,6,3,3), 1, 3, 1e5)

## The function is currently defined as
function (y, a, b, smc)
{
    theta <- rgamma(smc, a + sum(y), b + length(y))
    ret <- rpois(smc, theta)
    return(ret)
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
