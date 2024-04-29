\name{betaGeometricPosterior}
\alias{betaGeometricPosterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Posterior Estimates for a Beta Prior and Geometric Sampling Distribution
}
\description{
This function generates posterior estimates for the probability of success, when the prior distribution is Beta and the sampling distribution is Geometric.
}
\usage{
betaGeometricPosterior(y, a, b, smc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data numerical vector}
  \item{a}{a numerical parameter for the beta prior distribution (alpha parameter)}
  \item{b}{a numerical parameter for the beta prior distribution (beta parameter)}
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
betaGeometricPosterior(c(0.1,0.4,0.6,0.7,0.8), 1, 2, 1e5)

## The function is currently defined as
function (y, a, b, smc)
{
    ret <- rbeta(smc, a + length(y), b + sum(y))
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
