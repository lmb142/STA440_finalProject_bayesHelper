\name{uniformBernoulliPosterior}
\alias{uniformBernoulliPosterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo Posterior Estimates for a Uniform Prior and Bernoulli Sampling Distribution
}
\description{
This function generates posterior estimates for the probability of success, when the prior distribution is Uniform and the sampling distribution is Bernoulli.
}
\usage{
uniformBernoulliPosterior(y, smc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data numerical vector}
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
uniformBernoulliPosterior(c(0.1,0.4,0.6,0.3,0.8), 1e5)

## The function is currently defined as
function (y, smc)
{
    ret <- rbeta(smc, sum(y), length(y) - sum(y))
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
