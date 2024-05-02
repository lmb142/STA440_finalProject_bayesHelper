\name{simulationSummary}
\alias{simulationSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation Summary Statistics
}
\description{
This function takes a vector of values generated from a simulation and displays certain statistics / model diagnostics, based on whether the simulation comes from a univariate Monte Carlo, univariate Markov Chain Monte Carlo, or multivariate Markov Chain Monte Carlo simulation. Note, this function can be called directly from the Monte Carlo and Markov Chain Monte Carlo simulation functions, via the optional "summary" parameter.
}
\usage{
simulationSummary(values, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{numerical vector of simulated values
%%     ~~Describe \code{values} here~~
}
  \item{type}{string denoting where such simulated values came from. Options include "MC", "Gibbs" and "GibbsMvn"
%%     ~~Describe \code{type} here~~
}
}
\details{
The function displays expected values, variance values, Monte Carlo standard errors, distribution approximation plots, autocorrelation plots and the effective sample size. Below in the function declaration you can see for which argument "type" which statistics are outputted.
}
\value{ A series of strings with embedded integers and plots
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
simulationSummary(c(5,4,5,6,7,3,5,6,6,3,6,7,8,4,5,6,2,4), "MC")
simulationSummary(c(5,4,5,6,7,3,5,6,6,3,6,7,8,4,5,6,2,4), "Gibbs")
## The function is currently defined as
function (values, type) {
    if (type == "MC") {
        print(paste("Approximate expected value is: ", mean(values),
            sep = ""))
        print(paste("Approximate variance value is: ", var(values),
            sep = ""))
        print(paste("Monte Carlo standard error is: ", sd(values),
            sep = ""))
        print(hist(values, main = "Approxmate Distribution",
            xlab = "Value", ylab = "Frequency"))
    }
    else if (type == "Gibbs") {
        print("there")
        print(paste("Approximate expected value is: ", mean(values),
            sep = ""))
        print(paste("Approximate variance value is: ", var(values),
            sep = ""))
        print(hist(values, main = "Approxmate Distribution",
            xlab = "Value", ylab = "Frequency"))
        print(acf(values, plot = TRUE))
        print(paste("Effective sample size is: ", coda::effectiveSize(values),
            sep = ""))
    }
    else if (type == "GibbsMvn") {
        print(paste("Approximate expected value is: ", apply(values,
            2, mean), sep = ""))
        print(paste("Approximate variance value is: ", apply(values,
            2, var), sep = ""))
        print(apply(values, 2, acf))
        print(paste("Effective sample size is: ", apply(values,
            2, coda::effectiveSize), sep = ""))
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
