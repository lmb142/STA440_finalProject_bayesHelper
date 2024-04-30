\name{rwish}
\alias{rwish}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wishart Distribution Function
}
\description{
A function taken directly from Peter Hoff's textbook (cited below), which takes the specified model parameters and generates random samples from a Wishart distribution.
}
\usage{
rwish(n, nu0, S0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nu0}{
%%     ~~Describe \code{nu0} here~~
}
  \item{S0}{
%%     ~~Describe \code{S0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hoff, P. D. (2009). A first course in Bayesian statistical methods. Springer. https://doi.org/10.1007/978-0-387-92407-6.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
function (n, nu0, S0)
{
    sS0 <- chol(S0)
    S <- array(dim = c(dim(S0), n))
    for (i in 1:n) {
        Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) \%*\%
            sS0
        S[, , i] <- t(Z) \%*\% Z
    }
    S[, , 1:n]
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
