\name{rmvnorm}
\alias{rmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Normal Distribution Function
}
\description{
A function taken directly from Peter Hoff's textbook (cited below), which takes specified model parameters and generates random samples from a Normal distribution.
}
\usage{
rmvnorm(n, mu, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{Sigma}{
%%     ~~Describe \code{Sigma} here~~
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
function (n, mu, Sigma)
{
    p <- length(mu)
    res <- matrix(0, nrow = n, ncol = p)
    if (n > 0 & p > 0) {
        E <- matrix(rnorm(n * p), n, p)
        res <- t(t(E \%*\% chol(Sigma)) + c(mu))
    }
    res
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
