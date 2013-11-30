\name{pregibon}
\alias{dpregibon}
\alias{ppregibon}
\alias{qpregibon}
\alias{rpregibon}

\title{Pregibon Distribution}

\description{
  Density, distribution function, quantile function and random
  generation for the Pregibon distribution with parameters
  \code{a} and \code{b}. It is a special case of the generalized
  Tukey lambda distribution.
}

\usage{
dpregibon(x, a = 0, b = 0, log = FALSE, tol = 1e-12)
ppregibon(q, a = 0, b = 0, lower.tail = TRUE, log.p = FALSE, tol = 1e-12)
qpregibon(p, a = 0, b = 0, lower.tail = TRUE, log.p = FALSE)
rpregibon(n, a = 0, b = 0)
}

\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{a, b}{distribution parameters.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{tol}{numeric tolerance for computation of the distribution function.}
}

\value{
  \code{dpregibon} gives the probability density function,
  \code{ppregibon} gives the cumulative distribution function,
  \code{qpregibon} gives the quantile function, and
  \code{rpregibon} generates random deviates.
}

\details{
The distribution is a special case of the generalized Tukey lambda distribution
and is used by Pregibon (1980) for goodness-of-link testing. See Koenker (2006)
and Koenker and Yoon (2009) for more details.

The implementation is based on the corresponding functions for the
\link[gld]{GeneralisedLambdaDistribution} in the \pkg{gld} package (King 2013).

The corresponding link generator is available in the function \code{\link{pregibon}}.
}

\references{
  King R (2013). \dQuote{Estimation and Use of the Generalised (Tukey) Lambda Distribution.}
    R package version 2.0.1. \url{http://CRAN.R-project.org/package=gld}

  Koenker R (2006). \dQuote{Parametric Links for Binary Response.}
    \emph{R News}, \bold{6}(4), 32--34.

  Koenker R, Yoon J (2009). \dQuote{Parametric Links for Binary Choice Models: A Fisherian-Bayesian Colloquy.}
    \emph{Journal of Econometrics}, \bold{152}, 120--130.

  Pregibon D (1980). \dQuote{Goodness of Link Tests for Generalized Linear Models.}
    \emph{Journal of the Royal Statistical Society C}, \bold{29}, 15--23. 
}

\seealso{\link[gld]{GeneralisedLambdaDistribution}, \code{\link{pregibon}}}

\keyword{distribution}
