\name{plinks}
\alias{plinks}
\alias{t_alpha}
\alias{ao1}
\alias{ao2}
\alias{gj}
\alias{angular}
\alias{foldexp}
\alias{rocke}
\alias{Gosset}
\alias{Pregibon}
\concept{parametric link}
\concept{transformation}
\title{Parametric Links for Binomial Generalized Linear Models}
\description{
	Various symmetric and asymmetric parametric links for use as link function for binomial generalized linear models.
}
\usage{
	gj(phi)
	foldexp(phi)
	ao1(phi)
	ao2(phi)
	t_alpha(alpha, eps = 2*.Machine$double.eps, maxiter = 100)
	rocke(shape1, shape2 = shape1)
	Gosset(nu)
	Pregibon(a, b)
	angular()
	}
\arguments{
	\item{phi}{numeric parameter.}
	\item{alpha}{numeric parameter in \eqn{[0,2]}{[0,2]}.}
	\item{shape1}{numeric parameter, non-negative.}
	\item{shape2}{numeric parameter, non-negative.}
	\item{eps}{desired accuracy (convergence tolerance) for Newton algorithm.}
	\item{maxiter}{maximal number of steps for Newton algorithm.}
	\item{nu}{numeric parameter, non-negative.}
	\item{a}{numeric parameter.}
	\item{b}{numeric parameter.}
	}

\details{
	Symmetric and asymmetric families parametric link functions are available. Many families contain the logit for some value(s) of their parameter(s). 
	
	The symmetric Aranda-Ordaz (1981) transformation
	\deqn{y = \frac{2}{\phi}\frac{x^\phi-(1-x)^\phi}{x^\phi+(1-x)^\phi}}{y = \tfrac{2}{\phi}\tfrac{x^\phi-(1-x)^\phi}{x^\phi+(1-x)^\phi}}
	and the asymmetric Aranda-Ordaz (1981) transformation
	\deqn{y = \log([(1-x)^{-\phi}-1]/\phi)}{y = \log([(1-x)^{-\phi}-1]/\phi)}
	both contain the logit for \eqn{\phi = 0}{\phi = 0} and \eqn{\phi = 1}{\phi = 1} respectively, where the latter also includes the complementary log-log for \eqn{\phi = 0}{\phi = 0}.
	
	The Pregibon (1980) two parameter family is the link given by
	\deqn{y = \frac{x^{a-b}-1}{a-b}-\frac{(1-x)^{a+b}-1}{a+b}.}
	For \eqn{a = b = 0} it is the logit. For \eqn{b = 0} it is symmetric and \eqn{b} controls the skewness; the heavyness of the tails is controlled by \eqn{a}. The implementation uses the generalized lambda distribution \code{\link{gl}}.
	
	The Guerrero-Johnson (1982) family
	\deqn{y = \frac{1}{\phi}\left(\left[\frac{x}{1-x}\right]^\phi-1\right)}{y = \frac{1}{\phi}\left(\left[\frac{x}{1-x}\right]^\phi-1\right)}
	is symmetric and contains the logit for \eqn{\phi = 0}{\phi = 0}.
	
	The Rocke (1993) family of links is, modulo a linear transformation, the cumulative density function of the Beta distribution. If both parameters are set to \eqn{0} the logit link is obtained. If both parameters equal \eqn{0.5} the Rocke link is, modulo a linear transformation, identical to the angular transformation. Also for \code{shape1} = \code{shape2} \eqn{= 1}, the identiy link is botained. Note that the family can be used as a one and a two parameter family.
	
	The folded exponential family (Piepho, 2003) is symmetric and given by
	\deqn{y = \left\{\begin{array}{ll}
		\frac{\exp(\phi x)-\exp(\phi(1-x))}{2\phi} 	&(\phi \neq 0) \\
		x- \frac{1}{2} 								&(\phi = 0)	
		\end{array}\right.}	
	
	The \eqn{t_\alpha} family (Doebler, 2011) given by
	\deqn{y = \alpha\log(x)-(2-\alpha)\log(1-x)}{y = \alpha\log(x)-(2-\alpha)\log(1-x)} is asymmetric and contains the logit for \eqn{\phi = 1}{\phi = 1}.
	
	The Gosset family of links is given by the inverse of the cumulative distribution function of the t-distribution. The degrees of freedom \eqn{\nu} control the heavyness of the tails and is restricted to values \eqn{>0}. For \eqn{\nu = 1} the Cauchy link is obtained and for \eqn{\nu \to \infty} the link converges to the probit. The implementation builds on \code{\link{qf}} and is reliable for \eqn{\nu \geq 0.2}. Liu (2004) reports that the Gosset link approximates the logit well for \eqn{\nu = 7}.
	
	Also the (parameterless) angular (arcsine) transformation \eqn{y = \arcsin(\sqrt{x})}{y = \arcsin(\sqrt{x})} is available as a link function.
	}
	
\value{
	An object of the class \code{link-glm}, see the documentation of \code{\link{make.link}}.
	}

\author{Philipp Doebler \email{doebler@uni-muenster.de}, 
		Roger Koenker \email{rkoekner@illinois.edu}}
\references{
	Pregibon, D. (1980). Goodness of Link Tests for Generalized Linear Models. \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, \bold{29}, 15--23. 
	
	Aranda-Ordaz, F. (1981). On two families of transformations to additivity for binary response data. \emph{Biometrika}, \bold{68}, 357--363.

	Guerrero, V. & Johnson, R. (1982). Use of the Box-Cox transformation with binary response 
models. \emph{Biometrika}, \bold{69}, 309--314.

	Rocke, D.M. (1993). On the Beta Transformation Family. \emph{Technometrics}, \bold{35}, 73--81.
	
	Piepho, H. (2003). The Folded Exponential Transformation for Proportions. \emph{Journal of the Royal Statistical Society. Series D (The Statistician)}, \bold{52}, 575--589.
	
	Liu, C. (2004). Robit regression: a simple robust alternative to logistic and probit regression. In Gelman, A. and Meng, X.-L.  (Eds.), \emph{Applied Bayesian Modeling and Casual Inference from Incomplete-Data Perspectives}, Chapter 21, pp. 227-–238. Wiley.

	Doebler, P. (2011) An asymmetric logit transformation. in preparation.
	}
\seealso{
	\code{\link{family}}
	}

\examples{
	}
\keyword{regression}