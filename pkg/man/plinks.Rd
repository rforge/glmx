\name{plinks}
\alias{plinks}
\alias{talpha}
\alias{ao1}
\alias{ao2}
\alias{gj}
\alias{angular}
\alias{foldexp}
\alias{rocke}
\alias{nblogit}
\alias{gosset}
\alias{pregibon}
\alias{loglog}

\concept{parametric link}
\concept{transformation}

\title{Parametric Links for Binomial Generalized Linear Models}
\description{
  Various symmetric and asymmetric parametric links for use as link function for binomial generalized linear models.
}

\usage{
  gj(phi, verbose = FALSE)
  foldexp(phi, verbose = FALSE)
  ao1(phi, verbose = FALSE)
  ao2(phi, verbose = FALSE)
  talpha(alpha, verbose = FALSE, splineinv = TRUE,
    eps = 2*.Machine$double.eps, maxit = 100)
  rocke(shape1, shape2, verbose = FALSE)

  gosset(nu, verbose = FALSE)
  pregibon(a, b)

  angular(verbose = FALSE)
  loglog()
}

\arguments{
  \item{phi}{numeric parameter.}
  \item{alpha}{numeric parameter in \eqn{[0,2]}{[0,2]}.}
  \item{shape1}{numeric parameter, non-negative.}
  \item{shape2}{numeric parameter, non-negative.}
  \item{splineinv}{logical. Should a (quick and dirty) spline function be used
    for computing the inverse link function? Alternatively, a more precise but somewhat
    slower Newton algorithm is used.}
  \item{eps}{numeric. Desired convergence tolerance for Newton algorithm.}
  \item{maxit}{integer. Maximal number of steps for Newton algorithm.}
  \item{nu}{numeric parameter, non-negative.}
  \item{a}{numeric parameter.}
  \item{b}{numeric parameter.}
  \item{verbose}{logical. Should warnings about numerical issues be printed?}
}

\details{
  Symmetric and asymmetric families parametric link functions are available.
  Many families contain the logit for some value(s) of their parameter(s). 
  
  The symmetric Aranda-Ordaz (1981) transformation
  \deqn{y = \frac{2}{\phi}\frac{x^\phi-(1-x)^\phi}{x^\phi+(1-x)^\phi}}{y = \tfrac{2}{\phi}\tfrac{x^\phi-(1-x)^\phi}{x^\phi+(1-x)^\phi}}
  and the asymmetric Aranda-Ordaz (1981) transformation
  \deqn{y = \log([(1-x)^{-\phi}-1]/\phi)}{y = \log([(1-x)^{-\phi}-1]/\phi)}
  both contain the logit for \eqn{\phi = 0}{\phi = 0} and
  \eqn{\phi = 1}{\phi = 1} respectively, where the latter also includes the
  complementary log-log for \eqn{\phi = 0}{\phi = 0}.
  
  The Pregibon (1980) two parameter family is the link given by
  \deqn{y = \frac{x^{a-b}-1}{a-b}-\frac{(1-x)^{a+b}-1}{a+b}.}
  For \eqn{a = b = 0} it is the logit. For \eqn{b = 0} it is symmetric and
  \eqn{b} controls the skewness; the heavyness of the tails is controlled by
  \eqn{a}. The implementation uses the generalized lambda distribution
  \code{\link{gl}}.
  
  The Guerrero-Johnson (1982) family
  \deqn{y = \frac{1}{\phi}\left(\left[\frac{x}{1-x}\right]^\phi-1\right)}{y = \frac{1}{\phi}\left(\left[\frac{x}{1-x}\right]^\phi-1\right)}
  is symmetric and contains the logit for \eqn{\phi = 0}{\phi = 0}.
  
  The Rocke (1993) family of links is, modulo a linear transformation, the
  cumulative density function of the Beta distribution. If both parameters are
  set to \eqn{0} the logit link is obtained. If both parameters equal
  \eqn{0.5} the Rocke link is, modulo a linear transformation, identical to the
  angular transformation. Also for \code{shape1} = \code{shape2} \eqn{= 1}, the
  identity link is obtained. Note that the family can be used as a one and a two
  parameter family.
  
  The folded exponential family (Piepho, 2003) is symmetric and given by
  \deqn{y = \left\{\begin{array}{ll}
    \frac{\exp(\phi x)-\exp(\phi(1-x))}{2\phi}   &(\phi \neq 0) \\
    x- \frac{1}{2}                 &(\phi = 0)  
    \end{array}\right.}  
  
  The \eqn{t_\alpha} family (Doebler, 2011) given by
  \deqn{y = \alpha\log(x)-(2-\alpha)\log(1-x)}{y = \alpha\log(x)-(2-\alpha)\log(1-x)}
  is asymmetric and contains the logit for \eqn{\phi = 1}{\phi = 1}.
  
  The Gosset family of links is given by the inverse of the cumulative
  distribution function of the t-distribution. The degrees of freedom \eqn{\nu}
  control the heavyness of the tails and is restricted to values \eqn{>0}. For
  \eqn{\nu = 1} the Cauchy link is obtained and for \eqn{\nu \to \infty} the link
  converges to the probit. The implementation builds on \code{\link{qf}} and is
  reliable for \eqn{\nu \geq 0.2}. Liu (2004) reports that the Gosset link
  approximates the logit well for \eqn{\nu = 7}.
  
  Also the (parameterless) angular (arcsine) transformation
  \eqn{y = \arcsin(\sqrt{x})}{y = \arcsin(\sqrt{x})} is available as a link
  function.
}
  
\value{
  An object of the class \code{link-glm}, see the documentation of \code{\link{make.link}}.
}

\references{
  Pregibon D (1980). \dQuote{Goodness of Link Tests for Generalized Linear Models.}
    \emph{Journal of the Royal Statistical Society C}, \bold{29}, 15--23. 
  
  Aranda-Ordaz F (1981). \dQuote{On Two Families of Transformations to Additivity for Binary Response Data.}
    \emph{Biometrika}, \bold{68}, 357--363.

  Guerrero V, Johnson R (1982). \dQuote{Use of the Box-Cox Transformation with Binary Response Models.}
    \emph{Biometrika}, \bold{69}, 309--314.

  Rocke DM. (1993). \dQuote{On the Beta Transformation Family.}
    \emph{Technometrics}, \bold{35}, 73--81.
  
  Piepho H (2003). The Folded Exponential Transformation for Proportions.
    \emph{Journal of the Royal Statistical Society D}, \bold{52}, 575--589.
  
  Liu C (2004). \dQuote{Robit Regression: A Simple Robust Alternative to Logistic and Probit Regression.}
    In Gelman A, Meng X-L (Eds.),
    \emph{Applied Bayesian Modeling and Casual Inference from Incomplete-Data Perspectives}, Chapter 21, pp. 227--238.
    John Wiley and Sons.

  Doebler P (2011) An Asymmetric Logit Transformation. In preparation.
}

\seealso{\code{\link{family}}}

\examples{
## negative binomial and nblogit illustration
if(require("MASS") & require("pscl") & require("lmtest")) {

##########
## DATA ##
##########

## artificial data from geometric regression
set.seed(1)
n <- 200
cf <- c(0, 3)
x <- runif(n, -1, 1)
y <- rnbinom(n, mu = exp(cf[1] + cf[2] * x), size = 1)


###########################
## COUNT DATA REGRESSION ##
###########################

## negative binomial regression via MASS:::glm.nb
m_nb1 <- glm.nb(y ~ x)

## negative binomial regression via glmx
X <- model.matrix(m_nb1)
m_nb2 <- glmx.fit(X, y, family = negative.binomial, xlink = "log", xstart = 0)

## comparison
logLik(m_nb1)
logLik(m_nb2)
coeftest(m_nb1)
coeftest(m_nb2)
m_nb1$theta
exp(coef(m_nb2)[3])
m_nb1$SE.theta
exp(coef(m_nb2)[3]) * sqrt(diag(vcov(m_nb2)))[3]


############################################
## BINARY REGRESSION WITH PARAMETRIC LINK ##
############################################

## negative binomial hurdle regression via pscl:::hurdle
m_hnb1 <- hurdle(y ~ x, dist = "negbin", zero.dist = "negbin")

## negative binomial zero hurdle part via glmx
## (truncated count part will be available soon)
nbbin <- function(theta) binomial(link = nblogit(theta))
m_hnb2 <- glmx.fit(X, y > 0, family = nbbin, xlink = "log", xstart = 0)

## comparison (see only zero hurdle part)
summary(m_hnb1)
coeftest(m_hnb2)

}
}

\keyword{regression}
