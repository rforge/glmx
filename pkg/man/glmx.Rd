\name{glmx}
\alias{glmx}
\alias{glmx.fit}

\title{Parametric Links for Binomial Generalized Linear Models}

\description{
  Various symmetric and asymmetric parametric links for use as
  link function for binomial generalized linear models.
}

\usage{
glmx.fit(x, y, weights = NULL, offset = NULL, start = NULL,
  family = negative.binomial, xlink = "log", xstart = 0,
  profile = TRUE, nuisance = FALSE, control = glm.control(), xcontrol = list())
}

\arguments{
  \item{weights}{optional numeric vector of case weights.}
  \item{offset}{optional numeric vector(s) with an a priori known component to be
    included in the linear predictor(s).}
  \item{family}{family object (including the link function of the mean model).}
  \item{control}{a list of control arguments.}
  \item{y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.
    For \code{glmx.fit}, \code{x} should be a numeric regressor matrix
    and \code{y} should be the numeric response vector (with values in (0,1)).}
  \item{start, xlink, xstart, profile, nuisance, xcontrol}{further arguments}
}

\details{
  to do
}
  
\value{
  An object of the class \code{glmx}.
}

\seealso{\code{\link{hetglm}}}

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
