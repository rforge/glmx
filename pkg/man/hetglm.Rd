\name{hetglm}
\alias{hetglm}
\alias{hetglm.fit}

\alias{print.hetglm}
\alias{summary.hetglm}
\alias{print.summary.hetglm}
\alias{predict.hetglm}
\alias{coef.hetglm}
\alias{vcov.hetglm}
\alias{bread.hetglm}
\alias{estfun.hetglm}
\alias{coeftest.hetglm}
\alias{logLik.hetglm}
\alias{terms.hetglm}
\alias{model.frame.hetglm}
\alias{model.matrix.hetglm}
\alias{residuals.hetglm}
\alias{update.hetglm}

\title{Heteroskedastic Binary Response GLMs}

\description{
  Fit heteroskedastic binary response models via maximum likelihood.
}

\usage{
hetglm(formula, data, subset, na.action, weights, offset,
  family = binomial,
  link = c("probit", "logit", "cloglog", "cauchit", "log"),
  link.scale = c("log", "sqrt", "identity"),
  control = hetglm.control(\dots),
  model = TRUE, y = TRUE, x = FALSE, \dots)

hetglm.fit(x, y, z = NULL, weights = NULL, offset = NULL,
  family = binomial(), link.scale = "log", control = hetglm.control())
}

\arguments{
  \item{formula}{symbolic description of the model (of type \code{y ~ x}
    or \code{y ~ x | z}; for details see below).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of case weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor.}
  \item{family}{family object, currently only \code{\link[stats]{binomial}} is
    supported.}
  \item{link}{character specification of the link function in
    the binary response model.}
  \item{link.scale}{character specification of the link function in
    the latent scale model.}
  \item{control}{a list of control arguments specified via
    \code{\link{hetglm.control}}.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.
    For \code{hetglm.fit}, \code{x} should be a numeric regressor matrix
    and \code{y} should be the numeric response vector (with values in (0,1)).}
  \item{z}{numeric matrix. Regressor matrix for the precision model, defaulting
    to an intercept only.}
  \item{\dots}{arguments passed to \code{\link{hetglm.control}}.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"hetglm"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{terms}}, \code{\link[stats]{update}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}},
  \code{\link[sandwich]{estfun}} and \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package), and
  \code{\link[lmtest]{coeftest}} (from the \pkg{lmtest} package).  
}

\value{
  \code{hetglm} returns an object of class \code{"hetglm"}, i.e., a list with components as follows.
  \code{hetglm.fit} returns an unclassed list with components up to \code{converged}.
  \item{coefficients}{a list with elements \code{"mean"} and \code{"scale"}
    containing the coefficients from the respective models,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{output from the \code{optim} call for maximizing the log-likelihood,}
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call,}
  \item{weights}{the weights used (if any),}
  \item{offset}{the offset vector used (if any),}
  \item{n}{number of observations,}
  \item{nobs}{number of observations with non-zero weights,}
  \item{df.null}{residual degrees of freedom in the homoskedastic null model,}
  \item{df.residual}{residual degrees of freedom in the fitted model,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{loglik.null}{log-likelihood of the homoskedastic null model,}
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{family}{the family object used,}
  \item{link}{a list with elements \code{"mean"} and \code{"scale"}
    containing the link objects for the respective models,}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{a list with elements \code{"mean"}, \code{"scale"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{levels}{a list with elements \code{"mean"}, \code{"scale"} and
    \code{"full"} containing the levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"mean"} and \code{"scale"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"mean"} and \code{"scale"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}).}
}

\seealso{\code{\link[Formula]{Formula}}}

\examples{
# Heteroscedastic Probit Test 1
set.seed(48)
n <- 200
x <- exp(rnorm(n))
Y <- 5 + x + (1 + .3 * x) * rnorm(n)
y  <- factor(Y > 10) 
par(mfrow = c(1, 2))
plot(Y ~ x)
plot(y ~ x)

m0a <- glm(y ~ x, family = binomial(link = "probit"))
m0b <- hetglm(y ~ x | 1)
m <- hetglm(y ~ x)

summary(m0a)
summary(m0b)
summary(m)

if(require("lmtest")) {
lrtest(m0b, m)
}
}

\keyword{regression}
