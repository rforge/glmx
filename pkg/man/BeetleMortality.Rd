\name{BeetleMortality}
\docType{data}
\alias{bliss}
\alias{beetle}
\alias{BeetleMortality}
\title{Bliss (1935) Beetle Mortality Data}
\description{Mortality of adult flour beetle after five hours' exposure to gaseous carbon disulphide.}

\usage{data("BeetleMortality")}

\format{A data frame with 8 rows and 3 columns. 
	\describe{
		\item{l10d}{numeric. \eqn{\log_{10}} dose}		\item{died}{integer. number killed}
		\item{n}{integer. number exposed}
		}
	}
	
\details{The data originates from Bliss (1935) and has been re-analysed frequently.}	

\source{
	Bliss, C. I. (1935)  The calculation of the dosage-mortality curve. \emph{Ann. Appl. Biology}, \bold{22}, 134--167.
}
\references{
Prentice, R. L (1976) A generalization of the probit and logit methods for dose response curves. \emph{Biometrics}, \bold{38}, 761--768. 

Pregibon, D. (1980) Goodness of link tests for generalized linear models. \emph{Appl. Statist.}, \bold{29}, 15--24. 
	
	Aranda-Ordaz, F. (1981) On two families of transformations to additivity for binary response data. \emph{Biometrika}, \bold{68}, 357--363.

Hauck, W. (1990) Choice of Scale and Asymmetric Logistic Models. \emph{Biom. J.}, \bold{32}, 79--86 
}
\keyword{datasets}