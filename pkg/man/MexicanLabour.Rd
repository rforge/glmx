\name{MexicanLabour}
\docType{data}
\alias{participation}
\alias{MexicanLabour}
\title{Mexican Women's Labour-Force Participation}

\description{The data are originally from the National Survey of Household Income and Expenditures for 1977, Secretaria de Programacion y Presupuesto, Mexico.}

\usage{data("MexicanLabour")}

\format{A data frame with 16 rows and 6 columns. 
	\describe{
		\item{total}{integer. number of women older than 12}
		\item{labour}{integer. number of women in labour force}
		\item{locality}{factor. levels: \code{rural} and \code{urban}}
		\item{age}{factor. levels: \code{young} (\eqn{\leq 24} years) and \code{old} (\eqn{> 24} years)}
		\item{income}{factor. levels: \code{low} (household income less than $2626.8 = 0) and \code{high} (more than $2626.8)}
		\item{schooling}{factor. levels: \code{primary} (primary school or less) and \code{further} (more than primary school)} 		}
	}

\details{The data were first analysed by Guerrero and Johnson (1982) as an example of a highly asymmetric data set, i.e. the observed proportions are rather low.}

\source{
Guerrero, V. & Johnson, R. (1982). Use of the Box-Cox transformation with binary response models. \emph{Biometrika}, \bold{69}, 309--314.
}

\keyword{datasets}