\name{summary.multibasta}
\alias{summary.multibasta}
\alias{print.multibasta}
\alias{coef.multibasta}
\title{
Summarizing multiple Bayesian Survival Trajectory Analsys (BaSTA) model output.
}

\description{
These functions allow users visualize the output from function \code{multibasta}.}

\usage{
\method{summary}{multibasta}(object, \dots)
\method{print}{multibasta}(x, \dots)
\method{coef}{multibasta}(object, showAll = FALSE, \dots)
}

\arguments{
  
\item{object }{An object of class \dQuote{\code{multibasta}}.}
 
\item{x }{An object of class \dQuote{\code{multibasta}}.}

\item{showAll }{Logical whether the coefficients from all models should be shown.}

\item{\dots }{Additional arguments to be passed to function \code{print}, \code{summary}, and \code{coef} (see \code{details}).}

}

\value{  
Function \code{coef.multibasta} returns a \code{list} with the following objects:

\item{coefficients }{A list with the \code{basta} coefficients for the model with lowes DIC or, if \code{showAll = TRUE}, from each model and shape tested.}
	
\item{DICs }{A matrix with the DIC values for each model, sorted from the model with lowest DIC to the model with highest DIC.}
	
}

\author{
Fernando Colchero \email{colchero@imada.sdu.dk}, Owen R. Jones \email{jones@biology.sdu.dk}, and Maren Rebke \email{maren.rebke@avitec-research.de}
}

\seealso{
\code{\link{basta}}, \code{\link{multibasta}}.
}

\examples{
\donttest{
## Load data:
data("bastaCMRdat", package = "BaSTA")

## Run short version of BaSTA on the data:
multiout <- multibasta(bastaCMRdat, studyStart = 51, studyEnd = 70, 
                 models = c("GO", "WE"), shapes = "simple",
                 niter = 100, burnin = 11, thinning = 10, 
                 nsim = 1, updateJumps = FALSE, parallel = FALSE)

## Print results:
summary(multiout, digits = 3)

}
}

\keyword{ methods }
