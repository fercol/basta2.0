\name{summary.bastaCheckCMR}
\alias{summary.bastaCheckCMR}
\alias{summary.bastaCheckCens}
\alias{print.bastaCheckCMR}
\alias{print.bastaCheckCens}

\title{
 Summary of outputs from the data checking function in BaSTA.
}

\description{
 These function provide general information about the input data set to be used for BaSTA, as well as detailes on potential issues that could have been found on the data. 
}

\usage{

\method{summary}{bastaCheckCMR}(object, \dots)

\method{summary}{bastaCheckCens}(object, \dots)

\method{print}{bastaCheckCMR}(x, \dots)

\method{print}{bastaCheckCens}(x, \dots)
}

\arguments{
	\item{object }{Output from function \code{\link{DataCheck}}}
	
	\item{x }{Output from function \code{\link{DataCheck}}}
	
	\item{\dots }{Additional arguments passed to functions \code{print} and \code{summary} (see details).}

}

\details{
Both functions print to the screen the results of the \code{\link{DataCheck}} function.
}


\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{DataCheck}} for running a data check on the input data for function \code{\link{basta}}, \code{\link{FixCMRdata}} to fix potential issues for capture-mark-recapture data.}

\examples{
## CMR data:
## --------- #
## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)

## census data:
## ------------ #
## Load data:
data("bastaCensDat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(object = bastaCensDat, dataType = "census")

## Printed output:
## --------------- #
## Print DataCheck results:
print(checkedData)
}

\keyword{FILL UP}

