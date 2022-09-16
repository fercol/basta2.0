\name{print.bastaCheckCMR}
\alias{print.bastaCheckCMR}
\alias{print.bastaCheckCens}

\title{
 FILL UP
}

\description{
 FILL UP
}

\usage{
print(bastacheck)

\method{bastacheck}{bastaCheckCMR}(bastacheck)

\method{bastacheck}{bastaCheckCens}(bastacheck)
}

\arguments{
	\item{bastacheck }{Output from function \code{\link{DataCheck}}}

}

\details{
The function prints to the screen the results of the \code{\link{DataCheck}} function.
}


\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{\code{\link{DataCheck}} for running a data check on the input data for function \code{\link{basta}}, \code{\link{FixCMRdata}} to fix potential issues for capture-mark-recapture data.}

\examples{
## CMR data:
## --------- #
## Load data:
data("bastaCMRdat", package = "BaSTA")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)

## census data:
## ------------ #
## Load data:
data("bastaCensDat", package = "BaSTA")

## Check data consistency:
checkedData  <- DataCheck(object = bastaCensDat, dataType = "census")

## Printed output:
## --------------- #
## Print DataCheck results:
print(checkedData)
}

\keyword{FILL UP}

