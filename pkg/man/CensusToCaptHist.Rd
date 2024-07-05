\name{CensusToCaptHist}
\alias{CensusToCaptHist}

\title{
Constructs a capture-history matrix from repeated individual observations to be used in Bayesian Survival Trajectory Analysis (BaSTA).
}

\description{
This function takes a table of repeated observations per individual (e.g., discrete interval census), consisting of a pair of vectors with ID and observation date (e.g., year), and converts it into a capture-history matrix to be used as part of a BaSTA CMR data input.
}
\usage{
CensusToCaptHist(ID, d, dformat = "\%Y", timeInt = "Y")
}

\arguments{
  \item{ID}{
A vector with individual IDs.
}
  \item{d}{
A vector of dates when each individual was observed (see \code{details}). 
}
  \item{dformat}{
Defines the date format for \code{d} when \code{d} is of class \code{character} (see \code{details}).}
  \item{timeInt }{
A one character string specifying which time interval should be used between capture occasions. Arguments are \dQuote{\code{Y}} for years, \dQuote{\code{M}} for months, \dQuote{\code{W}} for weeks and \dQuote{\code{D}} for days (see \code{details}).. 
}
}
\details{
The \code{d} argument can be specified as an object of class \code{POSIXct} or \code{POSIXlt}, as a vector of integer time intervals or as a character string indicating the day, month and year. (e.g. \code{dd/mm/yyyy}, \code{mmddyyyy}, \code{mm-dd-yyyy} etc.). When \code{d} is of class \code{character} then argument \code{dformat} needs to be specified using the same conventions as in function \code{format.POSIXct} for objects of class \code{POSIXct} or \code{POSIXlt}.
}

\value{
Creates an object of class \dQuote{\code{data.frame}} with the reconstructed capture recapture matrix.
}

\author{
Owen R. Jones \email{jones@biology.sdu.dk}.
}

\seealso{
\code{\link{DataCheck}} for running a data check on the input data for function \code{\link{basta}}.\code{\link{FixCMRdata}} to fix potential issues for capture-mark-recapture data. \code{\link{bastaCMRdat}} for an example of input CMR datasets. 
}
\examples{
## Create a simulated vector of repeated IDs:
IDvec <- sort(sample(1:5, size = 15, replace = TRUE))

## Simulate dates (e.g., years) of observation per individual:
dVec <- rep(0, length(IDvec))
for(i in unique(IDvec)) {
  svec <- which(IDvec == i)
  dVec[svec] <- sort(sample(1990:1995, length(svec)))
}

## Construct the capture-recapture matrix:
Y <- CensusToCaptHist(ID = IDvec,  d = dVec)
}
