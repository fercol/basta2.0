\name{DataCheck}
\alias{DataCheck}

\title{
 Error checking for BaSTA input data.
}

\description{
 A function to check the input data file for a Bayesian Survival Trajectory Analysis (BaSTA) for capture-mark-recapture (CMR) or census data.
}

\usage{DataCheck (object, dataType = "CMR", studyStart = NULL, studyEnd = NULL,  silent = TRUE)}

\arguments{
	\item{ object }{A \code{data.frame} to be used as an input data file for BaSTA. Note: BaSTA can take two types of datasets, namely capture-mark-recapture (\code{CMR}) or census data.}

	\item{dataType }{A \code{character string} indicating if the data are capture-mark-recapture (CMR) or census. Options are \dQuote{\code{CMR}} (default) or \dQuote{\code{census}}.}

	\item{studyStart }{Only required for \code{dataType =} \dQuote{\code{CMR}}, an integer indicating the first year of the study.}

	\item{studyEnd }{Only required for \code{dataType =} \dQuote{\code{CMR}}, an integer indicating the last year of the study.}

	\item{ silent }{Logical to indicate whether the results should be printed to the console.}

}

\details{
The function checks for inconsistencies in the dataset and reports them back. See \code{value} section for details on the types of errors detected by the function. 

\code{DATA SPECIFICATIONS:}

\code{1) CMR data:}
The input data \code{object} requires the following structure: the first column should be a vector of individual unique IDs, the second and third columns are birth and death years respectively. Columns \eqn{4, \dots, T+3} represent the observation window (i.e., recapture matrix) of \eqn{T} years.  This is followed (optionally) by columns for categorical and continuous covariates.

\code{2) census data:}
The input data \code{object} requires at least five dates columns, namely \dQuote{Birth.Date}, \dQuote{Min.Birth.Date}, \dQuote{Max.Birth.Date}, \dQuote{Entry.Date}, and \dQuote{Depart.Date}. All dates need to be format as \dQuote{\%Y-\%m-\%d}. In addition, a \dQuote{Depart.Type} column is required with two types of departures \dQuote{C} for Censored and \dQuote{D} for dead.
}

\value{
  \code{1) CMR data:}
  \item{newData }{The original data frame (for consistency with previous versions of BaSTA).}
  \item{type1 }{A vector of row numbers in the original data frame where there are deaths occurring before the study starts.}
  \item{type2 }{A vector of row numbers in the original data frame where there are no birth/death AND no obervations.}
  \item{type3 }{A vector of row numbers in the original data frame where there are births recorded after death.}
  \item{type4 }{A vector of row numbers in the original data frame where there are observations (i.e. recaptures) after death.}
  \item{type5 }{A vector of row numbers in the original data frame where there are observations (i.e. recaptures) before birth.}
  \item{type6 }{A vector of row numbers in the original data frame where the year of birth is not a zero in the recapture matrix.}
 	\item{summary }{List with summary information, e.g., sample size, number of records with known birth, number of records with known death, etc. }
 	\item{stopExec }{Logical that indicates if the data are free of errors or not. i.e. \code{TRUE} = the data have no apparent errors, and \code{FALSE} = there is at leat one error.}
 	\item{probDescr }{Character vector explaining the six types of problems the \code{\link{DataCheck}} functions looks for.}
 	\item{dataType }{Type of dataset, i.e., \dQuote{\code{CMR}}.}
 	\item{studyStart }{Integer indicating the study start time.}
 	\item{studyEnd }{Integer indicating the study end time.}
 	
	\code{2) census data:}
	\item{n }{Integer for the number of rows (i.e., records) in the dataset.}
 	\item{stopExec }{Logical that indicates if the data are free of errors or not. i.e. \code{TRUE} = the data have no apparent errors, and \code{FALSE} = there is at leat one error.}
 	\item{nas }{List organised by column indicating whether \code{NA}s were detected in a given column.}
 	\item{DateRan }{Matrix of dates ranges (as character strings) for each date column in the dataset.}
 	\item{probDescr }{Character vector explaining the seven types of problems the \code{\link{DataCheck}} functions looks for.}
 	\item{MinBBirth }{Vector of indices of rows where \dQuote{\code{Min.Birth.Date}} was larger than \dQuote{\code{Birth.Date}}.}
 	\item{BirthMaxB }{Vector of indices of rows where \dQuote{\code{Birth.Date}} was larger than \dQuote{\code{Max.Birth.Date}}.}
 	\item{MinBMaxB }{Vector of indices of rows where \dQuote{\code{Min.Birth.Date}} was larger than  \dQuote{\code{Max.Birth.Date}}.}
 	\item{BirthEntr }{Vector of indices of rows where \dQuote{\code{Birth.Date}} was larger than \dQuote{\code{Endtry.Date}}.}
 	\item{MinBEntr }{Vector of indices of rows where \dQuote{\code{Min.Birth.Date}} was larger than  \dQuote{\code{Entry.Date}}.}
 	\item{MaxBEntr }{Vector of indices of rows where \dQuote{\code{Max.Birth.Date}} was larger than  \dQuote{\code{Entry.Date}}.}
 	\item{EntrDep }{Vector of indices of rows where \dQuote{\code{Entry.Date}} was larger than  \dQuote{\code{Depart.Date}}.}
 	\item{DepartType }{Vector of indices of rows where \dQuote{\code{Depart.Type}} does not fall within the \dQuote{\code{C}} (i.e., censored) or \dQuote{\code{D}} (i.e., uncensored or death) categories.}
 	\item{idUnCens }{Vector of indices of rows for uncensored (i.e., death) records.}
 	\item{nUnCens }{Integer indicating the number of uncensored records.}
 	\item{idCens }{Vector of indices of rows for censored records.}
 	\item{nCens }{Integer indicating the number of uncensored records.}
 	\item{idNoBirth }{Vector of indices of rows for records with uncertain birth date.}
 	\item{nNoBirth }{Integer indicating the number of records with uncertain birth date.}	
}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{\code{\link{FixCMRdata}} to fix potential issues for capture-mark-recapture data.}

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

