\name{stats}
\alias{stats}
\title{In-sample forecasting accuracy statistics for panel model}
\description{
  This function computes the R squared, MASE, MAPE and MAE for plm and pgmm objects.
}
\usage{
  stats(object,freq=NULL,...)
}
\arguments{
  \item{object}{an object of class \code{"plm"} or \code{"pgmm"},}
  \item{freq}{the frequency of data, default is 1,}
  \item{...}{other arguments to be passed to \code{fitted} function.}
}

\details{
  An argument \code{model} allows to define on which transformation of
  the data statistics have to be calculated, default is
  \code{model="pooling"}. See \code{\link{fitted.plm}} for details.   

  \code{R} square is computed using correlation mathod.

  \code{MASE} - mean absolute scaled error.

  \code{MAPE} - mean absolute percentage error.

  \code{MAE} - mean absolute error.
}
\value{
  a list of elements with statistics. Each element represents an individual
  of panel data.
}

\examples{
data("Produc", package = "EMplm")
zz <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, index = c("state","year"))
stats(zz)
}

