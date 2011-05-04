\name{fitted.plm}
\alias{fitted.plm}
\alias{predict.plm}
\alias{forecast.plm}
\title{Extract Fitted or Predict Panel Regression Model}
\description{Extract fitted or predict values based on panel regression object}
\usage{
\method{predict}{plm}(object,newdata,
        horizon=NULL,
        inverse=function(x)x,
        output=c("pseries","pdata.frame"),
        index=NULL,...)

\method{fitted}{plm}(object,model=NULL,output=c("pseries","pdata.frame"),...)

forecast.plm(object, newdata, horizon,
             inverse = function(x)x, 
             output = c("pseries","pdata.frame"),
             index=NULL,levels,...)
}

\arguments{
  \item{object}{an object of class \code{"plm"}}
  \item{model}{in what format the result is presented. See
  \code{details}.}
  \item{output}{how the result is presented. If \code{output="pseries"}
  the fitted values are returned as a \code{pseries} object with index
  of data supplied to pgmm. If \code{output="pdata.frame"} then
  \code{pdata.frame} is returned, with first two columns representing
  the index of data supplied to pgmm, and the third column the fitted
  values.}
  \item{newdata}{An optional data frame in which to look for variables
    with which to predict.}
  \item{horizon}{the time horizon for the forecast. Must be character
    with values present in pdim attribute of data.}
  \item{inverse}{if the left-hand side variable is transformed for
  estimation, the inverse transformation needs to be supplied}
  \item{index}{the indexes}
  \item{levels}{meaningful if a panel data model are estimated in first
  diffecences. If FALSE then results are returned in first differences
  and in levels otherwise.}
  \item{\dots}{further arguments.}
}
\details{Method predict works as a wrapper method. It calls fitted.plm if no newdata is supplied or \code{horizon=NULL} and forecast.plm if newdata is supplied with which to predict.

The function forecast.plm is only supported for ''pooling'', ''random''
and ''within'' panel data models.

An argument \code{model} is supplied only for ''random'' and ''within''
panel data models. It sets the format of results: should it be
returned in ''random'' or ''within'' transformation. For example, if we
have a ''within'' panel data model and an argument \code{model="within"}
then results are returned in ''within'' transformation. Otherwise
results are not transformed.

In case of dinamic panel data model(models with lagged right-hand side
endogenous variable) iterative procedure is used. The first iteration
produces the forecast of y_{t+1} which is used then as a variable in the
right-hand side for forecasting the value of y_{t+2} and so on.
}

\value{pseries or pdata.frame object containing fitted or predicted values}
\author{EM}
\seealso{\code{\link{plm}},\code{\link{predict.pgmm}}.}
\examples{
data("Produc", package = "plm")

zz <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data =
Produc, index = c("state","year"))

zz.pred<-predict(zz,newdata=Produc,horizon=as.character(1974:1975),output="pdata.frame")
}
