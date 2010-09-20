\name{fitted.pgmm}
\alias{fitted.pgmm}
\alias{predict.pgmm}
\alias{forecast.pgmm}
\title{Extract Fitted Values or Predict Dynamic Panel Regression Model}

\description{Extract fitted or predict values based on dynamic panel regression object}

\usage{
\method{predict}{pgmm}(object,data,
        horizon=NULL,
        inverse=function(x)x,
        output=c("pseries","pdata.frame"),
        index=NULL,...)

\method{fitted}{pgmm}(object,output=c("pseries","pdata.frame"),...)

forecast.pgmm(object, data, horizon,
             inverse = function(x)x, 
             output = c("pseries","pdata.frame"),
             index=NULL)
}

\arguments{
		\item{object}{an object of class \code{"pgmm"}}
		\item{output}{how the result is presented. If \code{output="pseries"} the fitted values are returned as a \code{pseries} object with index of data supplied to pgmm. If \code{output="pdata.frame"} then \code{pdata.frame} is returned, with first two columns representing the index of data supplied to pgmm, and the third column the fitted values.}
		\item{data}{data frame in which to look for variables with
          which to predict}
		\item{horizon}{the time horizon for the forecast. Must be
		character with values present in pdim attribute of data, when it
		is transformed to pdata.frame. First and the last values are
		only considered to form the time horizon without gaps for forecasting}
		\item{inverse}{if the left-hand side variable is transformed for
		estimation, the
		inverse transformation needs to be supplied}
		\item{index}{the indexes}
		\item{\dots}{further arguments.}
}
\value{pseries or pdata.frame object containing fitted or predicted values}	
	
\details{Method predict acts as a wrapper method. It calls fitted.pgmm if no data is supplied to retrieve fitted values and forecast.pgmm if data is supplied with which to predict.

Usual procedure for forecasting in the reggression is to use conditional
expectation. The forecast (or predicted value) of y_t is then
E(y_t|x_t,x_{t-1},...), where x_t are the left hand side variables of
level equation. 

Since dynamic panel does not estimate the individual effects, it
presents the problem how the predicted values should be retrieved. The
approach taken here is to use difference equation for forecasting. This
means that instead of predicting y_t we predict y_t-y_{t-1} and get the
predicted value of y_t by adding actual value of y_{t-1}.

When predicting with new data iterative procedure is used.
The first iteration produces
the forecast of y_{t+1} which is used then as a variable in the
right-hand side for forecasting the value of y_{t+2} and so on.
		
}

\author{Vaidotas Zemlys}

\examples{
data("EmplUK", package = "plm")

## Arellano and Bond (1991), table 4b 
z1 <- pgmm(log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1)
           + log(capital) + lag(log(output), 0:1) | lag(log(emp), 2:99),
            data = EmplUK, effect = "twoways", model = "twosteps")

fz1 <- fitted(z1)
fz1
femp <- cbind(EmplUK,fz1)

fz1 <- fitted(z1,result = "pdata.frame")
femp <- merge(EmplUK,fz1)

##Since we log the variable emp, we need to supply exp as an inverse
##horizon=as.character(c(1980,1982)) will produce similar results

fz1 <- predict(z1,EmplUK,horizon=as.character(1980:1982),inverse=exp,output="pdata.frame")

colnames(fz1)[3] <- "femp"
femp <- merge(EmplUK,fz1)

}