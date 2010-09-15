\name{fitted.pgmm}
\alias{fitted.pgmm}
\title{Extract Fitted Values of panel GMM model}

\description{Extract fitted values of panel GMM model for level equations}

\usage{
		\method{fitted}{pgmm}(object,intercept=c("none","global","individual"),result=c("pseries","pdata.frame"),...)
}

\arguments{
		\item{object}{an object of class \code{"pgmm"}}
		\item{intercept}{the type of intercept for level equation, see Details}
		\item{result}{how the result is presented. If \code{result="pseries"} the fitted values are returned as a \code{pseries} object with index of data supplied to pgmm. If \code{result="pdata.frame"} then \code{pdata.frame} is returned, with first two columns representing the index of data supplied to pgmm, and the third column the fitted values.}
		\item{\dots}{further arguments.}

}
\value{If \code{intercept="none"}, the result is either object of class \code{pseries} or \code{pdata.frame}, depending on argument \code{result}. If \code{intercept!="none"} then the list is returned with the following elements:	
	
\item{fitted.values}{pseries or pdata.frame object containing the fitted values}
\item{intercept}{one constant for all individuals, or list of individual constants}
}

\details{The fitted values of pgmm object is simply calculated by multiplying the right hand side variables by the estimated coefficients. Since the intercept is not estimated due to differencing used is GMM estimation, three options are presented to deal with this problem.

If \code{intercept="none"}, the fitted values are left as is. If \code{intercept="global"} it is asumed that we have one constant for all model. It is calculated by taking mean of the residuals. If \code{intercept="individual"}, it is assumed that we have individual constants for each individual. They are estimated as a mean of residuals for each individual. The list containing individual constants is returned.
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

fz1 <- fitted(z1,result = "pdata.frame",intercept="global")
fz1$intercept

fz1 <- fitted(z1,result = "pdata.frame",intercept="individual")
fz1$intercept

}