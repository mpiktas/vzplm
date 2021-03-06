coef.pgmm <- function(object,...){
  model <- ifelse(is.null(object$call$model), "onestep", object$call$model)
  if(model=="onestep"){
    coefficients <- object$coefficients
  }
  else{
    coefficients <- object$coefficients[[2]]
  }
  coefficients
}
summary.pgmm <- function(object, robust = FALSE, time.dummies = FALSE, ...){
  model <- describe(object, "model")
  effect <- describe(object, "effect")
  transformation <- describe(object, "transformation")
  if (robust){
    vv <- vcovHC(object)
  }
  else{
    vv <- vcov(object)
  }
  if (model == "onestep")   K <- length(object$coefficients)
  else  K <- length(object$coefficients[[2]])
  Kt <- length(object$args$namest)
  if (!time.dummies && effect == "twoways") rowsel <- -c((K-Kt+1):K)
  else rowsel <- 1:K
  std.err <- sqrt(diag(vv))
  b <- coef(object)
  z <- b / std.err
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable <- CoefTable[rowsel, , drop = FALSE]
  object$sargan <- sargan(object)
  object$m1 <- mtest(object, 1, vv)
  object$m2 <- mtest(object, 2, vv)
  object$wald.coef <- wald(object, "param", vv)
  if (describe(object, "effect") == "twoways") object$wald.td <- wald(object,"time",vv)
  class(object) <- "summary.pgmm"
  object
}

sargan <- function(object){
  model <- describe(object, "model")
  transformation <- describe(object, "transformation")
  if (model == "onestep") Ktot <- length(object$coefficient)
  else Ktot <- length(object$coefficient[[2]])
  z <- mapply(function(x,y) t(x) %*% y,
              object$W,
              residuals(object),
              SIMPLIFY = FALSE)
  z <- Reduce("+", z)
  
  p <- ncol(object$W[[1]])
  if (model == "onestep") A <- object$A1
  else A <- object$A2
  stat <- as.numeric(crossprod(z, t(crossprod(z, A))))
  parameter <- p-Ktot
  names(parameter) <- "df"
  names(stat) <- "chisq"
  method <- "Sargan test"
  pval <- pchisq(stat, df = parameter, lower.tail = FALSE)
  sargan <- list(statistic = stat,
                 p.value = pval,
                 parameter = parameter,
                 method = "Sargan Test")
  class(sargan) <- "htest"
  sargan
}


mtest <- function(object, order = 1, vcov = NULL){
  myvcov <- vcov
  if (is.null(vcov)) vv <- vcov(object)
  else if (is.function(vcov)) vv <- myvcov(object)
  else vv <- myvcov
  model <- describe(object, "model")
  transformation <- describe(object, "transformation")
  Kt <- length(object$args$namest)
  if (transformation == "d"){
    resid <- object$residuals
    residl <- lapply(resid,
                     function(x) c(rep(0,order), x[1:(length(x)-order)])
                     )
  }
  else{
    resid <- lapply(object$residuals,
                    function(x) c(x[-c(Kt:(2*Kt + 1))], rep(0, Kt)))
    residl <- lapply(object$residuals,
                     function(x) c(rep(0, order), x[1:(Kt-order-1)], rep(0, Kt)))
  }
  X <- lapply(object$model, function(x) x[,-1])
  W <- object$W
  if (model == "onestep") A <- object$A1
  else  A <- object$A2
  EVE <- Reduce("+",
                mapply(function(x, y) t(y) %*% x %*% t(x) %*%y, resid, residl, SIMPLIFY = FALSE))
  EX <- Reduce("+", mapply(crossprod, residl, X, SIMPLIFY = FALSE))
  XZ <- Reduce("+", mapply(crossprod, W, X, SIMPLIFY = FALSE))
  ZVE <- Reduce("+",
                mapply(function(x,y,z) t(x)%*%y%*%t(y)%*%z, W, resid, residl, SIMPLIFY = FALSE))

  denom <- EVE - 2 * EX %*% vcov(object) %*% t(XZ) %*% A %*% ZVE + EX %*% vv %*% t(EX)
  num <- Reduce("+", mapply(crossprod, resid, residl, SIMPLIFY = FALSE))
  stat <- num / sqrt(denom)
  names(stat) <- "normal"
  pval <- pnorm(abs(stat), lower.tail = FALSE)
  mtest <- list(statistic = stat,
                p.value = pval,
                method = paste("Autocorrelation test of degree",order))
  class(mtest) <- "htest"
  mtest
}

wald <- function(object, param = "coef", vcov = NULL){
  myvcov <- vcov
  if (is.null(vcov)) vv <- vcov(object)
  else if (is.function(vcov)) vv <- myvcov(object)
  else vv <- myvcov
  model <- describe(object, "model")
  effect <- describe(object, "effect")
  if (param == "time" && effect == "individual") stop("no time-dummies in this model")
  transformation <- describe(object, "transformation")
  if (model == "onestep") coefficients <- object$coefficients
  else coefficients <- object$coefficients[[2]]
  Ktot <- length(coefficients)
  Kt <- length(object$args$namest)
  if (param == "time"){
    start <- Ktot - Kt + ifelse(transformation == "ld", 2, 1)
    end <- Ktot
  }
  else{
    start <- 1
    if (effect == "twoways") end <- Ktot-Kt else end <- Ktot
  }
  coef <- coefficients[start:end]
  vv <- vv[start:end, start:end]
  stat <- t(coef)%*%solve(vv)%*%coef
  names(stat) <- "chisq"
  parameter <- length(coef)
  pval <- pchisq(stat, df = parameter, lower.tail = FALSE)
  wald <- list(statistic = stat,
               p.value = pval,
               parameter = parameter,
               method = "Wald test")
  class(wald) <- "htest"
  wald
}

print.summary.pgmm <- function(x, digits = max(3, getOption("digits") - 2),
                               width = getOption("width"),
                               ...){
  model <- describe(x, "model")
  transformation <- describe(x, "transformation")
  effect <- describe(x, "effect")

  pdim <- attr(x,"pdim")
  formula <- x$call$formula

  cat(paste(effect.pgmm.list[effect]," ",sep=""))
  cat(paste(model.pgmm.list[model],"\n",sep=""))
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  print(pdim)
  ntot <- sum(unlist(x$residuals) != 0)
  cat("\nNumber of Observations Used: ",ntot,"\n")
  
  cat("\nResiduals\n")
  print(summary(unlist(residuals(x))))
  cat("\nCoefficients\n")
  printCoefmat(x$CoefTable,digits=digits)

  cat("\nSargan Test: ",names(x$sargan$statistic),
      "(",x$sargan$parameter,") = ",x$sargan$statistic,
      " (p.value=",format.pval(x$sargan$p.value,digits=digits),")\n",sep="")

  cat("Autocorrelation test (1): ",names(x$m1$statistic),
      " = ",x$m1$statistic,
      " (p.value=",format.pval(x$m1$p.value,digits=digits),")\n",sep="")
  
  cat("Autocorrelation test (2): ",names(x$m2$statistic),
      " = ",x$m2$statistic,
      " (p.value=",format.pval(x$m2$p.value,digits=digits),")\n",sep="")
  cat("Wald test for coefficients: ",names(x$wald.coef$statistic),
      "(",x$wald.coef$parameter,") = ",x$wald.coef$statistic,
      " (p.value=",format.pval(x$wald.coef$p.value,digits=digits),")\n",sep="")
  
  
  if (describe(x, "effect") == "twoways"){
    cat("Wald test for time dummies: ",names(x$wald.td$statistic),
        "(",x$wald.td$parameter,") = ",x$wald.td$statistic,
        " (p.value=",format.pval(x$wald.td$p.value,digits=digits),")\n",sep="")
  }
  invisible(x)
}

fitted.pgmm <- function(object,output=c("pseries","pdata.frame"),...) {
    
    output <- match.arg(output)
     
    if(object$args$model=="twosteps") coeffs <- object$coefficients[[2]]
    else coeffs <- object$coefficients
    
    if(object$args$effect=="twoways") {
        notd <- length(object$arg$namest)
        ncoeff <- length(coeffs)
        td <- diag(1,notd)
        rownames(td) <- object$arg$namest
        
        prodXc <- mapply(function(x) {
            X <- cbind(x[,-1],matrix(NA,nrow=nrow(x),ncol=notd))
            tdX <- intersect(rownames(x),rownames(td))
            X[tdX,ncoeff-notd:1+1] <- td[tdX,]
            X <- diff(X)
            crossprod(t(X),coeffs)
        },object$data)
    }
    else {
        prodXc <- mapply(function(x)crossprod(t(diff(x[,-1])),coeffs),object$data)
    }

    fit <- mapply(function(x,y){
        yy <- y[rownames(x),1,drop=FALSE]
        yy <- rbind(NA,yy[-nrow(yy),1,drop=FALSE])
        r <- matrix(NA,nrow=nrow(y),ncol=1)
        rownames(r) <- rownames(y)
        r[rownames(x),1] <- x+yy
        r
    },prodXc,object$data)
    
    result <- ldply(fit,function(l)data.frame(time=rownames(l),value=l[,1]))
    
    index <- colnames(object$index)
    colnames(result)[1:2] <- index
    result <- pdata.frame(result)
 
    if(output == "pseries") result <- result[,3]

    result
}

predict.pgmm <- function(object,data,horizon=NULL,
        inverse=function(x)x,
        output=c("pseries","pdata.frame"),
        index=NULL,...) {
    if(missing(data) | is.null(horizon)) {
        fitted(object,output=output,...)
    }
    else {
        forecast.pgmm(object,data,horizon,inverse,output,index,...)
    }
    
}
