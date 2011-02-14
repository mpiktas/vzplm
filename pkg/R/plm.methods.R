summary.plm <- function(object,...){
  object$fstatistic <- Ftest(object, test = "F")
  model <- describe(object, "model")
  effect <- describe(object, "effect")
  object$r.squared <- c(rsq  = r.squared(object),
                        adjrsq = r.squared(object, dfcor = TRUE))
  # construct the table of coefficients
  std.err <- sqrt(diag(vcov(object)))
  b <- coefficients(object)
  z <- b/std.err
  p <- 2*pt(abs(z), df = object$df.residual, lower.tail=FALSE)
  object$coefficients <- cbind("Estimate"   = b,
                               "Std. Error" = std.err,
                               "t-value"    = z,
                               "Pr(>|t|)"   = p)
  class(object) <- c("summary.plm", "plm", "panelmodel")
  object
}

print.summary.plm <- function(x,digits= max(3, getOption("digits") - 2),
                              width=getOption("width"),...){
  formula <- formula(x)
  has.instruments <- (length(formula)[2] == 2)
  effect <- describe(x, "effect")
  model <- describe(x, "model")
  cat(paste(effect.plm.list[effect]," ",sep=""))
  cat(paste(model.plm.list[model]," Model",sep=""))
  
  if (model=="random"){
    ercomp <- describe(x, "random.method")
    cat(paste(" \n   (",
              random.method.list[ercomp],
              "'s transformation)\n",
              sep=""))
  }
  else{
    cat("\n")
  }
  if (has.instruments){
    ivar <- describe(x, "inst.method")
    cat(paste("Instrumental variable estimation\n   (",
              inst.method.list[ivar],
              "'s transformation)\n",
              sep=""))
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  pdim <- pdim(x)
  print(pdim)
  if (model == "random"){
    cat("\nEffects:\n")
    print(x$ercomp)
  }
  cat("\nResiduals :\n")
  save.digits <- unlist(options(digits = digits))
  on.exit(options(digits = save.digits))
  print(sumres(x))
  
  cat("\nCoefficients :\n")
  printCoefmat(coef(x), digits = digits)
  cat("\n")
  cat(paste("Total Sum of Squares:    ",signif(tss(x),digits),"\n",sep=""))
  cat(paste("Residual Sum of Squares: ",signif(deviance(x),digits),"\n",sep=""))
  cat(paste("R-Squared      : ", signif(x$r.squared[1], digits),"\n"))
  cat("      Adj. R-Squared : ", signif(x$r.squared[2], digits),"\n")
  fstat <- x$fstatistic
  if (names(fstat$statistic) == "F"){
    cat(paste("F-statistic: ",signif(fstat$statistic),
              " on ",fstat$parameter["df1"]," and ",fstat$parameter["df2"],
              " DF, p-value: ",format.pval(fstat$p.value,digits=digits),"\n",sep=""))
  }
  else{
    cat(paste("Chisq: ",signif(fstat$statistic),
              " on ",fstat$parameter,
              " DF, p-value: ",format.pval(fstat$p.value,digits=digits),"\n",sep=""))
    
  }
  invisible(x)
}

fitted.plm <- function(object, model = NULL,output=c("pseries","pdata.frame"),...){
  # there are two 'models' used ; the fitted model and the
  # transformation used for the fitted values
  fittedmodel <- describe(object, "model")
  if (is.null(model)) model <- fittedmodel
  effect <- describe(object, "effect")
  X <- model.matrix(object, model = model)
  y <- pmodel.response(object, model = model)
  beta <- coef(object)
  if (model == "within" & fittedmodel != "within"){
    Xw <- model.matrix(object, model = "within", effect = effect)
    varwith <- colnames(Xw)
    beta <- beta[varwith]
  }
  if (fittedmodel == "within"){
    if (model == "pooling"){
      if (has.intercept(object)) X <- X[,-1]
      index <- attr(model.frame(object), "index")
      if (effect != "time") id <- index[[1]]
      if (effect != "individual") time <- index[[2]]
      fe <- switch(effect,
                   individual = fixef(object, effect = "individual")[as.character(id)],
                   time = fixef(object, effect="time")[as.character(time)],
                   twoways = fixef(object, effect = "individual")[as.character(id)] +
                   fixef(object, effect = "time")[as.character(time)])
      fv <- as.numeric(crossprod(t(X), beta)) + fe
    }
    if (model == "between"){
      alpha <- mean(y) - crossprod(apply(X[, -1], 2, mean), beta)
      beta <- c(alpha, beta)
      fv <- as.numeric(crossprod(t(X), beta))
    }
    if (model == "within"){
      fv <- as.numeric(crossprod(t(X), beta))
    }
  }
  else{
    fv <- as.numeric(crossprod(t(X), beta))
  }
  fv<-data.frame(attributes(object$model)$index,value=fv)
  fv<-pdata.frame(fv)
  if(output=="pseries")fv<-fv[,"value"]
  fv
}

  
predict.plm <- function(object, newdata = NULL,horizon=NULL,
                        inverse=function(x)x,
                        output=c("pseries","pdata.frame"),
                        index=NULL,...){
  tt <- terms(object)
  result<-if (is.null(newdata) | is.null(horizon)){
    fitted(object,output=output, ...)
  }
  #else{
  #  Terms <- delete.response(tt)
  #  m <- model.frame(Terms, newdata)
  #  X <- model.matrix(Terms, m)
  #  beta <- coef(object)
  #  result <- as.numeric(crossprod(beta, t(X)))
  #}
 # else{
   # forecast.plm(object,newdata,horizon,inverse,output,index)
 # }
  result
}

deviance.panelmodel <- function(object, model = NULL, ...){
  if (is.null(model)) sum(resid(object)^2)
  else sum(residuals(object, model = model)^2)
} 

tss <- function(x, ...){
  UseMethod("tss")
}

tss.default <- function(x){
  var(x)*(length(x)-1)
}

tss.plm <- function(x, model = NULL){
  if (is.null(model)) model <- describe(x, "model")
  effect <- describe(x, "effect")
  if (model == "ht") model = "pooling"
  if (model == "random") theta <- x$ercomp$theta else theta <- NULL
  tss(pmodel.response(x, model = model, effect = effect, theta = theta))
}

r.squared <- function(object, model = NULL,
                      type = c('cor', 'rss', 'ess'), dfcor = FALSE){
  if (is.null(model)) model <- describe(object, "model")
  effect <- describe(object, "effect")
  type <- match.arg(type)
  if (type == 'cor'){
    y <- pmodel.response(object, model = model, effect = effect)
    haty <- fitted(object, model = model, effect = effect)
    R2 <- cor(y, haty)^2
  }
  if (type == 'rss'){
    R2 <- 1 - deviance(object, model = model) / tss(object, model = model)
  }
  if (type == 'ess'){
    haty <- fitted(object, model = model)
    mhaty <- mean(haty)
    ess <- sum( (haty - mhaty)^2)
    R2 <- ess / tss(object, model = model)
  }
  if (dfcor) R2 <- R2 * df.residual(object) / length(resid(object) - 1)
  R2
}
  

residuals.plm <- function(object, model = NULL, ...){
  fittedmodel <- describe(object, "model")
  if (is.null(model)) res <- object$residuals
  else{
    beta <- coef(object)
    effect <- describe(object, "effect")
    X <- model.matrix(object, model = model, effect = effect)
    y <- pmodel.response(object, model = model, effect = effect)
    if (model == "within" & fittedmodel != "within") beta <- beta[-1]
    if (model != "within" & fittedmodel == "within"){
      alpha <- mean(y) - crossprod(apply(X[, -1], 2, mean), beta)
      beta <- c(alpha, beta)
    }
    res <- y - as.numeric(crossprod(t(X), beta))
  }
  res
}

formula.plm <- function(x, ...){
  x$formula
}

# describe function: to extract the characteristics of the plm model
describe <- function(x,
                     what = c('model', 'effect', 'random.method',
                       'inst.method', 'transformation')){
  what <- match.arg(what)
  cl <- x$args
##   if (is.name(cl$effect)) cl$effect <- eval(cl$effect, parent.frame())
##   if (is.name(cl$model)) cl$model <- eval(cl$model, parent.frame())
##   if (is.name(cl$random.method)) cl$random.method <- eval(cl$random.method, parent.frame())
##   if (is.name(cl$inst.method)) cl$inst.method <- eval(cl$inst.method, parent.frame())
  switch(what,
         model  = ifelse(!is.null(cl$model), cl$model, "within"),
         effect = ifelse(!is.null(cl$effect), cl$effect, "individual"),
         random.method = ifelse(!is.null(cl$random.method),
           cl$random.method, "swar"),
         inst.method   = ifelse(!is.null(cl$inst.method),
           cl$inst.method, "bvk"),
         transformation = ifelse(!is.null(cl$transformation),
           cl$transformation, "d")
         )
}
         
  
  
forecast.pooling<-function(object,newdata,horizon,
                       inverse=function(x)x,
                       output=c("pseries","pdata.frame"),
                       levels=FALSE,    
                       index=NULL){

    output <- match.arg(output)
  
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("formula", "data", "subset", "na.action", "index"),names(mf),0)
    mf <- mf[c(1,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("plm")
    mf$model <- NA
    mf$formula <- formula(object$formula)
    mf$na.action <- "na.pass"
    mf$data <- as.name("newdata")
    mf$index <- eval(expression(index),list(index=index))

    #Format the data as pdata.frame   
    if (inherits(newdata, "pdata.frame") && !is.null(index)) 
        warning("the index argument is ignored because data is a pdata.frame")
    if (!inherits(newdata, "pdata.frame")) 
        newdata <- pdata.frame(newdata, index)

    index <- attr(newdata, "index")
    pdim <- pdim(newdata)
    pdim.model <- attr(object,"pdim")
    time <- pdim$panel.names$time.names

    #Check for non-existant times
    horizon <- factor(horizon,levels=pdim$panel.names$time.names)    
    if(sum(is.na(horizon))>0)
      stop("The forecast horizon containts times not present in supplied data")
    
    #Make horizon continuous
    start <- which(time == horizon[1])
    end <- which(time == horizon[length(horizon)])
    if(end<start) stop("The end of forecast horizon is earlier than the start")
    horizon <- time[start:end]
    #######
    max.lag <- max(sapply(dynterms(object$formula),max))+1
    time <- pdim$panel.names$time.names
    time.column <- colnames(index)[2]

    #Identifying if the set of explantory variables have a lagged dependent
    #variable
    
    lagged.var<-names(which(mapply(function(x)x>0,dynterms(object))==TRUE))
    endoname <- all.vars(object$formula[[2]])
    lagged.endo<-length(grep(endoname,lagged.var))>0
    endo<-strsplit(as.character(object$formula),"~")[[2]]
    diff.endo<-length(grep("diff",endo))>0

    i<-start
    while(i<=end){
      
      et <-if(lagged.endo|diff.endo)time[(i-max.lag):(i+1)] 
      else time[(start-max.lag):end]

      fdata <- eval(mf, list(newdata=as.data.frame(newdata[newdata[,time.column]%in%et])))
      attr(fdata,"formula") <- formula(object$formula)
      yX <- extract.data(fdata)
      coeffs<-coef(object)
      intercept<-has.intercept(object$formula)
      if(intercept)yX<- mapply(function(x){
                                 x<-cbind(x,rep(1,dim(x)[1]))
                                 colnames(x)[dim(x)[2]]<-"(Intercept)"
                                 temp<-colnames(x)[-c(1,dim(x)[2])]
                                 x<-cbind(x[,c(1,dim(x)[2])],x[,-c(1,dim(x)[2])])
                                 colnames(x)[-c(1,2)]<-temp
                                 x
                               },
                              yX,SIMPLIFY=FALSE)  
      prodXc <- mapply(function(x)crossprod(t(x[,-1]),coeffs),
                     yX,SIMPLIFY=FALSE)
      if(diff.endo){
        data.split <- split(as.data.frame(newdata), index[[1]])
        time.split <- split(index[[2]], index[[1]])
        data.split <- mapply(
                 function(x, y){
                   rownames(x) <- y
                   x
                 }
                 , data.split, time.split, SIMPLIFY = FALSE)
      
      
        fit <- mapply(function(x,y){
                                  yy <- y[rownames(x),endoname,drop=FALSE]
                                  yy <- rbind(NA,yy[-nrow(yy),1,drop=FALSE])
                                  z<-x+yy
                                  rownames(z)<-rownames(x)
                                  z
                                 },prodXc,data.split,SIMPLIFY=FALSE)
      }
      else fit<-prodXc
      fit <- ldply(fit,function(l)data.frame(time=rownames(l),l))

      if(lagged.endo|diff.endo){
        fit <- fit[fit[,2]==time[i],]
        newdata[newdata[,time.column]==time[i], endoname] <- if(diff.endo)
          fit[,3]
        else inverse(fit[,3])
        i<-i+1
      }
      else{
        newdata[newdata[,time.column]%in%horizon, endoname] <- inverse(fit[fit[,"time"]%in%horizon,3])
        i<-end+1
      }
    }
   
    result <- newdata[,c(colnames(index),endoname)]
    if(!levels) result[,endoname]<-diff(result[,endoname])
    result <- result[result[,2] %in% horizon,]

    if(output == "pseries") result <- result[,3]

    result
}
