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

#check
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
    max.lag <- max(sapply(dynterms(object$formula),max))#+1
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
      
      et <-if(lagged.endo)time[1:(i)]#+1)]i-max.lag 
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
      
      fit<-prodXc
      fit <- ldply(fit,function(l)data.frame(time=rownames(l),l))

      if(lagged.endo){
        fit <- fit[fit[,2]==time[i],]
        if(diff.endo){
          newdata[newdata[,time.column]==time[i], endoname] <- newdata[newdata[,time.column]==time[i-1], endoname] + inverse(fit[,3])
        }
        else
          newdata[newdata[,time.column]==time[i], endoname] <- inverse(fit[,3])
          i<-i+1
      }
      else{
        newdata[newdata[,time.column]%in%horizon, endoname] <- inverse(fit[fit[,"time"]%in%horizon,3])
        i<-end+1
      }
    }
    if(levels & diff.endo & !lagged.endo){
       data.split <- split(as.data.frame(newdata), index[[1]])
       time.split <- split(index[[2]], index[[1]])
       data.split <- mapply(
                 function(x, y){
                   rownames(x) <- y
                   x
                 }
                 , data.split, time.split, SIMPLIFY = FALSE)
       
      data.split <- mapply(function(x){
                                   first <- which(x[,time.column]%in%horizon[1]) - 1
                                   start.endo <-x[first,endoname]
                                   y<-x
                                   for(i in 1:length(horizon)){
                                      temp<-x[x[,time.column]%in%horizon[i],endoname]
                                     y[y[,time.column]%in%horizon[i],endoname]<-if(horizon[i]==horizon[1])start.endo + temp 
                                         else y[y[,time.column]%in%horizon[i-1],endoname]+temp 
                                   }
                                   y
                                 },data.split,SIMPLIFY=FALSE)
       newdata<-ldply(data.split,function(l)data.frame(time=rownames(l),l))
    }
    result <- newdata[,c(colnames(index),endoname)]
    if(!levels&diff.endo&lagged.endo) result[,endoname]<-diff(result[,endoname])
    result <- result[result[,2] %in% horizon,]

    if(output == "pseries") result <- result[,3]

    result
}

forecast.within<-function(object,newdata,horizon,
                       inverse=function(x)x,
                       model="pooling",   
                       output=c("pseries","pdata.frame"),
                       levels=FALSE,    
                       index=NULL){

    output <- match.arg(output)
    effect<-describe(object,"effect")

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
    if(lagged.endo) FIT<-NULL 
    
    while(i<=end){
      
      et <-if(lagged.endo)time[1:(i)]#+1)]i-max.lag 
      else time[(start-max.lag):end]

      fdata <- eval(mf, list(newdata=as.data.frame(newdata[newdata[,time.column]%in%et])))
      attr(fdata,"formula") <- formula(object$formula)
      yX <- extract.data(fdata)
      
      if(model=="within") {
        WITHIN<-trans(yX,fdata,effect,time.column)
        yX <- WITHIN$res
        time.mean<-WITHIN$time.mean
        ind.mean<-WITHIN$individual.mean
        var.mean<-WITHIN$var.mean
      }
      
      coeffs<-if((effect=="twoways"|effect=="time")& model!="within") c(coef(object),fixef(object,effect="time"))
              else coef(object)
    
      if(effect=="twoways"|effect=="time") {
            combined.time <- factor(sort(unique(c(names(fixef(object,effect="time")),pdim(newdata)$panel.names$time.names))))  
            notd <- length(names(fixef(object,effect="time")))
            ncoeff <- length(coeffs)
            td <- diag(1,notd)
            rownames(td) <- names(fixef(object,effect="time"))

            zerotd <- matrix(0,ncol=notd,nrow=length(combined.time))
            rownames(zerotd) <- levels(combined.time)
            colnames(zerotd) <- colnames(td)
            zerotd[rownames(td),] <- td
            td <- zerotd
       }     
       if(model!="within"&(effect=="twoways"|effect=="time")){
              prodXc <- mapply(function(x) {
                 X <- cbind(x[,-1],matrix(NA,nrow=nrow(x),ncol=notd))
                 tdX <- intersect(rownames(x),rownames(td))
                 X[tdX,ncoeff-notd:1+1] <- td[tdX,]
                 crossprod(t(X),coeffs)
              },yX,SIMPLIFY=FALSE)
      }
      else prodXc <- mapply(function(x)crossprod(t(x[,-1]),coeffs),yX,SIMPLIFY=FALSE)
      if(effect=="twoways"|effect=="individual"){

fix.effect<-as.list(fixef(object,effect="individual"))
        if(model!="within"){
                            fit<-mapply(function(x,y){
                                  x+y
                                 },prodXc,fix.effect,SIMPLIFY=FALSE)
                           }else fit<-prodXc
      }
      else fit<-prodXc

      if(lagged.endo){
         temp<-ldply(fit,function(l)data.frame(time=rownames(l),l))
         FIT<-rbind(FIT,temp[temp[,2]==time[i],])
         if(model=="within")fit=switch(effect,
                                        individual=mapply(function(x,y,z){
                                                    ind.part<-as.numeric(crossprod(y[-1],coeffs))
                                                    y<-x+ind.part+z
                                                    y
                                                   },fit,ind.mean,fix.effect,SIMPLIFY=FALSE),
                                        time=mapply(function(x){
                                              time.part<-as.numeric(crossprod(t(time.mean[,-1]),coeffs))
                                              fix.part<-as.numeric(crossprod(t(td[rownames(td)%in%rownames(x),]),fixef(object,"time")))
                                              y<-x+time.part+fix.part
                                              y
                                             },fit,SIMPLIFY=FALSE),
                                        twoways=mapply(function(x,y,z){
                                              ind.part<-as.numeric(crossprod(y[-1],coeffs))
                                              time.part<-as.numeric(crossprod(t(time.mean[,-1]),coeffs))
                                              var.part<-as.numeric(crossprod(t(var.mean[,-1]),coeffs))
                                              fix.part<-as.numeric(crossprod(t(td[rownames(td)%in%rownames(x),]),fixef(object,"time")))     
                                              
                                              y<-x+ind.part+time.part-var.part+fix.part+z
                                              y
                                             },fit,ind.mean,fix.effect,SIMPLIFY=FALSE)
                                         
                                )
      }
      fit <- ldply(fit,function(l)data.frame(time=rownames(l),l))

     if(lagged.endo){
        fit <- fit[fit[,2]==time[i],]
        if(diff.endo){
          newdata[newdata[,time.column]==time[i], endoname] <- newdata[newdata[,time.column]==time[i-1], endoname] + inverse(fit[,3])
        }
        else
          newdata[newdata[,time.column]==time[i], endoname] <- inverse(fit[,3])
          i<-i+1
      }
      else{
        newdata[newdata[,time.column]%in%horizon, endoname] <- inverse(fit[fit[,"time"]%in%horizon,3])
        i<-end+1
      }
    }
  
    if(levels & diff.endo & model!="within" & !lagged.endo){
       data.split <- split(as.data.frame(newdata), index[[1]])
       time.split <- split(index[[2]], index[[1]])
       data.split <- mapply(
                 function(x, y){
                   rownames(x) <- y
                   x
                 }
                 , data.split, time.split, SIMPLIFY = FALSE)
       
      data.split <- mapply(function(x){
                                   first <- which(x[,time.column]%in%horizon[1]) - 1
                                   start.endo <-x[first,endoname]
                                   y<-x
                                   for(i in 1:length(horizon)){
                                      temp<-x[x[,time.column]%in%horizon[i],endoname]
                                     y[y[,time.column]%in%horizon[i],endoname]<-if(horizon[i]==horizon[1])start.endo + temp 
                                         else y[y[,time.column]%in%horizon[i-1],endoname]+temp 
                                   }
                                   y
                                 },data.split,SIMPLIFY=FALSE)
       newdata<-ldply(data.split,function(l)data.frame(time=rownames(l),l))
    }
    result <- newdata[,c(colnames(index),endoname)]
    if(lagged.endo & !levels){
      result <- merge(result[,c(1,2)],FIT,by=c(1,2))
      colnames(result)[3]<-endoname
    }
    result <- result[result[,2] %in% horizon,]

    if(output == "pseries") result <- result[,3]

    result
}


trans<-function(X1,X2,effect,time.column,theta=1,intercept=FALSE){
  if(effect=="twoways" & class(theta)=="list"){
    theta.time<-theta$time
    theta.id<-theta$id
    theta.total<-theta$total
  }
  else{
    theta.time<-theta
    theta.id<-theta
    theta.total<-theta
  }
    
  X2<-cbind(attributes(X2)$index,as.data.frame(X2))
  if(intercept) {
    X2<-cbind(X2,rep(1,dim(X2)[1]))
    colnames(X2)[dim(X2)[2]]<-"(Intercept)"
    temp<-colnames(X2)[-c(1,2,dim(X2)[2])]
    X2<-cbind(X2[,c(1,2,3,dim(X2)[2])],X2[,-c(1,2,3,dim(X2)[2])])
    colnames(X2)[-c(1,2,3)]<-temp
  }  
  time.mean<-aggregate(X2[,-c(1,2)],by=list(X2[,time.column]),function(l)mean(l,na.rm=TRUE))
  colnames(time.mean)<-c(time.column,colnames(X1[[1]]))
  ind.mean<-mapply(function(x){
                             apply(x,2,function(l)mean(l,na.rm=TRUE))
                           },X1,SIMPLIFY=FALSE)
  var.mean<-matrix(apply(X2[,-c(1,2)],2,function(l)mean(l,na.rm=TRUE)),ncol=dim(X2)[2]-2)
  colnames(var.mean)<-colnames(X1[[1]])
  res <- switch(effect,
                individual=mapply(function(x){
                             apply(x,2,function(l)l-theta.id*mean(l,na.rm=TRUE))
                           },X1,SIMPLIFY=FALSE),
                time=mapply(function(x){
                             y<-as.matrix(x-theta.time*time.mean[,-1])
                             rownames(y)<-rownames(x)
                             y
                     },X1,SIMPLIFY=FALSE),
                twoways=mapply(function(x){
                             y<-apply(x,2,function(l)l-theta.id*mean(l,na.rm=TRUE))
                             y<-y-theta.time*time.mean[,-1]
                             z<-as.matrix(y+theta.total*apply(var.mean,2,function(l)rep(l,dim(y)[1])))
                             rownames(z)<-rownames(x)
                             z
                         },X1,SIMPLIFY=FALSE)
                )
  RES<-list(res=res,time.mean=time.mean[,colnames(time.mean)%in%colnames(X1[[1]])],individual.mean=ind.mean,var.mean=var.mean)
  RES
}


forecast.random<-function(object,newdata,horizon,
                       inverse=function(x)x,
                       output=c("pseries","pdata.frame"),
                       levels=FALSE,    
                       index=NULL){

    output <- match.arg(output)
    effect<-describe(object,"effect")

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

    # theta for random transformation

    theta<-object$ercomp$theta
    if(lagged.endo) FIT<-NULL
    i<-start
    while(i<=end){
      
      #et <-if(lagged.endo|diff.endo)time[(i-max.lag):(i+1)] 
      #else time[(start-max.lag):end]
      et <-if(lagged.endo)time[1:(i)]#+1)]i-max.lag 
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
      RANDOM<-trans(yX,fdata,effect,time.column,theta,intercept)
      yX<-RANDOM$res
      time.mean<-RANDOM$time.mean
      ind.mean<-RANDOM$individual.mean
      var.mean<-RANDOM$var.mean
      
      prodXc <- mapply(function(x)crossprod(t(x[,-1]),coeffs),
                     yX,SIMPLIFY=FALSE)
      fit<-prodXc
      
      if(lagged.endo){
         temp<-ldply(fit,function(l)data.frame(time=rownames(l),l))
         FIT<-rbind(FIT,temp[temp[,2]==time[i],])
         fit=switch(effect,
                    individual=mapply(function(x,y,z){
                      ind.part<-as.numeric(crossprod(y[-1],coeffs))
                      y<-x+theta*ind.part
                      y
                    },fit,ind.mean,SIMPLIFY=FALSE),
                    time=mapply(function(x){
                       time.part<-as.numeric(crossprod(t(time.mean[,-1]),coeffs))
                       y<-x+theta*time.part
                       y
                    },fit,SIMPLIFY=FALSE),
                    twoways=mapply(function(x,y){
                      ind.part<-as.numeric(crossprod(y[-1],coeffs))
                      time.part<-as.numeric(crossprod(t(time.mean[,-1]),coeffs))
                      var.part<-as.numeric(crossprod(var.mean[,-1],coeffs))
                      y<-x+theta$id*ind.part+theta$time*time.part-theta$total*var.part
                      y
                    },fit,ind.mean,SIMPLIFY=FALSE)
                    )
      }
      fit <- ldply(fit,function(l)data.frame(time=rownames(l),l))

       if(lagged.endo){
        fit <- fit[fit[,2]==time[i],]
        if(diff.endo){
          newdata[newdata[,time.column]==time[i], endoname] <- newdata[newdata[,time.column]==time[i-1], endoname] + inverse(fit[,3])
        }
        else
          newdata[newdata[,time.column]==time[i], endoname] <- inverse(fit[,3])
          i<-i+1
      }
      else{
        newdata[newdata[,time.column]%in%horizon, endoname] <- inverse(fit[fit[,"time"]%in%horizon,3])
        i<-end+1
      }
    }
   if(levels & diff.endo & !lagged.endo){
       data.split <- split(as.data.frame(newdata), index[[1]])
       time.split <- split(index[[2]], index[[1]])
       data.split <- mapply(
                 function(x, y){
                   rownames(x) <- y
                   x
                 }
                 , data.split, time.split, SIMPLIFY = FALSE)
       
      data.split <- mapply(function(x){
                                   first <- which(x[,time.column]%in%horizon[1]) - 1
                                   start.endo <-x[first,endoname]
                                   y<-x
                                   for(i in 1:length(horizon)){
                                      temp<-x[x[,time.column]%in%horizon[i],endoname]
                                     y[y[,time.column]%in%horizon[i],endoname]<-if(horizon[i]==horizon[1])start.endo + temp 
                                         else y[y[,time.column]%in%horizon[i-1],endoname]+temp 
                                   }
                                   y
                                 },data.split,SIMPLIFY=FALSE)
       newdata<-ldply(data.split,function(l)data.frame(time=rownames(l),l))
    }
    result <- newdata[,c(colnames(index),endoname)]
    if(lagged.endo & !levels){
      result <- merge(result[,c(1,2)],FIT,by=c(1,2))
      colnames(result)[3]<-endoname
    }
    result <- result[result[,2] %in% horizon,]

    if(output == "pseries") result <- result[,3]

    result
}



