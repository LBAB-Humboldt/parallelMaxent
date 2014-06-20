#This uses k-fold partitioning to divide up the presences, and uses a train and test
#background/pseudoabsence without partitioning.
EvaluatePOModel <- function(folds, covs.pres, covs.bkg.train, covs.bkg.test, mxnt.args){
  results<-data.frame(n.train=rep(0,folds), n.test=0, nparams=0, train.auc=0,
                      test.auc=0, stringsAsFactors=FALSE)
  kvector <- kfold(covs.pres, folds)
  
  for (k in 1:folds){
    n.train <- length(which(kvector!=k))
    n.test <- length(which(kvector==k))
    train.df <- rbind(covs.pres[kvector!=k, ], covs.bkg.train)
    test.df <- rbind(covs.pres[kvector==k, ], covs.bkg.test)
    y.train <- c(rep(1, n.train),rep(0,nrow(covs.bkg.train)))
    y.test <- c(rep(1, n.test),rep(0,nrow(covs.bkg.test)))
    mxnt.obj <- maxent(x=train.df, p=y.train, removeDuplicates=FALSE, args=mxnt.args)
    pocc.train <- predict(mxnt.obj, train.df)
    pocc.test <- predict(mxnt.obj, test.df)
    auc.train <- evaluate(pocc.train[y.train==1], pocc.train[y.train==0])@auc
    auc.test <- evaluate(pocc.test[y.test==1], pocc.test[y.test==0])@auc
    nparams <- sum(getLambdaTable(mxnt.obj@lambdas)[, 2] != 0)
    results[k, ]<-c(n.train, n.test, nparams, auc.train, auc.test)
  }
  return(results)
}

getLambdaTable<-function(lambdas){
  lambdas.list <- strsplit(lambdas,",")
  nparams = length(lambdas) - 4
  varnames=rep("NA",nparams)
  result<-data.frame(lambdas=rep(0,nparams))
  for (i in 1:nparams){
    varnames[i]<-lambdas.list[[i]][1]
    result[i,1]<-as.numeric(lambdas.list[[i]][2])
  }
  result<-data.frame(varnames,result,stringsAsFactors=F)
  return(result)
}


OptimizeLambda <- function(folds, covs.pres, covs.bkg.train, covs.bkg.test, mxnt.args, wd=getwd(), sp.prefix="species"){
  lambda.vector <- c(0.02,0.05,0.1,0.22,0.46,1,2.2,4.6)
  results <-data.frame()
  for(lambda in lambda.vector){
    mxnt.args <- c(mxnt.args,paste0("betamultiplier=",lambda))
    results <- rbind(results, 
                     EvaluatePOModel(folds, covs.pres, covs.bkg.train, covs.bkg.test, mxnt.args))
  }
  results <- cbind(lambda=rep(lambda.vector,each=folds), results)
  write.csv(results, paste0(wd, "/", sp.prefix, "_lambda.optimization.csv"), row.names=FALSE)
  lambda.params <- FindBestLambda(results)
  return(lambda.params)
}

FindBestLambda<-function(df){
  summary.df <- ddply(df, "lambda", summarise, mean.auc=mean(test.auc),
                      median.auc=median(test.auc),mean.nparams=mean(nparams))
  best.lambda <- summary.df$lambda[which.max(summary.df$median.auc)]
  optimum.lambda <- NA
  opt.idx <- NA
  start.idx <- which.max(summary.df$median.auc) + 1
  if(start.idx > nrow(summary.df)){
    return(best.lambda)
  }
  for(i in start.idx:nrow(summary.df)){
    pval <- with(df,
                 wilcox.test(test.auc[lambda == best.lambda], test.auc[lambda == summary.df$lambda[i]],
                             alternative="greater", paired=F)$p.value)
    if(pval<0.05){
      opt.idx <- (i-1)
      optimum.lambda <- summary.df$lambda[opt.idx]
      break
    }
  }
  if(is.na(opt.idx)){
    result=c(best.lambda = best.lambda, 
             best.nparams = summary.df$mean.nparams[(start.idx-1)],
             best.median.auc = max(summary.df$median.auc),
             optimum.lambda = NA,
             optimum.nparams = NA,
             optimum.median.auc = NA)
  } else {
    result=c(best.lambda = best.lambda,
             best.nparams = summary.df$mean.nparams[(start.idx-1)],
             best.median.auc = max(summary.df$median.auc),
             optimum.lambda = optimum.lambda,
             optimum.nparams = summary.df$mean.nparams[opt.idx],
             optimum.median.auc = summary.df$median.auc[opt.idx])
  }
  return(result)
}

 
  
  
