EvaluateBRTModel<-function(folds, covs.pres, covs.bkg.train, covs.bkg.test,brt.params){
  results<-data.frame(n.train=rep(0,folds), n.test=0, ntrees=0, lr=brt.params[2], train.auc=0,
                      test.auc=0, stringsAsFactors=FALSE)
  kvector <- kfold(covs.pres, folds)

  for (k in 1:folds){
    lr.tmp=brt.params[2]
    n.train <- length(which(kvector!=k))
    n.test <- length(which(kvector==k))
    train.df <- rbind(covs.pres[kvector!=k, ], covs.bkg.train)
    test.df <- rbind(covs.pres[kvector==k, ], covs.bkg.test)
    y.train <- c(rep(1, n.train),rep(0,nrow(covs.bkg.train)))
    y.test <- c(rep(1, n.test),rep(0,nrow(covs.bkg.test)))
    
    df <- data.frame(y.train, train.df)
    
    brt.obj <- gbm.step(data=df, gbm.x = 2:ncol(df), gbm.y = 1,
                   family = "bernoulli", tree.complexity = brt.params[1], 
                   learning.rate = lr.tmp, bag.fraction = brt.params[3],prev.stratify=FALSE,
                   site.weights=c(rep(1,n.train), rep(n.train/nrow(covs.bkg.train),nrow(covs.bkg.train))))
    
    while(is.null(brt.obj)){
      lr.tmp=lr.tmp*0.5
      brt.obj <- gbm.step(data=df, gbm.x = 2:ncol(df), gbm.y = 1,
                          family = "bernoulli", tree.complexity = brt.params[1], 
                          learning.rate = lr.tmp, bag.fraction = brt.params[3],prev.stratify=FALSE,
                          site.weights=c(rep(1,n.train), rep(n.train/nrow(covs.bkg.train),nrow(covs.bkg.train))))
    }
    
    pocc.train <- predict(brt.obj, train.df, n.trees=brt.obj$gbm.call$best.trees, type="response")
    pocc.test <- predict(brt.obj, test.df, n.trees=brt.obj$gbm.call$best.trees, type="response")
    
    auc.train <- evaluate(pocc.train[y.train==1], pocc.train[y.train==0])@auc
    auc.test <- evaluate(pocc.test[y.test==1], pocc.test[y.test==0])@auc
    ntrees <- brt.obj$n.trees
    results[k, ]<-c(n.train, n.test, ntrees, lr.tmp, auc.train, auc.test)
  }
  return(results)
}