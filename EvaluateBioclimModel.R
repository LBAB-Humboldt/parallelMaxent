#covs.pres=occs.covs[sp.idx, ]
#covs.bkg.train=train.bkg
#covs.bkg.test=test.bkg
#EvaluateBRTModel(10,covs.pres=occs.covs[sp.idx, ],covs.bkg.train=train.bkg,covs.bkg.test=test.bkg,lr=0.01)

EvaluateBioclimModel<-function(folds, covs.pres, covs.bkg.train, covs.bkg.test){
  results<-data.frame(n.train=rep(0,folds), n.test=0, train.auc=0,
                      test.auc=0, stringsAsFactors=FALSE)
  kvector <- kfold(covs.pres, folds)

  for (k in 1:folds){
    n.train <- length(which(kvector!=k))
    n.test <- length(which(kvector==k))
    train.df <- rbind(covs.pres[kvector!=k, ], covs.bkg.train)
    test.df <- rbind(covs.pres[kvector==k, ], covs.bkg.test)
    y.train <- c(rep(1, n.train),rep(0,nrow(covs.bkg.train)))
    y.test <- c(rep(1, n.test),rep(0,nrow(covs.bkg.test)))
    
    bc.obj <- bioclim(covs.pres[kvector!=k, ])
     
    pocc.train <- predict(bc.obj, train.df)
    pocc.test <- predict(bc.obj, test.df)
    
    auc.train <- evaluate(pocc.train[y.train==1], pocc.train[y.train==0])@auc
    auc.test <- evaluate(pocc.test[y.test==1], pocc.test[y.test==0])@auc
    results[k, ]<-c(n.train, n.test, auc.train, auc.test)
  }
  return(results)
}