library(MASS)
library(stabs)
library(glmnet)
library(SuperLearner)
library(randomForest)
library(sparseSVM)
library(pROC)
library(dplyr)
library(tidyverse)
library(readr)
library(made4)
library(parallel)
library(doParallel)
library(foreach)
library(purrr)
library(Boruta)

#n, sample size
#p0=5 (fixed), number of nonzero beta
#p1, number of covariates
#rho=0 if X is independent, rho=0.5 if X is correlated
#threshold=0.65 (fixed), threshold probability for stability selection
#iter=100 (fixed), number of iteration
#B=50 (fixed), 50*2, total 100 subsamples
#train.pcntg=0.6 (fixed), percentage of data in the training set

glm_wrapper=function(n, p0, p1, rho, threshold, iter, B, train.pcntg){
  cl <- parallel::makeCluster(detectCores()-1, setup_strategy = "sequential")
  # Activate cluster for foreach library
  registerDoParallel(cl)
  
  unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  
  stabel_iter=foreach(i = 1:iter, .packages = c("MASS","stabs", "sparseSVM", "Boruta", "made4", "SuperLearner", "pROC", "glmnet")) %dopar%{
    set.seed(1000+i)
    sigma=rho^abs(matrix(c(1:p1),p1,p1,byrow = T) - matrix(c(1:p1),p1,p1,byrow = F))#AR(1) covariance
    x=matrix(mvrnorm(n,rep(0,p1),sigma),n,p1) #without intercept
    #b = matrix(0,p1,1)
    #pos_b=sample(1:p1,p0, replace=FALSE)
    #b[pos_b] = runif(p0,-2,2)
    b <- if (rho == 0) { #if rho=0, the true betas are fixed
      b_matrix <- matrix(0, p1, 1)
      b_matrix[1:5] <- c(1.8, 1.2, 0.5, -1.1, -1.9)
      b_matrix
    } else { #if rh0 !=0, then position of betas are random and true betas are generated from uniform dist 
      b_matrix <- matrix(0, p1, 1)
      pos_b <- sample(1:p1, p0, replace = FALSE)
      b_matrix[pos_b] <- runif(p0, -2, 2)
      b_matrix
    }
    prob=(exp(x%*%b))/(1+exp(x%*%b))
    y=rbinom(n,1,prob)
    df=data.frame(x,y)
    train=sample(1:n, n*train.pcntg, replace=FALSE)
    test=(-train)
    df.train=df[train,]
    df.test=df[-train,]
    
    ############## stability selection with LASSSO
    set.seed(2000+i)
    stab.rconcave=stabsel(x=df[train,-(p1+1)], y=df[train,p1+1],fitfun=glmnet.lasso, cutoff=threshold, PFER = 1, sampling.type = "SS", assumption="r-concave")
    sel.lasso=data.frame(stab.rconcave$selected)$stab.rconcave.selected
    
    ##############stability selection with Sparse SVM
    sel.var.svm.vec=c() 
    for (j in 1:B){
      set.seed(2000+j)
      train.svm=sample(1:length(train), length(train)*0.5, replace=FALSE) # subsampling the training dataset
      yy.svm1=df.train[train.svm, (p1+1)]
      xx.svm1=as.matrix(df.train[train.svm, -(p1+1)])
      svm.mod1=sparseSVM(xx.svm1,yy.svm1, alpha=1, gamma=0.6,lambda=0.20, max.iter=1000, dfmax=p1+1) #dfmax=p+1
      d1= as.vector(coef(svm.mod1)[-1]) #all estimated coefficients, excluding intercept
      sel.var.svm1=which(d1!=0) #non-zero coefficients
      
      yy.svm2=df.train[-train.svm,p1+1]
      xx.svm2=as.matrix(df.train[-train.svm, -(p1+1)])
      svm.mod2=sparseSVM(xx.svm2,yy.svm2, alpha=1, gamma=0.6,lambda=0.20, max.iter=1000, dfmax=p1+1)
      d2= as.vector(coef(svm.mod2)[-1]) #all estimated coefficients, excluding intercept
      sel.var.svm2=which(d2!=0) #non-zero coefficients
      sel.var.svm=c(sel.var.svm1, sel.var.svm2)
      sel.var.svm.vec=c(sel.var.svm.vec, sel.var.svm) # keeping the record of selected variables (column)
    }
    sel.var.svm.tab=table(sel.var.svm.vec) #frequency of each selected variable with corresponding column number
    sel.svm=as.numeric(names(sel.var.svm.tab[sel.var.svm.tab > threshold*B*2])) # selected variables (columns)
    
    ##############stability selection with RF
    sel.var.rf.vec=c()
    for (k in 1:B){
      set.seed(2000+k)
      train.rf=sample(1:length(train), length(train)*0.5, replace=FALSE) # subsampling the training dataset
      df.rf1=df.train[train.rf,]
      rf.mod1=Boruta(y ~ ., data = df.rf1, doTrace = 0, ntree = 300, pValue = 0.01)
      sel.var.rf1 = which(rf.mod1$finalDecision == "Confirmed") # keeping the record of selected variables (column)
      
      df.rf2=df.train[-train.rf,]
      rf.mod2=Boruta(y ~ ., data = df.rf2, doTrace = 0, ntree = 300, pValue = 0.01)
      sel.var.rf2 = which(rf.mod2$finalDecision == "Confirmed") # keeping the record of selected variables (column)
      
      sel.var.rf=c(sel.var.rf1, sel.var.rf2)
      sel.var.rf.vec=c(sel.var.rf.vec, sel.var.rf) # combining the selected variables (columns) from he above two models
    }
    
    sel.var.rf.tab=table(sel.var.rf.vec) #frequency of each selected variable with corresponding column number
    sel.rf=as.numeric(names(sel.var.rf.tab[sel.var.rf.tab >threshold*B*2])) # selected variables (columns)
    
    sel.lasso.svm=union(sel.lasso, sel.svm) #taking the union of the selected variables by lasso and sparsesvm
    final.set=union(sel.lasso.svm, sel.rf) # final subset of selected variables by STABEL
    
    ##################################
    true.b=which(b != 0)
    
    tp.stabel <- length(intersect(final.set, true.b)) #True Positives
    fp.stabel <- length(setdiff(final.set, true.b)) #False Positives 
    fn.stabel <- length(setdiff(true.b, final.set)) #False Negatives
    tn.stabel <- length(setdiff(1:p1, union(true.b, final.set))) #True negatives
    
    newX = data.frame(df[test, final.set, drop = FALSE])
    true.out=df[test,p1+1]
    set.seed(2000+i)
    #######################random forest###########################
    sl.model.rf = SuperLearner(Y = df[train, (p1+1)], X = data.frame(df[train,final.set, drop = FALSE]),
                               method="method.NNLS", family = binomial(), SL.library ="SL.randomForest", cvControl = list(V=10))
    sl.pred.rf=predict(sl.model.rf, newX, type="response")
    ROC.sl.rf <- roc(as.vector(true.out), as.vector(sl.pred.rf$pred))
    auc.rf <- auc(ROC.sl.rf)
    ci.rf=ci.auc(ROC.sl.rf, conf.level=0.95, method="bootstrap", boot.n=1000)
    cord.rf=coords(ROC.sl.rf, x="best", input="threshold", best.method="youden")$threshold
    pred.out.rf=ifelse(sl.pred.rf$pred>cord.rf,1,0)
    cutoffs <- c(0.985, 0.95) #for specificity
    myData.rf <- with(ROC.sl.rf, data.frame(specificities, sensitivities, thresholds))
    sens.rf.spec=lapply(cutoffs, function(cutoff) myData.rf$sensitivities[which.min(abs(myData.rf$specificities-cutoff))])
    t.rf=table(true.out,pred.out.rf)
    sens.rf=t.rf[2,2]/(t.rf[2,2]+t.rf[2,1])
    spec.rf=t.rf[1,1]/(t.rf[1,2]+t.rf[1,1])
    pred.acc.rf <- 1-mean(pred.out.rf != true.out)
    
    #######################SVM###########################
    sl.model.svm = SuperLearner(Y = df[train,p1+1], X = data.frame(df[train, final.set, drop = FALSE]),
                                method="method.NNLS", family = binomial(), SL.library ="SL.svm", cvControl = list(V=10))
    sl.pred.svm=predict(sl.model.svm, newX, type="response")
    ROC.sl.svm <- roc(as.vector(true.out), as.vector(sl.pred.svm$pred))
    auc.svm <- auc(ROC.sl.svm)
    ci.svm=ci.auc(ROC.sl.svm, conf.level=0.95, method="bootstrap", boot.n=1000)
    cord.svm=coords(ROC.sl.svm, x="best", input="threshold", best.method="youden")$threshold
    pred.out.svm=ifelse(sl.pred.svm$pred>cord.svm,1,0)
    myData.svm <- with(ROC.sl.svm, data.frame(specificities, sensitivities, thresholds))
    sens.svm.spec=lapply(cutoffs, function(cutoff) myData.svm$sensitivities[which.min(abs(myData.svm$specificities-cutoff))])
    t.svm=table(true.out,pred.out.svm)
    sens.svm=t.svm[2,2]/(t.svm[2,2]+t.svm[2,1])
    spec.svm=t.svm[1,1]/(t.svm[1,2]+t.svm[1,1])
    pred.acc.svm <- 1-mean(pred.out.svm != true.out)
    
    #########################Ensemble learning####################
    sl.model.com = SuperLearner(Y = df[train,p1+1], X = data.frame(df[train, final.set, drop = FALSE]),
                                method="method.NNLS", family = binomial(), SL.library =c( "SL.randomForest","SL.glm", "SL.svm", "SL.lda"), cvControl = list(V=10))
    sl.pred.com=predict(sl.model.com, newX, type="response")
    ROC.sl.com <- roc(as.vector(true.out), as.vector(sl.pred.com$pred))
    auc.com <- auc(ROC.sl.com)
    ci.com=ci.auc(ROC.sl.com, conf.level=0.95, method="bootstrap", boot.n=1000)
    cord.com=coords(ROC.sl.com, x="best", input="threshold", best.method="youden")$threshold
    pred.out.com=ifelse(sl.pred.com$pred>cord.com,1,0)
    myData.com <- with(ROC.sl.com, data.frame(specificities, sensitivities, thresholds))
    sens.com.spec=lapply(cutoffs, function(cutoff) myData.com$sensitivities[which.min(abs(myData.com$specificities-cutoff))])
    t.com=table(true.out,pred.out.com)
    sens.com=t.com[2,2]/(t.com[2,2]+t.com[2,1])
    spec.com=t.com[1,1]/(t.com[1,2]+t.com[1,1])
    pred.acc.com <- 1-mean(pred.out.com != true.out)
    
    
    ###########################################################
    #####################LASSO with Superlearner###############
    ###########################################################
    set.seed(200+i)
    cvfit=cv.glmnet(as.matrix(df.train[, -(p1+1)]), as.factor(df.train[,(p1+1)]),family="binomial",alpha=1) #finding the best value of lambda
    bestlam=cvfit$lambda.min
    final.set.lasso=which(coef(cvfit, lambda=bestlam)[-1] !=0)
    
    tp.lasso <- length(intersect(final.set.lasso, true.b)) #True Positives
    fp.lasso <- length(setdiff(final.set.lasso, true.b)) #False Positives
    fn.lasso <- length(setdiff(true.b, final.set.lasso)) #False Negatives
    tn.lasso <- length(setdiff(1:p1, union(true.b, final.set.lasso))) #True negatives
    
    
    
    newX.lasso = data.frame(df[test, final.set.lasso, drop = FALSE])
    
    lasso.model.com = SuperLearner(Y = df[train,p1+1], X = data.frame(df[train, final.set.lasso, drop = FALSE]),
                                   method="method.NNLS", family = binomial(), SL.library =c( "SL.randomForest","SL.glm", "SL.svm", "SL.lda"), cvControl = list(V=10))
    lasso.pred.com=predict(lasso.model.com, newX.lasso, type="response")
    ROC.lasso.com <- roc(as.vector(true.out), as.vector(lasso.pred.com$pred))
    lasso.auc.com <- auc(ROC.lasso.com)
    lasso.ci.com=ci.auc(ROC.lasso.com, conf.level=0.95, method="bootstrap", boot.n=1000)
    lasso.cord.com=coords(ROC.lasso.com, x="best", input="threshold", best.method="youden")$threshold
    lasso.pred.out.com=ifelse(lasso.pred.com$pred>lasso.cord.com,1,0)
    myData.lasso <- with(ROC.lasso.com, data.frame(specificities, sensitivities, thresholds))
    sens.lasso.spec=lapply(cutoffs, function(cutoff) myData.lasso$sensitivities[which.min(abs(myData.lasso$specificities-cutoff))])
    lasso.t.com=table(true.out,lasso.pred.out.com)
    lasso.sens.com=lasso.t.com[2,2]/(lasso.t.com[2,2]+lasso.t.com[2,1])
    lasso.spec.com=lasso.t.com[1,1]/(lasso.t.com[1,2]+lasso.t.com[1,1])
    lasso.pred.acc.com <- 1-mean(lasso.pred.out.com != true.out)
    
    ######################################################################
    ###################Superlearner with all variables####################
    ######################################################################
    newX.all = data.frame(df[test, 1:p1])
    
    all.model.com = SuperLearner(Y = df[train,p1+1], X = data.frame(df[train, 1:p1]),
                                 method="method.NNLS", family = binomial(), SL.library =c( "SL.randomForest","SL.glm", "SL.svm", "SL.lda"), cvControl = list(V=10))
    all.pred.com=predict(all.model.com, newX.all, type="response")
    ROC.all.com <- roc(as.vector(true.out), as.vector(all.pred.com$pred))
    all.auc.com <- auc(ROC.all.com)
    all.ci.com=ci.auc(ROC.all.com, conf.level=0.95, method="bootstrap", boot.n=1000)
    all.cord.com=coords(ROC.all.com, x="best", input="threshold", best.method="youden")$threshold
    all.pred.out.com=ifelse(all.pred.com$pred>all.cord.com,1,0)
    myData.all <- with(ROC.all.com, data.frame(specificities, sensitivities, thresholds))
    sens.all.spec=lapply(cutoffs, function(cutoff) myData.all$sensitivities[which.min(abs(myData.all$specificities-cutoff))])
    all.t.com=table(true.out,all.pred.out.com)
    all.sens.com=all.t.com[2,2]/(all.t.com[2,2]+all.t.com[2,1])
    all.spec.com=all.t.com[1,1]/(all.t.com[1,2]+all.t.com[1,1])
    all.pred.acc.com <- 1-mean(all.pred.out.com != true.out)
    
    c(pred.acc.rf, auc.rf,sens.rf,spec.rf,sens.rf.spec[[2]], sens.rf.spec[[1]],
      pred.acc.svm,auc.svm,sens.svm,spec.svm, sens.svm.spec[[2]], sens.svm.spec[[1]],
      pred.acc.com,auc.com,sens.com,spec.com, sens.com.spec[[2]], sens.com.spec[[1]],
      tp.stabel,fp.stabel,fn.stabel,tn.stabel,
      lasso.pred.acc.com, lasso.auc.com, lasso.sens.com, lasso.spec.com, sens.lasso.spec[[2]], sens.lasso.spec[[1]],
      tp.lasso,fp.lasso,fn.lasso,tn.lasso,
      all.pred.acc.com, all.auc.com, all.sens.com, all.spec.com, sens.all.spec[[2]], sens.all.spec[[1]])
    
  }
  pred.sum=as.data.frame(do.call(rbind, stabel_iter))
}

result=glm_wrapper(n=600, p0=5, p1=200, rho=0, threshold=0.65, iter=100, B=50, train.pcntg=0.6)

colnames(result) <- c("model.acc.rf","model.auc.rf","model.sens.rf","model.spec.rf", "model.sens.95.spec.rf", "model.sens.985.spec.rf",
                      "model.acc.svm","model.auc.svm","model.sens.svm","model.spec.svm", "model.sens.95.spec.svm", "model.sens.985.spec.svm",
                      "model.acc.com","model.auc.com","model.sens.com","model.spec.com", "model.sens.95.spec.com", "model.sens.985.spec.com",
                      "tp.stabel","fp.stabel","fn.stabel","tn.stabel",
                      "lasso.pred.acc.com", "lasso.auc.com", "lasso.sens.com", "lasso.spec.com", "lasso.sens.95.spec.com", "lasso.sens.985.spec.com",
                      "tp.lasso","fp.lasso","fn.lasso","tn.lasso",
                      "all.pred.acc.com", "all.auc.com", "all.sens.com", "all.spec.com", "all.sens.95.spec.com", "all.sens.985.spec.com")


