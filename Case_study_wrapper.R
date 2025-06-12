t1=Sys.time()
library(stabs)
library(glmnet)
library(randomForest)
library(e1071)
library(sparseSVM)
library(VSURF)
library(SuperLearner)
library(caret)
library(ggplot2)
library(pROC)
library(gtools)
library(Boruta)
setwd("...")
dat=read.csv("end_ovar_lgdata3_with_ca125.csv", header=T)
attach(dat)

xxx=as.matrix(data.frame(lgccl21_pg_ml, lgadiponectin_pg_ml, lgadipsin_pg_ml, lgbca_1_pg_ml, lgccl19_mip3b_pg_ml,
                         lgccl20_mip3a_pg_ml, lgcrp_pg_ml, lgctack_pg_ml, lgcxcl11_i_tac_pg_ml, lgcxcl6_gcp2_pg_ml,
                         lgcxcl9_mig_pg_ml, lgegf_pg_ml,lgena_78_pg_ml, lgeotaxin_pg_ml, lgeotaxin_2_pg_ml, 
                         lgfgf_2_pg_ml, lggro_pg_ml, lgg_csf_pg_ml, lgil_16_pg_ml, lgil_1ra_pg_ml, lgil_29_ifnl1_pg_ml,
                         lgil_33_pg_ml, lgil_7_pg_ml, lgil_8_pg_ml, lgip_10_pg_ml, lglipocalin_2_pg_ml, lgmcp_1_pg_ml,
                         lgmcp_2_pg_ml, lgmcp_4_pg_ml, lgmdc_pg_ml, lgmip_1b_pg_ml, lgmip_1d_pg_ml, lgpai_1_pg_ml,
                         lgresistin_pg_ml, lgsaa_pg_ml,lgsap_pg_ml, lgscf_pg_ml, lgsdf_1a_b_pg_ml, lgsegfr_pg_ml,
                         lgsgp130_pg_ml, lgsilrii_pg_ml, lgsil_4r_pg_ml, lgsil_6r_pg_ml, lgstnfri_pg_ml,
                         lgstnfrii_pg_ml, lgsvegfr2_pg_ml, lgsvegfr3_pg_ml, lgtarc_pg_ml, lgtgf_a_pg_ml, lgtgf_b1_pg_ml,
                         lgtnfa_pg_ml, lgtnf_b_pg_ml, lgtpo_pg_ml, lgtrail_pg_ml, lgtslp_pg_ml, lgvegf_pg_ml,
                         lgamylin_pg_ml, lgc_peptide_pg_ml, lgflt_3l_pg_ml, lgfractalkine_pg_ml,
                         lggip_pg_ml, lgglp_1_pg_ml, lgglucagon_pg_ml, lggm_csf_pg_ml, lgifna2_pg_ml, lgifng_pg_ml,
                         lgil_10_pg_ml, lgil_12p40_pg_ml, lgil_12p70_pg_ml, lgil_15_pg_ml, lgil_17_pg_ml, lgil_1a_pg_ml,
                         lgil_1b_pg_ml, lgil_2_pg_ml, lgil_3_pg_ml, lgil_4_pg_ml, lgil_5_pg_ml, lgil_6_pg_ml,
                         lginsulin_pg_ml, lgleptin_pg_ml,lgmcp_3_pg_ml, lgmip_1a_pg_ml, lgpp_pg_ml, lgpyy_pg_ml,
                         lgscd30_pg_ml, lgscd40l_pg_ml,lgsil_1ri_pg_ml, 
                         lgsil_2ra_pg_ml, lgsrage_pg_ml, lgsvegfr1_pg_ml, lgca125, ovar_cancer, endo_cancer),nrow=869, ncol=92)

##checking percentage of missing values
many_na=(colMeans(is.na(xxx)))*100
##removing variables with more than 60% missing values
df.m.p=data.frame(names(many_na), many_na) #variable names and missing percentage
colnames(df.m.p) <- c("variable", "percentage")
new.df.m.p=subset(df.m.p, df.m.p$percentage<33) #27 biomarkers left
df1=xxx[, new.df.m.p$variable]
df2=na.omit(df1) #868 samples left

df3=as.matrix(subset(data.frame(df2),endo_cancer==0)) #584 rows, 28 columns, making it a matrix

x=df3[,-29] #removing endo_cancer variable
rownames(x)=1:nrow(x)
#B=number of subsamples
#threshold=proabability for stability selection
#p=27, if lgca125 is included, p=26, if lgca125 is excluded

oc_wrapper=function(p, B, threshold){
  
  if (p == 27) {
    x <- x
  } else if (p == 26) {
    x <- x[, -(p + 1)]
  }
  
  nfold = 5
  set.seed(6578)
  s0=which(x[,(p+1)]==0);s1=which(x[,(p+1)]==1) #separating cases and controls for stratified cross-validation
  fold1 = c(sample(s0,length(s0)*0.2, replace=FALSE),sample(s1,length(s1)*0.2, replace=FALSE))
  fold1=unname(fold1)
  fold2 = c(sample(setdiff(s0, fold1), size=length(s0)*0.2, replace=FALSE), sample(setdiff(s1, fold1), size=length(s1)*0.2, replace=FALSE))
  fold3 = c(sample(setdiff(s0, c(fold1, fold2)), size=length(s0)*0.2, replace=FALSE), sample(setdiff(s1, c(fold1, fold2)), size=length(s1)*0.2, replace=FALSE))
  fold4 = c(sample(setdiff(s0, c(fold1, fold2, fold3)), size=length(s0)*0.2, replace=FALSE), sample(setdiff(s1, c(fold1, fold2, fold3)), size=length(s1)*0.2, replace=FALSE))
  fold5 <- c(setdiff(s0, c(fold1, fold2, fold3, fold4)), setdiff(s1, c(fold1, fold2, fold3, fold4)))
  
  
  all.fold <- list(fold1,fold2,fold3,fold4,fold5)
  combinations <- combn(all.fold, 4, simplify = FALSE) # Generate all combinations of 4 vectors
  combined.index <- lapply(combinations, function(comb) do.call(c, comb)) # Perform the combination operation (e.g., concatenate)
  
  true.out1=x[-combined.index[[1]],p+1]
  true.out2=x[-combined.index[[2]],p+1]
  true.out3=x[-combined.index[[3]],p+1]
  true.out4=x[-combined.index[[4]],p+1]
  true.out5=x[-combined.index[[5]],p+1]
  true.out=c(true.out1,true.out2,true.out3,true.out4,true.out5) #true outcome
  
  
  
  pred.prob.rf=vector("numeric")
  pred.prob.lr=vector("numeric")
  pred.prob.lda=vector("numeric")
  pred.prob.svm=vector("numeric")
  pred.prob.com=vector("numeric")
  lasso.pred.prob.com=vector("numeric")
  all.pred.prob.com=vector("numeric")
  selected.variables=list()
  selected.variables.lasso=list()
  
  for (i in 1:nfold){
    set.seed(200*i)
    ####Lasso
    stab.rconcave=stabsel(x=x[combined.index[[i]], -(p+1)], y=x[combined.index[[i]],"ovar_cancer"], fitfun=glmnet.lasso, cutoff=threshold, PFER = 1, 
                          sampling.type = "SS", assumption="r-concave")
    sel.lasso=data.frame(stab.rconcave$selected)$stab.rconcave.selected
    ###Sparse SVM
    sel.var.svm.vec=c()
    for (k in 1:B){
      set.seed(20*k)
      train.svm=sample(combined.index[[i]], length(combined.index[[i]])*0.5, replace=FALSE)
      yy.svm1=x[train.svm,p+1]
      xx.svm1=as.matrix(x[train.svm, -(p+1)])
      svm.mod1=sparseSVM(xx.svm1,yy.svm1, alpha=1, gamma=0.3,lambda=0.15, max.iter=1000, dfmax=p+1) #dfmax=p+1
      d1= as.vector(coef(svm.mod1)[-1])
      sel.var.svm1=which(d1 !=0)
      
      yy.svm2=x[-train.svm,p+1]
      xx.svm2=as.matrix(x[-train.svm, -(p+1)])
      svm.mod2=sparseSVM(xx.svm2, yy.svm2, alpha=1, gamma=0.3,lambda=0.15, max.iter=1000, dfmax=p+1)
      d2= as.vector(coef(svm.mod2) [-1])
      sel.var.svm2=which(d2 !=0)
      
      sel.var.svm=c(sel.var.svm1, sel.var.svm2)
      sel.var.svm.vec=c(sel.var.svm.vec, sel.var.svm)
    }
    
    sel.var.svm.tab=table(sel.var.svm.vec) #frequency of each selected variable with corresponding column number
    sel.svm=as.numeric(names(sel.var.svm.tab[sel.var.svm.tab > threshold*B*2])) # selected variables (columns)
    
    ####Random forest
    sel.var.rf.vec=c()
    for (j in 1:B){
      set.seed(100*j)
      train.rf=sample(combined.index[[i]], length(combined.index[[i]])*0.5, replace=FALSE)
      
      df.rf1=data.frame(x[train.rf, ])
      rf.mod1=Boruta(ovar_cancer ~ ., data = df.rf1, doTrace = 0, ntree = 300, pValue = 0.01)
      sel.var1.rf = which(rf.mod1$finalDecision == "Confirmed")
      
      df.rf2=data.frame(x[-train.rf, ])
      rf.mod2=Boruta(ovar_cancer ~ ., data = df.rf2, doTrace = 0, ntree = 300, pValue = 0.01)
      sel.var2.rf = which(rf.mod2$finalDecision == "Confirmed")
      
      sel.var.rf=c(sel.var1.rf, sel.var2.rf)
      sel.var.rf.vec=c(sel.var.rf.vec, sel.var.rf) # keeping the record of selected variables (column)
    }
    sel.var.rf.tab=table(sel.var.rf.vec) #frequency of each selected variable with corresponding column number
    sel.rf=as.numeric(names(sel.var.rf.tab[sel.var.rf.tab >threshold*B*2])) # selected variables (columns)
    
    sel.lasso.svm=union(sel.lasso, sel.svm) #taking the union of the selected variables by lasso and sparsesvm
    final.set=union(sel.lasso.svm, sel.rf) # final subset of selected variables by STABEL
    
    selected.variables[[i]]=final.set #keeping a record of selected biomarkers
    ########################RF###############################
    sl.model.rf = SuperLearner(Y = x[combined.index[[i]],p+1], X = data.frame(x[combined.index[[i]], final.set]), method="method.NNLS", 
                               family = binomial(), SL.library = "SL.randomForest", cvControl = list(V=10))
    newX=data.frame(x[-combined.index[[i]], final.set])
    sl.pred.rf=predict(sl.model.rf, newX, type="response")
    pred.prob.rf=c(pred.prob.rf, sl.pred.rf$pred)
    
    # ###################logistic################################
    # sl.model.lr = SuperLearner(Y = x[combined.index[[i]],p+1], X = data.frame(x[combined.index[[i]], final.set]), method="method.NNLS", 
    #                            family = binomial(), SL.library = "SL.glm", cvControl = list(V=10))
    # sl.pred.lr=predict(sl.model.lr, newX, type="response")
    # pred.prob.lr=c(pred.prob.lr, sl.pred.lr$pred)
    # 
    # #####################LDA####################################
    # sl.model.lda = SuperLearner(Y = x[combined.index[[i]],p+1], X = data.frame(x[combined.index[[i]], final.set]), method="method.NNLS", 
    #                             family = binomial(), SL.library = "SL.lda", cvControl = list(V=10))
    # sl.pred.lda=predict(sl.model.lda, newX, type="response")
    # pred.prob.lda=c(pred.prob.lda, sl.pred.lda$pred)
    
    ###########################SVM##############################
    sl.model.svm = SuperLearner(Y = x[combined.index[[i]],p+1], X = data.frame(x[combined.index[[i]], final.set]), method="method.NNLS", 
                                family = binomial(), SL.library = "SL.ksvm", cvControl = list(V=10))
    sl.pred.svm=predict(sl.model.svm, newX, type="response")
    pred.prob.svm=c(pred.prob.svm, sl.pred.svm$pred)
    
    ######################Combined#######################
    sl.model.com= SuperLearner(Y = x[combined.index[[i]],p+1], X = data.frame(x[combined.index[[i]], final.set]), method="method.NNLS", 
                               family = binomial(), SL.library = c("SL.randomForest","SL.glm", "SL.ksvm", "SL.lda"), cvControl = list(V=10))
    sl.pred.com=predict(sl.model.com, newX, type="response")
    pred.prob.com=c(pred.prob.com, sl.pred.com$pred)
    
    ##########################LASSO with Superlearner######################
    set.seed(6578)
    cv.fit=  cv.glmnet(as.matrix(x[combined.index[[i]], -(p+1)]), as.factor(x[combined.index[[i]], p+1]), family="binomial",alpha=1)
    bestlam=cv.fit$lambda.min
    final.set.lasso=which(coef(cv.fit, s=bestlam)[-1] !=0)
    selected.variables.lasso[[i]]=final.set.lasso
    newX.lasso = data.frame(x[-combined.index[[i]], final.set.lasso])
    lasso.model.com= SuperLearner(Y = x[combined.index[[i]],p+1], X = data.frame(x[combined.index[[i]], final.set.lasso]), method="method.NNLS", 
                                  family = binomial(), SL.library = c("SL.randomForest","SL.glm", "SL.ksvm", "SL.lda"), cvControl = list(V=10))
    lasso.pred.com=predict(lasso.model.com, newX.lasso, type="response")
    lasso.pred.prob.com=c(lasso.pred.prob.com, lasso.pred.com$pred)
    
    #########################All variables with Superlearner######################
    newX.all = data.frame(x[-combined.index[[i]], -(p+1)])
    
    all.model.com = SuperLearner(Y = x[combined.index[[i]], p+1], X = data.frame(x[combined.index[[i]], -(p+1)]), method="method.NNLS", 
                                 family = binomial(), SL.library = c("SL.randomForest","SL.glm", "SL.ksvm", "SL.lda"), cvControl = list(V=6))
    all.pred.com=predict(all.model.com, newX.all, type="response")
    all.pred.prob.com=c(all.pred.prob.com, all.pred.com$pred)
    
  }
  
  ############random forest###################
  ROC.sl.rf <- roc(as.vector(true.out), as.vector(pred.prob.rf))
  auc.rf <- auc(ROC.sl.rf)
  ci.rf=ci.auc(ROC.sl.rf, conf.level=0.95, method="bootstrap", boot.n=1000)
  cord.rf=coords(ROC.sl.rf, x="best", input="threshold", best.method="youden")
  pred.out.rf=ifelse(pred.prob.rf>cord.rf$threshold,1,0)
  cutoffs <- c(0.985, 0.95) #for specificity
  myData.rf <- with(ROC.sl.rf, data.frame(specificities, sensitivities, thresholds))
  sens.rf.spec=lapply(cutoffs, function(cutoff) myData.rf$sensitivities[which.min(abs(myData.rf$specificities-cutoff))])
  t.rf=table(true.out,pred.out.rf)
  sens.rf=t.rf[2,2]/(t.rf[2,2]+t.rf[2,1])
  spec.rf=t.rf[1,1]/(t.rf[1,2]+t.rf[1,1])
  pred.acc.rf <- 1-mean(pred.out.rf != true.out)
  
  #########################SVM####################################
  ROC.sl.svm <- roc(as.vector(true.out), as.vector(pred.prob.svm))
  auc.svm <- auc(ROC.sl.svm)
  ci.svm=ci.auc(ROC.sl.svm, conf.level=0.95, method="bootstrap", boot.n=1000)
  cord.svm=coords(ROC.sl.svm, x="best", input="threshold", best.method="youden")
  pred.out.svm=ifelse(pred.prob.svm>cord.svm$threshold,1,0)
  myData.svm <- with(ROC.sl.svm, data.frame(specificities, sensitivities, thresholds))
  sens.svm.spec=lapply(cutoffs, function(cutoff) myData.svm$sensitivities[which.min(abs(myData.svm$specificities-cutoff))])
  t.svm=table(true.out,pred.out.svm)
  sens.svm=t.svm[2,2]/(t.svm[2,2]+t.svm[2,1])
  spec.svm=t.svm[1,1]/(t.svm[1,2]+t.svm[1,1])
  pred.acc.svm <- 1-mean(pred.out.svm != true.out)
  
  ######################combined (Superlearner)-STABEL####################################
  ROC.sl.com <- roc(as.vector(true.out), as.vector(pred.prob.com))
  auc.com <- auc(ROC.sl.com)
  ci.com=ci.auc(ROC.sl.com, conf.level=0.95, method="bootstrap", boot.n=1000)
  cord.com=coords(ROC.sl.com, x="best", input="threshold", best.method="youden")
  pred.out.com=ifelse(pred.prob.com>cord.com$threshold,1,0)
  myData.com <- with(ROC.sl.com, data.frame(specificities, sensitivities, thresholds))
  sens.com.spec=lapply(cutoffs, function(cutoff) myData.com$sensitivities[which.min(abs(myData.com$specificities-cutoff))])
  t.com=table(true.out,pred.out.com)
  sens.com=t.com[2,2]/(t.com[2,2]+t.com[2,1])
  spec.com=t.com[1,1]/(t.com[1,2]+t.com[1,1])
  pred.acc.com <- 1-mean(pred.out.com != true.out)
  
  #################LASSO with superlearner#########################
  ROC.lasso.com <- roc(as.vector(true.out), as.vector(lasso.pred.prob.com))
  lasso.auc.com <- auc(ROC.lasso.com)
  lasso.ci.com=ci.auc(ROC.lasso.com, conf.level=0.95, method="bootstrap", boot.n=1000)
  lasso.cord.com=coords(ROC.lasso.com, x="best", input="threshold", best.method="youden")
  lasso.pred.out.com=ifelse(lasso.pred.prob.com>lasso.cord.com$threshold,1,0)
  myData.lasso <- with(ROC.lasso.com, data.frame(specificities, sensitivities, thresholds))
  sens.lasso.spec=lapply(cutoffs, function(cutoff) myData.lasso$sensitivities[which.min(abs(myData.lasso$specificities-cutoff))])
  lasso.t.com=table(true.out,lasso.pred.out.com)
  lasso.sens.com=lasso.t.com[2,2]/(lasso.t.com[2,2]+lasso.t.com[2,1])
  lasso.spec.com=lasso.t.com[1,1]/(lasso.t.com[1,2]+lasso.t.com[1,1])
  lasso.pred.acc.com <- 1-mean(lasso.pred.out.com != true.out)
  
  ######################All variables with superlearner#############
  ROC.all.com <- roc(as.vector(true.out), as.vector(all.pred.prob.com))
  all.auc.com <- auc(ROC.all.com)
  all.ci.com=ci.auc(ROC.all.com, conf.level=0.95, method="bootstrap", boot.n=1000)
  all.cord.com=coords(ROC.all.com, x="best", input="threshold", best.method="youden")
  all.pred.out.com=ifelse(all.pred.prob.com>all.cord.com$threshold,1,0)
  myData.all <- with(ROC.all.com, data.frame(specificities, sensitivities, thresholds))
  sens.all.spec=lapply(cutoffs, function(cutoff) myData.all$sensitivities[which.min(abs(myData.all$specificities-cutoff))])
  all.t.com=table(true.out, all.pred.out.com)
  all.sens.com=all.t.com[2,2]/(all.t.com[2,2]+all.t.com[2,1])
  all.spec.com=all.t.com[1,1]/(all.t.com[1,2]+all.t.com[1,1])
  all.pred.acc.com <- 1-mean(all.pred.out.com != true.out)
  
  ###################################################################################
  ##################################2-way interactions###############################
  ###################################################################################
  
  outcome=x[ ,(p+1)]
  main_effect=x[ ,-(p+1)]
  variable_combinations <- combn(colnames(main_effect), 2)
  
  # Create a data frame to store the interactions
  inter_effect <- matrix(NA, nrow=nrow(main_effect), dim(variable_combinations)[2])
  
  # Loop through the variable combinations and calculate interactions
  for (i in 1:ncol(variable_combinations)) {
    int1 <- variable_combinations[1, i]
    int2 <- variable_combinations[2, i]
    interaction <- main_effect[, int1] * main_effect[, int2]
    inter_effect[, i]=interaction
  }
  
  colnames(inter_effect) <- paste0("interaction", 1:dim(variable_combinations)[2])
  
  data=cbind(main_effect, inter_effect,outcome)
  
  pval_mat=matrix(NA, nrow= (dim(data)[2]-1), nfold)
  for (j in 1:nfold){
    #pval_mat=matrix(NA, nrow= (dim(data)[2]-1), 1)
    pval_vec=c()
    for(i in 1:(dim(data)[2]-1)){
      mylogit <- glm(data[combined.index[[j]],dim(data)[2]] ~ data[combined.index[[j]],i], data=data.frame(data), family="binomial")
      pval=summary(mylogit)$coefficients[, "Pr(>|z|)"]
      pval_vec[i]=pval[2]
    }
    pval_mat[,j]=pval_vec
  }
  colnames(pval_mat) <- paste0("pvalue", 1:nfold) #378*5, variables=378, fold=5
  
  
  # Function to select the 92=n/log(n) lowest values for each column with row numbers
  select_lowest_values <- function(matrix) {
    selected_values <- lapply(1:ncol(matrix), function(col_index) {
      col <- matrix[, col_index]
      sorted_indices <- order(col)
      lowest_indices <- sorted_indices[1:92]
      lowest_values <- col[lowest_indices]
      data.frame(Row = lowest_indices, Value = lowest_values)
    })
    return(selected_values)
  }
  
  # Use the function to select the lowest values for each column
  lowest_values <- select_lowest_values(pval_mat)
  
  # cov_set1=unname(unlist(lowest_values[[1]][1])) #extracting the covariates from previous matrix
  # cov_set2=unname(unlist(lowest_values[[2]][1]))
  # cov_set3=unname(unlist(lowest_values[[3]][1]))
  # cov_set4=unname(unlist(lowest_values[[4]][1]))
  # cov_set5=unname(unlist(lowest_values[[5]][1]))
  
  int.pred.prob.com=vector("numeric")
  for (m in 1:nfold){
    int.model.com = SuperLearner(Y = data[combined.index[[m]],dim(data)[2]], X = data.frame(data[combined.index[[m]], unname(unlist(lowest_values[[m]][1]))]), method="method.NNLS", 
                                family = binomial(), SL.library = c("SL.randomForest","SL.glm", "SL.ksvm", "SL.lda"), cvControl = list(V=10))
    newX.int=data.frame(data[-combined.index[[m]], unname(unlist(lowest_values[[m]][1]))])
    int.pred.com=predict(int.model.com, newX.int, type="response")
    int.pred.prob.com=c(int.pred.prob.com, int.pred.com$pred)
  }
  ROC.int.com <- roc(as.vector(true.out), as.vector(int.pred.prob.com))
  int.auc.com <- auc(ROC.int.com)
  int.ci.com=ci.auc(ROC.int.com, conf.level=0.95, method="bootstrap", boot.n=1000)
  int.cord.com=coords(ROC.int.com, x="best", input="threshold", best.method="youden")
  int.pred.out.com=ifelse(int.pred.prob.com>int.cord.com$threshold,1,0)
  myData.int <- with(ROC.int.com, data.frame(specificities, sensitivities, thresholds))
  sens.int.spec=lapply(cutoffs, function(cutoff) myData.int$sensitivities[which.min(abs(myData.int$specificities-cutoff))])
  int.t.com=table(true.out,int.pred.out.com)
  int.sens.com=int.t.com[2,2]/(int.t.com[2,2]+int.t.com[2,1])
  int.spec.com=int.t.com[1,1]/(int.t.com[1,2]+int.t.com[1,1])
  int.pred.acc.com <- 1-mean(int.pred.out.com != true.out)
  
  list(selected.variables=selected.variables, selected.variables.lasso=selected.variables.lasso,
    ROC.sl.com=ROC.sl.com,ROC.sl.rf=ROC.sl.rf,ROC.sl.svm=ROC.sl.svm,ROC.lasso.com=ROC.lasso.com,ROC.all.com=ROC.all.com,ROC.int.com=ROC.int.com,
    
    Accuracy_COM = pred.acc.com, AUC_COM = auc.com, CI_COM = paste0(ci.com[1], "-", ci.com[3]), Sensitivity_COM = sens.com, Specificity_COM = spec.com,
    SENS_95_SPEC_COM=sens.com.spec[[2]], SENS_98.5_SPEC_COM=sens.com.spec[[1]],
    
    Accuracy_RF = pred.acc.rf, AUC_RF = auc.rf, CI_RF = paste0(ci.rf[1], "-", ci.rf[3]), Sensitivity_RF = sens.rf, Specificity_RF = spec.rf,
    SENS_95_SPEC_RF=sens.rf.spec[[2]], SENS_98.5_SPEC_RF=sens.rf.spec[[1]],
    
    Accuracy_SVM = pred.acc.svm, AUC_SVM = auc.svm, CI_SVM = paste0(ci.svm[1], "-", ci.svm[3]), Sensitivity_SVM = sens.svm, Specificity_SVM = spec.svm, 
    SENS_95_SPEC_SVM=sens.svm.spec[[2]], SENS_98.5_SPEC_SVM=sens.svm.spec[[1]],
    
    Accuracy_LASSO = lasso.pred.acc.com, AUC_LASSO = lasso.auc.com, CI_LASSO = paste0(lasso.ci.com[1], "-", lasso.ci.com[3]), 
    Sensitivity_LASSO = lasso.sens.com, Specificity_LASSO = lasso.spec.com, SENS_95_SPEC_LASSO=sens.lasso.spec[[2]], SENS_98.5_SPEC_LASSO=sens.lasso.spec[[1]], 
    
    Accuracy_ALL = all.pred.acc.com, AUC_ALL = all.auc.com, CI_ALL = paste0(all.ci.com[1], "-", all.ci.com[3]), Sensitivity_ALL = all.sens.com, 
    Specificity_ALL = all.spec.com, SENS_95_SPEC_ALL=sens.all.spec[[2]], SENS_98.5_SPEC_ALL=sens.all.spec[[1]],
    
    Accuracy_INT = int.pred.acc.com, AUC_INT = int.auc.com, CI_INT = paste0(int.ci.com[1], "-", int.ci.com[3]), Sensitivity_INT = int.sens.com, 
    Specificity_INT = int.spec.com, SENS_95_SPEC_INT=sens.int.spec[[2]], SENS_98.5_SPEC_INT=sens.int.spec[[1]])

}

result=oc_wrapper(p=27, B=50, threshold=0.75) ####Wrapper
#p=26 excludes Ca125 biomarker, p=27 includes all biomarkers
time_taken=Sys.time()-t1

result$selected.variables #selected variables from STABEL
result$selected.variables.lasso #selected variables from LASSO

# Remove selected items from the list
filtered_result <- result[!names(result) %in% c("selected.variables", "selected.variables.lasso",
                                                "ROC.sl.com", "ROC.sl.rf", "ROC.sl.svm", "ROC.lasso.com", "ROC.all.com", "ROC.int.com")]

# Convert to a dataframe
result_df <- as.data.frame(matrix(unlist(filtered_result), ncol=7, byrow=TRUE))
colnames(result_df) <- c("Accuracy", "AUC", "CI", "Sensitivity", "Specificity", "SENS_95_SPEC", "SENS_98.5_SPEC")
rownames(result_df) <- c(  "STABEL", "STABEL-RF", "STABEL-SVM", "LASSO-EL", "AV-EL", "INT-EL")



plot(result$ROC.sl.com, legacy.axes=TRUE, col = "black", xlab="1-Specificity", lty=1,type="l", lwd=1,cex.lab=1.5) 
lines(result$ROC.sl.rf, col = "chartreuse2",lty=2, type="l", lwd=1)
lines(result$ROC.sl.svm, col = "blue",lty=3, type="l", lwd=1)
lines(result$ROC.lasso.com, col = "red",lty=4, type="l", lwd=1)
lines(result$ROC.all.com, col = "darkorchid",lty=5, type="l", lwd=1)


legend(0.60,0.25, legend=c("STABEL, AUC=0.851(0.823-0.881)", 
                           "STABEL-RF, AUC=0.843(0.815-0.874)",
                           "STABEL-SVM, AUC=0.839(0.806-0.872)",
                           "LASSO-EL, AUC=0.841(0.810-0.871)",
                           "AV-EL, AUC=0.843(0.811-0.871)"), 
       col=c("black", "chartreuse2", "blue", "red", "darkorchid"), 
       lty=c(1,2,3,4,5,6), cex=1, lwd = 2)
