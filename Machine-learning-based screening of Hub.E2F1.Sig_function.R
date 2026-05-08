#' 循环进行lasso分析
#' 
#' @param genes(characters) 选择要筛选的基因
#' @param group(factor or Surv) 选择分组信息，诊断模型是二分类变量，预后模型是生存信息（包括time和event）event需要转变为0和1
#' @param data(data.frame) 包含基因和分组信息的数据
#' @param lambda_choose(characters) 'lambda.1se', 'lambda.min'。lambda.min下的模型交叉验证误差最小。lambda.1se给出了最正则化的模型，使得交叉验证误差在最小值的一个标准误差内。
#' @param nfold(numeric) K-fold交叉验证，数据量少可以选5
#' @param iter.time(numeric) 循环次数
#' 
#' @return (data.frame) lasso筛选的基因，基因数和种子数
#' 
#' @export
lasso_iter <- function(genes, group, data, lambda_choose = c('lambda.1se', 'lambda.min'), nfold=10, iter.times=100){
  #cv.glmnet输入变量：x, y
  mat_x <- as.matrix(data[, genes])
  y <- data[, group]
  lasso_methods <- "binomial"
  #group如果为两个：预后lasso
  if(length(group) > 1){
    if(length(group) == 2){
      #判断group变量
      len <- apply(y, 2, function(x) length(unique(x)))
      event_col <- which(len == 2)
      time_col <- which(len != 2)
      if(identical(event_col, integer(0))){
        stop("No event variables were found\n")
      } else{
        cat(paste0("use ", colnames(y)[event_col], " as event variable\n"),
            paste0("use ", colnames(y)[time_col], " as time variable\n"))
        y <- data.matrix(Surv(time=as.double(y[, time_col]), event=as.double(y[, event_col])))
        lasso_methods <- "cox"
      }
    }else{
      stop("The length of group must be 1 or 2\n")
    }
  }
  lambda_choose <- match.arg(lambda_choose)
  #cv.glmnet循环iter.times次，选择cv.glmnet结果的最优解
  lasso_fea_list <- lapply(seq_len(iter.times), function(x){ # 
    set.seed(x)
    cvfit <- cv.glmnet(mat_x, y, family= lasso_methods, nfolds = nfold)
    # 取出最优lambda
    fea <- rownames(coef(cvfit, s = lambda_choose))[coef(cvfit, s = lambda_choose)[, 1]!= 0]
    if(is.element("(Intercept)", fea)) {
      fea <- fea[-1] # 去掉截距项并排序
    } else {
      fea <- fea
    }
    lasso_res <- data.frame(sed = x,
                            len = length(fea),
                            genes = paste0(sort(fea), collapse = " | "),
                            stringsAsFactors = F)
  })
  dat_res <- do.call(rbind, lasso_fea_list)
  dat_res$marker <- factor(dat_res$genes,labels = LETTERS[1:length(unique(dat_res$genes))])
  dat_res$marker <- paste(dat_res$len, "genes", dat_res$marker, sep = " ")
  return(dat_res)
}

machine_res <- function(df){
  ls_res <- list()
  cl <- makePSOCKcluster(10)
  clusterEvalQ(cl, .libPaths("/work/reference/Rlibrary/4.2"))
  registerDoParallel(cl)
  set.seed(1234)

  #LQV
  control_lqv <- trainControl(method="repeatedcv", number=10, repeats=3)
  ## train the model
  lqv_res <- train(response~., data=df, method="lvq", preProcess="scale", trControl=control_lqv)
  importance <- varImp(lqv_res, scale=FALSE)
  
  #bagged trees
  # define the control using a random forest selection function
  control_bag <- rfeControl(functions=treebagFuncs, method="cv", number=30)
  ## run the RFE algorithm
  treebag_res <- rfe(df[,2:ncol(df)], df[,1], sizes=c(2:ncol(df)), rfeControl=control_bag)
  
  #boruta
  boruta_res <- Boruta(response~., data = df, doTrace = 3)
  
  #random forest
  ## define the control using a random forest selection function
  control_rfe <- rfeControl(functions=rfFuncs, method="cv", number=30)
  ## run the RFE algorithm
  rfe_res <- rfe(df[,2:ncol(df)], df[,1], sizes=c(2:ncol(df)), rfeControl=control_rfe)
  
  #Bayesian
  control_bayesian <- rfeControl(functions=nbFuncs, method="cv", number=30)
  ## run the RFE algorithm
  bayesian_res <- rfe(df[,2:ncol(df)], df[,1], sizes=c(2:ncol(df)), rfeControl=control_bayesian)
  
  stopCluster(cl)
  
  #hubgenes
  lqv_hubgenes <- rownames(importance$importance)[which(importance$importance$R > cutoff_line)]
  treebag_hubgenes <- predictors(treebag_res)
  boruta_hubgenes <- TentativeRoughFix(boruta_res)
  boruta_hubgenes <- getSelectedAttributes(boruta_hubgenes, withTentative = F)
  ref_hubgenes <- predictors(rfe_res)
  bayesian_hubgenes <- predictors(bayesian_res)
  ls_res[["res"]] <- list(lqv_res = importance, treebag_res = treebag_res,
                        boruta_res = boruta_res, rfe_res = rfe_res,
                        bayesian_res = bayesian_res)
  ls_res[["hubgenes"]] <- list(lqv = lqv_hubgenes,
                             treebag = treebag_hubgenes,
                             boruta = boruta_hubgenes,
                             ref = ref_hubgenes,
                             bayesian = bayesian_hubgenes)
  return(ls_res)
}












