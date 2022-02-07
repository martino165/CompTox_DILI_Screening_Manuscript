#dili model_suite


dili_model_suite <- function(DATA,
                             models2run = c("xgbTree","gbm","cforest","ranger"),
                             sampling = c("up", "down"),
                             cvnumber = 5,
                             cvrepeat = 5) {
  
  require(xgboost)
  require(gbm)
  require(ranger)
  require(caret)
  require(party)

  twoClassSummary2 <- function (data, lev = NULL, model = NULL) {
    lvls <- levels(data$obs)
    if (length(lvls) > 2) 
      stop(paste("Your outcome has", length(lvls), "levels. The twoClassSummary() function isn't appropriate."))
    caret:::requireNamespaceQuietStop("ModelMetrics")
    if (!all(levels(data[, "pred"]) == lvls)) 
      stop("levels of observed and predicted data do not match")
    rocAUC <- ModelMetrics::auc(ifelse(data$obs == lev[2], 0, 
                                       1), data[, lvls[1]])
    sens <- sensitivity(data[, "pred"], data[, "obs"], lev[1])
    spec <- specificity(data[, "pred"], data[, "obs"], lev[2])
    ba <- (sens+spec)/2
    
    cm <- confusionMatrix(data = data[, "pred"],
                          reference = data[, "obs"])
    cm$overall
    
    out <- c(rocAUC, sens, spec, ba, cm$overall[c(1,2,7)])
    names(out) <- c("ROC", "Sens", "Spec", "BA", "Accuracy", "Kappa", "McPval")
    out
  }
  
  ctrl <- trainControl(method = "repeatedcv",
                       number = cvnumber,
                       repeats = cvrepeat,
                       summaryFunction = twoClassSummary2,
                       classProbs = TRUE,
                       savePredictions = T,
                       verboseIter = FALSE)
  
  if(!"xgbTree" %in% models2run){
    xgb_up_model <- NULL
    xgb_dn_model <- NULL
  } else {
    print("Running xgbTree")
    tGrid <- expand.grid(nrounds = c(5, 10, 25, 50),
                         max_depth = 1:3,
                         eta = 0.1, 
                         gamma = 0,
                         colsample_bytree = 0.7,
                         min_child_weight = 1,
                         subsample = 1)
    
    
    if("up" %in% sampling){
      ctrl$sampling <- "up"
      xgb_up_model <- train(Y ~ .,
                            data = DATA,
                            method = "xgbTree",
                            metric = "ROC",
                            trControl = ctrl,
                            tuneGrid = tGrid)
    } else {
      xgb_up_model <- NULL
    }
    
    if("down" %in% sampling){
      ctrl$sampling <- "down"
      xgb_dn_model <- train(Y ~ .,
                            data = DATA,
                            method = "xgbTree",
                            metric = "ROC",
                            trControl = ctrl,
                            tuneGrid = tGrid)
    } else {
      xgb_dn_model <- NULL
    }
  }
  
  
  if(!"gbm" %in% models2run){
    gbm_up_model <- NULL
    gbm_dn_model <- NULL
  } else {
    
    print("Running gbm")
    
    tGrid <- expand.grid(interaction.depth = 1:3,
                         n.trees = c(5, 10, 25, 50),
                         shrinkage = 0.1,
                         n.minobsinnode = 10)
    
    if("up" %in% sampling){
      ctrl$sampling <- "up"
      gbm_up_model <- train(Y ~ .,
                            data = DATA,
                            method='gbm',
                            trControl=ctrl,
                            metric = "ROC",
                            preProc = c("center", "scale"),
                            tuneGrid = tGrid
      )
    } else {
      gbm_up_model <- NULL
    }
    
    if("down" %in% sampling){
      ctrl$sampling <- "down"
      gbm_dn_model <- train(Y ~ .,
                            data = DATA,
                            method='gbm',
                            trControl=ctrl,
                            metric = "ROC",
                            preProc = c("center", "scale"),
                            tuneGrid = tGrid
      )
    } else {
      gbm_dn_model <- NULL
    }
  } 
  
  if(!"cforest" %in% models2run){
    cfr_up_model <- NULL
    cfr_dn_model <- NULL
  } else {
    
    print("Running cforest")
    
    if("up" %in% sampling){
      ctrl$sampling <- "up"
      cfr_up_model <- train(Y ~ .,
                            data = DATA,
                            method='cforest',
                            controls = cforest_unbiased(ntree = 1000),
                            trControl= ctrl,
                            metric = "ROC"
      )
    } else {
      cfr_up_model <- NULL
    }
    if("down" %in% sampling){
      ctrl$sampling <- "down"
      cfr_dn_model <- train(Y ~ .,
                            data = DATA,
                            method='cforest',
                            controls = cforest_unbiased(ntree = 1000),
                            trControl= ctrl,
                            metric = "ROC"
      )
    } else {
      cfr_dn_model <- NULL
    }
  }
  
  
  
  if(!"ranger" %in% models2run){
    rrf_up_model <- NULL
    rrf_dn_model <- NULL
  } else {
    
    print("Running ranger")  
    
    tgrid <- expand.grid(
      .mtry = unique(c(floor(sqrt(ncol(DATA)-1)), ceiling(sqrt(ncol(DATA)-1)))),
      .splitrule = "gini",
      .min.node.size = c(1)
    )
    if("up" %in% sampling){
      ctrl$sampling <- "up" 
      rrf_up_model <- train(Y ~ .,
                            data = DATA,
                            method = "ranger",
                            verbose = FALSE,
                            metric = "ROC",
                            maximize = TRUE,
                            trControl = ctrl,
                            tuneGrid = tgrid,
                            num.trees = 1000,
                            importance = "impurity",
                            respect.unordered.factors = F)
    } else {
      rrf_up_model <- NULL
    }
    if("down" %in% sampling){
      ctrl$sampling <- "down" 
      rrf_dn_model <- train(Y ~ .,
                            data = DATA,
                            method = "ranger",
                            verbose = FALSE,
                            metric = "ROC",
                            maximize = TRUE,
                            trControl = ctrl,
                            tuneGrid = tgrid,
                            num.trees = 1000,
                            importance = "impurity",
                            respect.unordered.factors = F)
    } else {
      rrf_dn_model <- NULL
    }
  }
  
  return(list(
    xgb_up_model = xgb_up_model,
    xgb_dn_model = xgb_dn_model,
    gbm_up_model = gbm_up_model,
    gbm_dn_model = gbm_dn_model,
    cfr_up_model = cfr_up_model,
    cfr_dn_model = cfr_dn_model,
    rrf_up_model = rrf_up_model,
    rrf_dn_model = rrf_dn_model
  ))
  
}