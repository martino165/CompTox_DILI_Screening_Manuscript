#DILI DATA FOR MODEL PREPARATION

data_prep <- function(dat,
                      dili_positive = c("Most"),
                      dili_negative = c("None", "Less"),
                      col_sets = c("dose", "physchem"),
                      required_non_missing = c("dili", "dose"),
                      allowed_fraction_missing = 0.2,
                      fill_missing = T,
                      remove_correlated_cutoff = 0.8
){
  
  require(data.table)
  require(caret)
  require(car)
  
  #subset the data to just the dili annotations that match either the
  #positive or negative designation
  dat <- dat[dili %in% c(dili_positive, dili_negative)]
  
  dat[ , dili := factor(dili, levels = c(dili_positive, dili_negative))]
  dat[ , dili := droplevels(dili)]
  dat[ , Y := as.factor(ifelse(dili %in% dili_negative,
                               "dili_negative", "dili_positive"))]
  
  
  #subsets dat to just the fields that must have no missing values and then
  #pulls the row indexes for all rows that have complete cases (no NA values)
  row_subset <- which(complete.cases(dat[ , (required_non_missing), with = F]))
  
  col_set <- list("dose" = "dose",
                  "cmax" = "cmax",
                  "physchem" = c("mw","fsp3","clogp", "clogd", "psa", "ion", "eccs"),
                  "bddcs" = c("bddcs"),
                  "insilico_pk" = c("cFuHPPB", "cFuBrain"),
                  "invitro_pk" = c("fu_liver", "fu_bsa", "kp", "kpuu"),
                  "invitro_pk_exp" = c("fu_liver_exp", "fu_bsa_exp", "kp", "kpuu_exp"), 
                  "invitro_safety" = c("bsep", "thle", "hepg2_glu", "hepg2_gal",
                                       "hepg2_72h", "mrp2", "mrp3", "mrp4"))
  
  cols <- unlist(col_set[col_sets])
  
  cids <- dat[row_subset, cid]
  Y <- dat[row_subset, Y]
  X <- dat[row_subset, cols, with = F]
  log_cols <- c("dose", "cmax", "kp", "kpuu", "Cmax","Cavg", "CLp", "CLr",
                "Rbp", "CLh", "Chmax", "mw", 
                "fu_liver", "fu_bsa", "fu_plasma", "Fup", "cFuHPPB", "cFuBrain") #previous logit transformations (removing for simplicity)
  
  log_cols <- log_cols[which(log_cols %in% colnames(X))]
  if(length(log_cols)>0) X[ , (log_cols) := lapply(.SD, log10), .SDcols = log_cols]
  
  # logit_cols <- c("fu_liver", "fu_bsa", "fu_plasma", "Fup", "cFuHPPB", "cFuBrain")
  # logit_cols <- logit_cols[which(logit_cols %in% colnames(X))]
  # if(length(logit_cols)>0) X[ , (logit_cols) := lapply(.SD, logit), .SDcols = logit_cols]
  # 
  DATA <- as.data.frame(cbind(Y,X)) #gets repeated once missing data is imputed
  rownames(DATA) <- cids
  
  num_cols <- names(which(sapply(DATA, typeof)=="double"))
  
  if(length(num_cols) > 1) {
    cormat <- cor(as.matrix(DATA[,num_cols]), use = "complete.obs")
    
    remove.correlated.var.index <- caret::findCorrelation(cormat,
                                                          cutoff = remove_correlated_cutoff)
    
    if(length(remove.correlated.var.index) >  0){
      remove.cor.col <- colnames(cormat)[remove.correlated.var.index]
      DATA[ , remove.cor.col] <- list(NULL)
    }
    
  } #only run for data with continuous values
  
  FRACTION_MISSING <- allowed_fraction_missing
  keep_col <- names(which(apply(DATA, 2, function(x) length(which(is.na(x))))/nrow(X) <= FRACTION_MISSING))
  DATA <- DATA[ , keep_col]
  
  #if only Y value left remaining then just return NULL
  if(length(keep_col)==1){
    return(NULL)
  }
  
  if(!nrow(DATA)==length(which(complete.cases(DATA))) & fill_missing == T){
    DATA1 <- preProcess(DATA, method = "medianImpute")
    DATA <-  predict(DATA1, DATA)
  }
  
  set.seed(123456)
  trainIndex <- createDataPartition(DATA$Y, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  TRAIN <- DATA[trainIndex,]
  TEST  <- DATA[-trainIndex,]
  
  return(list(
    DATA = DATA,
    TRAIN = TRAIN,
    TEST = TEST
  ))
  
}