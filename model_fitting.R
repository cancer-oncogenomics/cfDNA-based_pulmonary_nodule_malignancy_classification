# the module takes in data and parameters to fit models

### load packages
library(dplyr)
library(tidyr)
library(stringr)
library(pROC)
library(h2o)

fit_model <- function(feature_data = feature_data,
                      feature_type = "FSP,CNV", # feature_type is a list of feature types separated by commas
                      sample_info = sample_info, # sample_contains train/valid group information
                      l2_stacking = TRUE, # l2_stacking controls 2nd layer stacking
                      output_path = "./output_FSP/"){
  
  ### initiate h2o
  # h2o.init(port = 18181 + sample(1:1000,1), nthreads = 16, max_mem_size = "256G")
  
  ### create output directory
  dir.create(output_path, recursive = T)
  
  ### label feature data; extract train/valid frames based on sample_info
  feature_data_labeled <- feature_data %>%
    inner_join(select(sample_info, ID, train, fold_column), by = c("ID")) %>%
    mutate(Train_Group = factor(Train_Group, levels = c("BLN","Cancer")))
    
  train_frame <- feature_data_labeled %>%
    filter(train == "training")
  
  valid_frame <- feature_data_labeled %>%
    filter(train == "validation")
  
  ### create score frames for export
  train_score <- select(train_frame, ID, Train_Group) 
  valid_score <- select(valid_frame, ID, Train_Group)
  
  ### start fitting base models by specified feature type
  model_list <- list()
  predictor_list <- list()
  leaderboard_top <- list()
  cols <- colnames(feature_data)
  leaderboard_top_output <- data.frame()
  
  for(feature_type_tmp in unlist(str_split(feature_type,","))){
    
    predictor_list[[feature_type_tmp]] <- cols[grepl(paste0(feature_type_tmp,"_"), cols)]
  
    model_list[[feature_type_tmp]] <- h2o.automl(x = predictor_list[[feature_type_tmp]],
                                             y = "Train_Group",
                                             
                                             project_name = feature_type_tmp,
                                             
                                             # nfolds = 5,
                                             fold_column = "fold_column", # pre-specified fold column to allow for reproducibility
                                             balance_classes = T,
                                             
                                             keep_cross_validation_predictions = T,
                                             keep_cross_validation_fold_assignment = T,
                                             
                                             stopping_metric = "AUC", 
                                             # stopping_tolerance = 0.01, 
                                             max_runtime_secs_per_model=1200, 
                                             
                                             include_algos = c("XgBoost","GLM","DRF","DeepLearning"),
                                             # exclude_algos = c("DRF"),
                                             
                                             training_frame = as.h2o(train_frame), 

                                             max_models = 50, 
                                             seed = 99)
    
    ##### get leaderboard
    leaderboard_top[[feature_type_tmp]] <- model_list[[feature_type_tmp]]@leaderboard %>% as.data.frame() %>%
      mutate(Train_AUC = NA, Valid_AUC = NA, 
             Train_cutoff95 = NA, Train_sens95 = NA, Train_spec95 = NA, Valid_sens95 = NA, Valid_spec95 = NA,
             Train_cutoff98 = NA, Train_sens98 = NA, Train_spec98 = NA, Valid_sens98 = NA, Valid_spec98 = NA) %>%
      head(5)
    
    ##### save top5 base models by AUC in training frame (5XCV)
    for(i in 1:5){
      
      model_i <- h2o.getModel(leaderboard_top[[feature_type_tmp]]$model_id[i])
      varimp_i <- h2o.varimp(model_i) %>% as.data.frame()
      h2o.saveModel(model_i, path = paste0(output_path,"/Base_model_",feature_type_tmp,"_",i), 
                    export_cross_validation_predictions = T, force = T)
      write.csv(varimp_i, paste0(output_path,"/Base_model_",feature_type_tmp,"_",i,
                                 "/varimp_",feature_type_tmp,"_",i,".csv"), row.names = F)
      
      ####### save performance metrics
      leaderboard_top[[feature_type_tmp]]$Train_AUC[i] = h2o.auc(h2o.performance(model_i, xval = T))

      if(nrow(valid_frame) > 0){
        leaderboard_top[[feature_type_tmp]]$Valid_AUC[i] = h2o.auc(h2o.performance(model_i, newdata = as.h2o(valid_frame)))
      }
      
      ####### get predict score
      model_name <- paste0(feature_type_tmp,"_",i)
      train_score[, model_name] = as.numeric(as.data.frame(h2o.getFrame(model_i@model$cross_validation_holdout_predictions_frame_id$name))[,"Cancer"])
      train_score$fold_column = train_frame$fold_column
      
      ####### get roc metrics of train frame
      roc_tmp <- roc(train_score$Train_Group,train_score[,model_name],plot = F,levels = c("BLN", "Cancer"))
      coords95_tmp <- coords(roc = roc_tmp, x = 'all', transpose = F, as.matrix = T) %>%
        as.data.frame() %>%
        filter(specificity >= 0.95) %>%
        filter(specificity <= min(specificity)) %>%
        arrange(desc(sensitivity))
      leaderboard_top[[feature_type_tmp]]$Train_cutoff95[i] <- coords95_tmp$threshold[1]
      leaderboard_top[[feature_type_tmp]]$Train_sens95[i] <- coords95_tmp$sensitivity[1]
      leaderboard_top[[feature_type_tmp]]$Train_spec95[i] <- coords95_tmp$specificity[1]
      
      coords98_tmp <- coords(roc = roc_tmp, x = 'all', transpose = F, as.matrix = T) %>%
        as.data.frame() %>%
        filter(specificity >= 0.98) %>%
        filter(specificity <= min(specificity)) %>%
        arrange(desc(sensitivity))
      leaderboard_top[[feature_type_tmp]]$Train_cutoff98[i] <- coords98_tmp$threshold[1]
      leaderboard_top[[feature_type_tmp]]$Train_sens98[i] <- coords98_tmp$sensitivity[1]
      leaderboard_top[[feature_type_tmp]]$Train_spec98[i] <- coords98_tmp$specificity[1]
      
      ####### get roc metrics of valid frame
      if(nrow(valid_frame) > 0){
        valid_score[, model_name] = as.numeric(as.data.frame(h2o.predict(model_i, newdata=as.h2o(valid_frame)))[,"Cancer"])
        roc_valid_tmp <- roc(valid_score$Train_Group,valid_score[,model_name],plot = F,levels = c("BLN", "Cancer"))
        
        coords95_valid_tmp <- coords(roc = roc_valid_tmp, x = coords95_tmp$threshold[1], input = "threshold", ret = "all")
        leaderboard_top[[feature_type_tmp]]$Valid_sens95[i] <- coords95_valid_tmp$sensitivity[1]
        leaderboard_top[[feature_type_tmp]]$Valid_spec95[i] <- coords95_valid_tmp$specificity[1]
        
        coords98_valid_tmp <- coords(roc = roc_valid_tmp, x = coords98_tmp$threshold[1], input = "threshold", ret = "all")
        leaderboard_top[[feature_type_tmp]]$Valid_sens98[i] <- coords98_valid_tmp$sensitivity[1]
        leaderboard_top[[feature_type_tmp]]$Valid_spec98[i] <- coords98_valid_tmp$specificity[1]
      }
    } # end of top5 base model iteration
    
    tmp_output <- leaderboard_top[[feature_type_tmp]] %>% mutate(feature_type = feature_type_tmp)
    if(nrow(leaderboard_top_output)){
      leaderboard_top_output <- rbind(leaderboard_top_output, tmp_output)
    }else{
      leaderboard_top_output <- tmp_output
    }
  
  } # end of feature type iteration
  
  ##### export layer 1 results
  write.csv(leaderboard_top_output, paste0(output_path,"Layer1_metrics.csv"),row.names = F)
  write.csv(train_score, paste0(output_path,"Layer1_train_score.csv"),row.names = F)
  if(nrow(valid_frame) > 0){
    write.csv(valid_score, paste0(output_path,"Layer1_valid_score.csv"),row.names = F)
  }
  
  h2o.removeAll()
  
  ##### start layer 2 stacking using layer 1 scores
  if(l2_stacking){
    predictor_l2 <- colnames(train_score) %>% setdiff(.,c("ID","Train_Group","fold_column"))
    predictor_l2 <- predictor_l2[grepl("_[1-5]$",predictor_l2)] #keep top 5 models
    response_l2 <- "Train_Group"
    

    print("------------------  starting L2 autoML  ------------------")
    
    model_l2 <- h2o.automl(x = predictor_l2,
                           y = "Train_Group",
                           
                           project_name = "Layer2",
                           
                           # nfolds = 5, 
                           fold_column = "fold_column",
                           
                           balance_classes = T,
                           
                           keep_cross_validation_predictions =T,
                           keep_cross_validation_fold_assignment = T,
                           
                           stopping_metric = "AUC", 
                           # stopping_tolerance = 0.01, 
                           max_runtime_secs_per_model=600, 
                           
                           include_algos = c("GLM","XGBoost","DRF"),
                           
                           training_frame = as.h2o(train_score), 

                           max_models = 50, 
                           seed = 99)
    
    ##### get leaderboard - layer 2 
    leaderboard_top_l2 <- model_l2@leaderboard %>% as.data.frame() %>%
      mutate(Train_AUC = NA, Valid_AUC = NA, 
             Train_cutoff95 = NA, Train_sens95 = NA, Train_spec95 = NA, Valid_sens95 = NA, Valid_spec95 = NA,
             Train_cutoff98 = NA, Train_sens98 = NA, Train_spec98 = NA, Valid_sens98 = NA, Valid_spec98 = NA) %>%
      head(3)
    
    for(i in 1:min(3,nrow(leaderboard_top_l2))){
      model_l2_i <- h2o.getModel(leaderboard_top_l2$model_id[i])
      h2o.saveModel(model_l2_i, path = paste0(output_path,"/","l2_",i), 
                    export_cross_validation_predictions = T, force = T)
      
      l2_model_name <- paste0("l2_",i)
      
      ####### get predict score and metrics - layer 2 train frame
      leaderboard_top_l2$Train_AUC[i] = h2o.auc(h2o.performance(model_l2_i, xval = T))
      
      train_score[, l2_model_name] = as.numeric(as.data.frame(h2o.getFrame(model_l2_i@model$cross_validation_holdout_predictions_frame_id$name))[,"Cancer"])
      
      roc_tmp_l2 <- roc(train_score$Train_Group,train_score[,l2_model_name],plot = F,levels = c("BLN", "Cancer"))
      coords95_tmp_l2 <- coords(roc = roc_tmp_l2, x = 'all', transpose = F, as.matrix = T) %>%
        as.data.frame() %>%
        filter(specificity >= 0.95) %>%
        filter(specificity <= min(specificity)) %>%
        arrange(desc(sensitivity))
      leaderboard_top_l2$Train_cutoff95[i] <- coords95_tmp_l2$threshold[1]
      leaderboard_top_l2$Train_sens95[i] <- coords95_tmp_l2$sensitivity[1]
      leaderboard_top_l2$Train_spec95[i] <- coords95_tmp_l2$specificity[1]
      
      coords98_tmp_l2 <- coords(roc = roc_tmp_l2, x = 'all', transpose = F, as.matrix = T) %>%
        as.data.frame() %>%
        filter(specificity >= 0.98) %>%
        filter(specificity <= min(specificity)) %>%
        arrange(desc(sensitivity))
      leaderboard_top_l2$Train_cutoff98[i] <- coords98_tmp_l2$threshold[1]
      leaderboard_top_l2$Train_sens98[i] <- coords98_tmp_l2$sensitivity[1]
      leaderboard_top_l2$Train_spec98[i] <- coords98_tmp_l2$specificity[1]
      
      ####### get predict score - layer 2 valid frame and metrics
      if(nrow(valid_frame) > 0){
        
        leaderboard_top_l2$Valid_AUC[i] = h2o.auc(h2o.performance(model_l2_i, newdata = as.h2o(valid_score)))
        
        valid_score[, l2_model_name] = as.numeric(as.data.frame(h2o.predict(model_l2_i, newdata=as.h2o(valid_score)))[,"Cancer"])
        roc_tmp_valid_l2 <- roc(valid_score$Train_Group,valid_score[,l2_model_name],plot = F,levels = c("BLN", "Cancer"))      
      
        coords95_valid_tmp_l2 <- coords(roc = roc_tmp_valid_l2, x = coords95_tmp_l2$threshold[1], input = "threshold", ret = "all")
        leaderboard_top_l2$Valid_sens95[i] <- coords95_valid_tmp_l2$sensitivity[1]
        leaderboard_top_l2$Valid_spec95[i] <- coords95_valid_tmp_l2$specificity[1]
        
        coords98_valid_tmp_l2 <- coords(roc = roc_tmp_valid_l2, x = coords98_tmp_l2$threshold[1], input = "threshold", ret = "all")
        leaderboard_top_l2$Valid_sens98[i] <- coords98_valid_tmp_l2$sensitivity[1]
        leaderboard_top_l2$Valid_spec98[i] <- coords98_valid_tmp_l2$specificity[1]
      }
      
    } # end of top 3 L2 model iteration
    write.csv(leaderboard_top_l2, paste0(output_path,"Layer2_metrics.csv"),row.names = F)
    
    ##### export layer 2 results
    write.csv(leaderboard_top_l2, paste0(output_path,"Layer2_metrics.csv"),row.names = F)
    write.csv(train_score, paste0(output_path,"Layer2_train_score.csv"),row.names = F)
    if(nrow(valid_frame) > 0){
      write.csv(valid_score, paste0(output_path,"Layer2_valid_score.csv"),row.names = F)
    }    
  } # end of if l2_stacking
  h2o.removeAll()
  # h2o.shutdown(prompt = F) 
} # end of model fitting function

