
#### assemble 5X5 CV results and export the final score for each sample

assemble_cv_results <- function(target_path = "./cfDNA_internal_validation/"){
    round_path <- list.files(target_path)
    round_path <- round_path[grepl("round[0-9]",round_path)]
    validfold_score <- data.frame()
    for(i in 1:length(round_path)){
      
      fold_list <- paste0("fold",c(1:5))
      validfold_tmp_rbind <- NULL
      for(j in 1:length(fold_list)){
        
        path_tmp <- paste0(target_path,round_path[i],"/",fold_list[j],"/")
        
        validfold_tmp <- read.csv(paste0(path_tmp,"Layer2_valid_score.csv")) %>%
          mutate(fold = fold_list[j]) %>%
          mutate(l2_mean_score = (l2_1 + l2_2 + l2_3)/3) %>%
          select(ID, Train_Group, fold, l2_mean_score) %>%
          rename(!!paste0(round_path[i],"_score") := "l2_mean_score") %>%
          rename(!!paste0(round_path[i],"_fold") := "fold")
        
        validfold_tmp_rbind <- rbind(validfold_tmp_rbind, validfold_tmp)
      }
      if(nrow(validfold_score)){
        
        validfold_score <- left_join(validfold_score, validfold_tmp_rbind, multiple = "all",
                                     by = c("ID"="ID","Train_Group"="Train_Group"))
      }else{
        validfold_score <- validfold_tmp_rbind
      }
    }
    validfold_score[,"Final_score"] = rowMeans(validfold_score[,paste0(round_path,"_score")])
    
    write.csv(validfold_score, paste0(target_path,"Final_score.csv"),row.names = F)
    
    return(validfold_score)
    
} # end of function
