# functions
# build reference model ####
build.ref.model <- function(ref_vars, dataset){
  ref_model <- glm(as.formula(paste('disease ~ ', paste(ref_vars, sep="", collapse=" + ") )), data = dataset, family = binomial) 
  return(ref_model)
}
#
# build new models ####
build.new.models <- function(ref_vars, new_vars, dataset){
  new_models <- list()
  for(i in seq_along(new_vars)){
    new_models[[i]] <- glm(as.formula(paste0(paste('disease ~ ', paste(ref_vars, sep="", collapse=" + ")), '+', new_vars[i])),
                           data = dataset, family = binomial)
  }
  names(new_models) <- new_vars
  return(new_models)
}
#
# calculate raw results with repetitions for training and test datasets ####
# prev - prevalence of event class / percentage of the event class / imbalance 
# rep - number of repetitions

calculate.raw.results.w.test <- function(prev, rep, 
                                         ds_nonev, 
                                         ds_event, 
                                         ref_vars,
                                         new_vars,
                                         n = 300, n_test = 100, cutoff = 0.5){
  
  n_nonev <- nrow(ds_nonev)
  n_event <- nrow(ds_event)
  
  n0 <- n - prev * n
  n1 <- prev * n
  n0_test <- n_test - prev * n_test
  n1_test <- prev * n_test
  
  # data frame for convergence check
  convergence <- data.frame()
  
  # data frame for indices for each repetition 
  indices <- data.frame()
  
  # data frames for results for reference and new models
  raw_results_ref <- data.frame()
  raw_results_new <- data.frame()
  
  for(i in 1:rep){ 
    ## reference model ####
    ind_nonev <- sample(1:n_nonev, n0, replace = T)
    ind_event <- sample(1:n_event, n1, replace = T)
    ds <- rbind(ds_nonev[ind_nonev,], 
                ds_event[ind_event,])
    
    ref_model  <- build.ref.model(ref_vars, ds)
    y <- ref_model$y
    p_ref <- ref_model$fitted.values
    res_sq_ref <- (y - p_ref)^2
    
    BS_ref   <- sum(res_sq_ref) / n
    BS_0_ref <- sum(res_sq_ref[y == 0]) / n0
    BS_1_ref <- sum(res_sq_ref[y == 1]) / n1
    ROC_ref <- roc(y, p_ref, direction = '<')
    AUC_ref <- ROC_ref$auc
    
    TN_ref <- sum(y == 0 & p_ref <  cutoff)
    FP_ref <- sum(y == 0 & p_ref >= cutoff)
    FN_ref <- sum(y == 1 & p_ref <  cutoff)
    TP_ref <- sum(y == 1 & p_ref >= cutoff)
    MCC_denominator_ref <- (TP_ref+FP_ref) * (TP_ref+FN_ref) * (TN_ref+FP_ref) * (TN_ref+FN_ref)
    
    MCC_ref <- ifelse(MCC_denominator_ref == 0, NA, (TN_ref*TP_ref - FN_ref*FP_ref) / sqrt(MCC_denominator_ref))
    F1_ref <- (2 * TP_ref) / (2 * TP_ref + FP_ref + FN_ref)
    
    ## test reference model ####
    ds_nonev_test <- ds_nonev[-ind_nonev,]
    ds_event_test <- ds_event[-ind_event,]
    ind_nonev_test <- sample(1:nrow(ds_nonev_test), n0_test, replace = T)
    ind_event_test <- sample(1:nrow(ds_event_test), n1_test, replace = T)
    ds_test <- rbind(ds_nonev_test[ind_nonev_test,],
                     ds_event_test[ind_event_test,])
    
    y_test <- ds_test$disease
    p_ref_test <- predict(ref_model, ds_test, type = 'r')
    res_sq_ref_test <- (y_test - p_ref_test)^2
    
    BS_ref_test <- sum(res_sq_ref_test) / n_test
    BS_0_ref_test <- sum(res_sq_ref_test[y_test == 0]) / n0_test
    BS_1_ref_test <- sum(res_sq_ref_test[y_test == 1]) / n1_test
    ROC_ref_test <- roc(y_test, p_ref_test, direction = '<')
    AUC_ref_test <- ROC_ref_test$auc
    
    TN_ref_test <- sum(y_test == 0 & p_ref_test <  cutoff)
    FP_ref_test <- sum(y_test == 0 & p_ref_test >= cutoff)
    FN_ref_test <- sum(y_test == 1 & p_ref_test <  cutoff)
    TP_ref_test <- sum(y_test == 1 & p_ref_test >= cutoff)
    
    MCC_denominator_ref_test <- (TP_ref_test+FP_ref_test) * (TP_ref_test+FN_ref_test) * (TN_ref_test+FP_ref_test) * (TN_ref_test+FN_ref_test)
    MCC_ref_test <-  ifelse(MCC_denominator_ref_test == 0, NA, (TN_ref_test*TP_ref_test - FN_ref_test*FP_ref_test) / sqrt(MCC_denominator_ref_test) )
    F1_ref_test <- (2 * TP_ref_test) / (2 * TP_ref_test + FP_ref_test + FN_ref_test)
    
    temp_ref <- data.frame(iteration = i,
                           dataset = 'train',
                           prevalence = prev,
                           AUC_ref,
                           BS_ref,
                           BS_0_ref,
                           BS_1_ref,
                           TN_ref,
                           FP_ref,
                           FN_ref,
                           TP_ref,
                           MCC_ref,
                           F1_ref  ) 
    
    temp_ref_test <- data.frame(iteration = i,
                                dataset = 'test',
                                prevalence = prev,
                                AUC_ref = AUC_ref_test,
                                BS_ref = BS_ref_test,
                                BS_0_ref = BS_0_ref_test,
                                BS_1_ref = BS_1_ref_test,
                                TN_ref = TN_ref_test,
                                FP_ref = FP_ref_test,
                                FN_ref = FN_ref_test,
                                TP_ref = TP_ref_test,
                                MCC_ref = MCC_ref_test,
                                F1_ref = F1_ref_test )
    
    raw_results_ref <- rbind(raw_results_ref, temp_ref, temp_ref_test) 
    
    temp_ind <- data.frame(iteration = rep(i, n + n_test), 
                           dataset = c(rep('train', n), rep('test', n_test)),
                           prevalence = rep(prev, n + n_test), 
                           class = c(rep(0, n0), rep(1, n1), rep(0, n0_test), rep(1, n1_test)),
                           index = c(ind_nonev, ind_event, ind_nonev_test, ind_event_test) )
    indices <- rbind(indices, temp_ind)
    
    # # # # # # # # # # #  # #
    
    ## new models ####
    new_models <- build.new.models(ref_vars, new_vars, ds)
    
    for(j in seq_along(new_models)){
      conv_temp <- data.frame(iteration = i, 
                              prevalence = prev,
                              model = names(new_models)[j], 
                              conv_ref = ref_model$converged, 
                              conv = new_models[[j]]$converged)
      convergence <- rbind(convergence, conv_temp)
      ## BARBI ####
      p <- new_models[[j]]$fitted.values
      delta <- p - p_ref
      flag <- ifelse(delta > 0, 'up', 
                     ifelse(delta < 0, 'dw', 'c'))
      subclass <- ifelse(y == 0 & flag == 'dw',                      'nonev_be',
                         ifelse(y == 0 & flag == 'up',               'nonev_wo',
                                ifelse(y == 1 & flag == 'dw',        'event_wo',
                                       ifelse(y == 1 & flag == 'up', 'event_be', 'unknown'))))
      
      res_sq <- (y - p)^2
      
      SS_nonev_dw_ref<- sum(res_sq_ref[subclass == 'nonev_be'])
      SS_nonev_up_ref<- sum(res_sq_ref[subclass == 'nonev_wo'])
      SS_event_dw_ref<- sum(res_sq_ref[subclass == 'event_wo'])
      SS_event_up_ref<- sum(res_sq_ref[subclass == 'event_be'])
      SS_nonev_dw    <- sum(res_sq    [subclass == 'nonev_be'])
      SS_nonev_up    <- sum(res_sq    [subclass == 'nonev_wo'])
      SS_event_dw    <- sum(res_sq    [subclass == 'event_wo'])
      SS_event_up    <- sum(res_sq    [subclass == 'event_be'])
      
      BA_nonev_be <-  (SS_nonev_dw_ref - SS_nonev_dw) / n0
      BA_nonev_wo <- -(SS_nonev_up_ref - SS_nonev_up) / n0
      BA_event_wo <- -(SS_event_dw_ref - SS_event_dw) / n1
      BA_event_be <-  (SS_event_up_ref - SS_event_up) / n1
      
      BA_nonev <- BA_nonev_be - BA_nonev_wo
      BA_event <- BA_event_be - BA_event_wo
      
      BA_overall <- n0/n * BA_nonev + n1/n * BA_event
      
      RB_nonev_be <-  (SS_nonev_dw_ref - SS_nonev_dw) / (SS_nonev_dw_ref + SS_nonev_up_ref)
      RB_nonev_wo <- -(SS_nonev_up_ref - SS_nonev_up) / (SS_nonev_dw_ref + SS_nonev_up_ref)
      RB_event_wo <- -(SS_event_dw_ref - SS_event_dw) / (SS_event_dw_ref + SS_event_up_ref)
      RB_event_be <-  (SS_event_up_ref - SS_event_up) / (SS_event_dw_ref + SS_event_up_ref)
      
      RB_nonev <- RB_nonev_be - RB_nonev_wo
      RB_event <- RB_event_be - RB_event_wo
      
      RB_overall <- n0/n * RB_nonev + n1/n * RB_event
      
      II_nonev_be  <- sum(subclass == 'nonev_be') / n0
      II_nonev_wo  <- sum(subclass == 'nonev_wo') / n0
      II_event_wo  <- sum(subclass == 'event_wo') / n1
      II_event_be  <- sum(subclass == 'event_be') / n1
      
      II_nonev <- II_nonev_be - II_nonev_wo
      II_event <- II_event_be - II_event_wo
      
      II_overall <- n0/n * II_nonev + n1/n * II_event
      
      ## brier LRT ROC ####
      BS   <- sum(res_sq) / n
      BS_0 <- sum(res_sq[y == 0]) / n0
      BS_1 <- sum(res_sq[y == 1]) / n1
      
      delta_BS <- BS_ref - BS
      BSS <- 1 - BS / BS_ref
      
      LRT_statistic <- lrtest(ref_model, new_models[[j]])$Chisq[2]
      LRT_p_value   <- lrtest(ref_model, new_models[[j]])$`Pr(>Chisq)`[2]
      
      ROC <- roc(y, p, direction = '<')
      
      AUC <- ROC$auc
      AUC_statistic <- roc.test(ROC_ref, ROC)$statistic
      AUC_p_value   <- roc.test(ROC_ref, ROC)$p.value
      delta_AUC <- AUC - AUC_ref
      ## MCC F1 ####
      TN <- sum(y == 0 & p <  cutoff)
      FP <- sum(y == 0 & p >= cutoff)
      FN <- sum(y == 1 & p <  cutoff)
      TP <- sum(y == 1 & p >= cutoff)
      MCC_denominator <- (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)
      
      MCC <- ifelse(MCC_denominator == 0, NA, (TN*TP - FN*FP) / sqrt(MCC_denominator))
      F1 <- (2 * TP) / (2 * TP + FP + FN)
      
      delta_MCC <- MCC - MCC_ref
      delta_F1  <- F1  - F1_ref
      
      #
      ## test BARBI #### 
      p_test <- predict(new_models[[j]], ds_test, type = 'r')
      delta_test <- p_test - p_ref_test
      flag_test <- ifelse(delta_test > 0, 'up', 
                          ifelse(delta_test < 0, 'dw', 'c'))
      subclass_test <- ifelse(y_test == 0 & flag_test == 'dw',                 'nonev_be',
                              ifelse(y_test == 0 & flag_test == 'up',               'nonev_wo',
                                     ifelse(y_test == 1 & flag_test == 'dw',        'event_wo',
                                            ifelse(y_test == 1 & flag_test == 'up', 'event_be', 'unknown'))))
      
      res_sq_test <- (y_test - p_test)^2
      
      SS_nonev_dw_ref_test<- sum(res_sq_ref_test[subclass_test == 'nonev_be'])
      SS_nonev_up_ref_test<- sum(res_sq_ref_test[subclass_test == 'nonev_wo'])
      SS_event_dw_ref_test<- sum(res_sq_ref_test[subclass_test == 'event_wo'])
      SS_event_up_ref_test<- sum(res_sq_ref_test[subclass_test == 'event_be'])
      SS_nonev_dw_test    <- sum(res_sq_test    [subclass_test == 'nonev_be'])
      SS_nonev_up_test    <- sum(res_sq_test    [subclass_test == 'nonev_wo'])
      SS_event_dw_test    <- sum(res_sq_test    [subclass_test == 'event_wo'])
      SS_event_up_test    <- sum(res_sq_test    [subclass_test == 'event_be'])
      
      BA_nonev_be_test <-  (SS_nonev_dw_ref_test - SS_nonev_dw_test) / n0_test
      BA_nonev_wo_test <- -(SS_nonev_up_ref_test - SS_nonev_up_test) / n0_test
      BA_event_wo_test <- -(SS_event_dw_ref_test - SS_event_dw_test) / n1_test
      BA_event_be_test <-  (SS_event_up_ref_test - SS_event_up_test) / n1_test
      
      BA_nonev_test <- BA_nonev_be_test - BA_nonev_wo_test
      BA_event_test <- BA_event_be_test - BA_event_wo_test
      
      BA_overall_test <- n0/n * BA_nonev_test + n1/n * BA_event_test
      
      RB_nonev_be_test <-  (SS_nonev_dw_ref_test - SS_nonev_dw_test) / (SS_nonev_dw_ref_test + SS_nonev_up_ref_test)
      RB_nonev_wo_test <- -(SS_nonev_up_ref_test - SS_nonev_up_test) / (SS_nonev_dw_ref_test + SS_nonev_up_ref_test)
      RB_event_wo_test <- -(SS_event_dw_ref_test - SS_event_dw_test) / (SS_event_dw_ref_test + SS_event_up_ref_test)
      RB_event_be_test <-  (SS_event_up_ref_test - SS_event_up_test) / (SS_event_dw_ref_test + SS_event_up_ref_test)
      
      RB_nonev_test <- RB_nonev_be_test - RB_nonev_wo_test
      RB_event_test <- RB_event_be_test - RB_event_wo_test
      
      RB_overall_test <- n0/n * RB_nonev_test + n1/n * RB_event_test
      
      II_nonev_be_test  <- sum(subclass_test == 'nonev_be') / n0_test
      II_nonev_wo_test  <- sum(subclass_test == 'nonev_wo') / n0_test
      II_event_wo_test  <- sum(subclass_test == 'event_wo') / n1_test
      II_event_be_test  <- sum(subclass_test == 'event_be') / n1_test
      
      II_nonev_test <- II_nonev_be_test - II_nonev_wo_test
      II_event_test <- II_event_be_test - II_event_wo_test
      
      II_overall_test <- n0/n * II_nonev_test + n1/n * II_event_test
      
      ## test Brier ROC ####
      BS_test   <- sum(res_sq_test) / n_test
      BS_0_test <- sum(res_sq_test[y_test == 0]) / n0_test
      BS_1_test <- sum(res_sq_test[y_test == 1]) / n1_test
      
      delta_BS_test <- BS_ref_test - BS_test
      BSS_test <- 1 - BS_test / BS_ref_test
      
      ROC_test <- roc(y_test, p_test, direction = '<')
      AUC_test <- ROC_test$auc
      delta_AUC_test <- AUC_test - AUC_ref_test
      
      ## test MCC F1 ####
      TN_test <- sum(y_test == 0 & p_test <  cutoff)
      FP_test <- sum(y_test == 0 & p_test >= cutoff)
      FN_test <- sum(y_test == 1 & p_test <  cutoff)
      TP_test <- sum(y_test == 1 & p_test >= cutoff)
      MCC_denominator_test <- (TP_test+FP_test) * (TP_test+FN_test) * (TN_test+FP_test) * (TN_test+FN_test)
      
      MCC_test <- ifelse(MCC_denominator_test == 0, NA, (TN_test*TP_test - FN_test*FP_test) / sqrt(MCC_denominator_test))
      F1_test <- (2 * TP_test) / (2 * TP_test + FP_test + FN_test)
      
      delta_MCC_test <- MCC_test - MCC_ref_test
      delta_F1_test  <- F1_test  - F1_ref_test
      
      ## temp results ####
      temp <- data.frame(iteration = i,
                         dataset = 'train',
                         prevalence = prev,
                         model = names(new_models)[j], 
                         BA_nonev_be,
                         BA_nonev_wo,
                         BA_event_wo,
                         BA_event_be,
                         BA_nonev,
                         BA_event,
                         BA_overall,
                         RB_nonev_be,
                         RB_nonev_wo,
                         RB_event_wo,
                         RB_event_be,
                         RB_nonev,
                         RB_event,
                         RB_overall,
                         II_nonev_be,
                         II_nonev_wo,
                         II_event_wo,
                         II_event_be,
                         II_nonev,
                         II_event,
                         II_overall,
                         BS  , 
                         BS_0, 
                         BS_1, 
                         delta_BS,
                         BSS,
                         LRT_statistic, # NA for test
                         LRT_p_value  , # NA for test
                         AUC          ,
                         AUC_statistic, # NA for test
                         AUC_p_value  , # NA for test
                         delta_AUC    , 
                         TN,
                         FP,
                         FN,
                         TP,
                         MCC,
                         F1,
                         delta_MCC,
                         delta_F1 
      )
      
      temp_test <- data.frame(iteration = i,
                              dataset = 'test',
                              prevalence = prev,
                              model = names(new_models)[j],
                              BA_nonev_be = BA_nonev_be_test,
                              BA_nonev_wo = BA_nonev_wo_test,
                              BA_event_wo = BA_event_wo_test,
                              BA_event_be = BA_event_be_test,
                              BA_nonev = BA_nonev_test,
                              BA_event = BA_event_test,
                              BA_overall = BA_overall_test,
                              RB_nonev_be = RB_nonev_be_test,
                              RB_nonev_wo = RB_nonev_wo_test,
                              RB_event_wo = RB_event_wo_test,
                              RB_event_be = RB_event_be_test,
                              RB_nonev = RB_nonev_test,
                              RB_event = RB_event_test,
                              RB_overall = RB_overall_test,
                              II_nonev_be = II_nonev_be_test,
                              II_nonev_wo = II_nonev_wo_test,
                              II_event_wo = II_event_wo_test,
                              II_event_be = II_event_be_test,
                              II_nonev = II_nonev_test,
                              II_event = II_event_test,
                              II_overall = II_overall_test,
                              BS = BS_test  ,
                              BS_0 = BS_0_test,
                              BS_1 = BS_1_test,
                              delta_BS = delta_BS_test,
                              BSS = BSS_test,
                              LRT_statistic = NA,
                              LRT_p_value = NA,
                              AUC = AUC_test,
                              AUC_statistic = NA,
                              AUC_p_value = NA,
                              delta_AUC = delta_AUC_test, 
                              TN = TN_test,
                              FP = FP_test,
                              FN = FN_test,
                              TP = TP_test,
                              MCC = MCC_test,
                              F1 = F1_test,
                              delta_MCC = delta_MCC_test,
                              delta_F1 = delta_F1_test 
      )
      
      raw_results_new <- rbind(raw_results_new, temp, temp_test)
    } # end inner for loop
    
  } # end outer for loop 
  
  result <- list(raw_results_ref = raw_results_ref,
                 raw_results_new = raw_results_new,
                 convergence = convergence,
                 indices = indices)
  
  return(result)
  
} # end function
#
