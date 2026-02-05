
#################################################
############### PREDICTION MODEL ################
#################################################
library(RANN)
library('dplyr')
library("tidyr")
library("caret")
library("doParallel")
library("Metrics")
library("ggplot2")
library("TSrepr") #mape

setwd("/Users/ahhyun/Desktop/SI Predction/abct")
# To do:
# Remove Participants with more than 75% of missing data 
      # only include PTs over than 75%
# Decide if missing EMAs will be filled with placeholders 
      # Done
# Are first 7 days already excluded?
      # Done
# Check if EMA predictors are lagged (not from the same time point)
      # Done!
# Check MAPE (probably sth is wrong)

# Add parallelized version (Anna can do that)

# Am I running for different type of model
# H1: Each individualized model (i.e., EMA, passive sensing, combined) will predict SI more accurately than a simple baseline model that uses each individuals mean SI levels from the training set as the predictor.
# Idiogrpahic model vs nomothetic model
# H2: We expect the combined model to outperform the ema only and passive sensing only model.
# Within idiographic model, combined model will outperform 

set.seed(12361488)

############# trim down the data set ############
# remove PT with less than 75% of compliance rate # already did it when I merge the data
# Put placeholders for missing EMA data
# Exclude the first 7 days of data
# Q9 lag
# simplify the names
data <- read.csv("data.csv") 

ema_levels <- c("morning", "afternoon", "evening")

# Put placeholders for missing EMA data
data_filtered <- data %>%
  # parse "6/15/2021" style dates correctly
  mutate(
    date_rep        = as.Date(date_rep, format = "%m/%d/%Y"),
    ema_time_of_day = factor(ema_time_of_day, levels = ema_levels)
  ) %>%
  
  group_by(uid) %>%
  # skip any participants where date_rep is entirely NA (avoids seq() error)
  filter(!all(is.na(date_rep))) %>%
  
  # create consecutive dates + 3 daily EMA slots
  complete(
    date_rep = seq(
      min(date_rep, na.rm = TRUE),
      max(date_rep, na.rm = TRUE),
      by = "day"
    ),
    ema_time_of_day = ema_levels   # morning/afternoon/evening
  ) %>%
  arrange(uid, date_rep, ema_time_of_day) %>%
  ungroup() 
  

View(data_filtered)

# Exclude the first 7 days of data
data <- data_filtered %>%
  group_by(uid) %>%
  mutate(first_date = min(date_rep, na.rm = TRUE)) %>%
  filter(date_rep >= first_date + 7) %>%
  select(-first_date) %>%
  ungroup() %>%
  mutate(q9_lag = lag(q9, 1))%>%
  group_by(uid, date_rep) %>%
  mutate(
    # get morning q3 (can be NA), or NA if morning row doesn't exist
    q3_morning = dplyr::first(q3[ema_time_of_day == "morning"], 
                              default = NA_real_),
    # overwrite all q3 for that day with the morning value (including NA)
    q3 = q3_morning
  ) %>%
  select(-q3_morning) %>%
  ungroup() 
View(data)
data <- data %>%
  filter(compliance >= 0.75) %>%
  filter(!is.na(q9_lag))

View(data)

length(unique(data$uid)) #227
summary(data)
write.csv(data, "data.csv")


########## Choose predictors & Outcome  ########
#load data
data <- read.csv("data.csv")
data <- data_raw %>%
  rename(day = ema_time_of_day) %>%       # ensure 'day' exists
  filter(!is.na(q9_lag)) %>%              # outcome must exist
  select(
    -any_of(c("date_rep", "ema_time"))    # drop only if present
  ) %>%
  group_by(uid) %>%                       # scale per participant
  mutate(
    across(
      .cols = where(is.numeric) & !matches("^uid$|^x$"),
      .fns  = ~ as.numeric(scale(.)),
      .names = "{.col}"
    )
  ) %>%
  ungroup()
View(data)

predictors <- "onlyESM" # "all", onlypassive", "onlyESM"(select included set of predictor variables)
outcome <- "q9_lag" 

############ Load and Prepare Data ###############
if(predictors == "onlypassive") {
  
  data <- data %>%
    select(
     uid, day,
      mean_1h_hr, MAD_1h_hr, P90P10_pos_1h_hr,
      mean_4h_hr, MAD_4h_hr, P90P10_pos_4h_hr,
      mean_1d_hr, MAD_1d_hr, P90P10_pos_1d_hr,
      mean_7d_hr, MAD_7d_hr, P90P10_pos_7d_hr,
      
      mean_1h_gps, MAD_1h_gps, P90P10_pos_1h_gps,
      mean_4h_gps, MAD_4h_gps, P90P10_pos_4h_gps,
      mean_1d_gps, MAD_1d_gps, P90P10_pos_1d_gps,
      mean_7d_gps, MAD_7d_gps, P90P10_pos_7d_gps,
      
      mean_1h_conv, MAD_1h_conv, P90P10_pos_1h_conv,
      mean_4h_conv, MAD_4h_conv, P90P10_pos_4h_conv,
      mean_1d_conv, MAD_1d_conv, P90P10_pos_1d_conv,
      mean_7d_conv, MAD_7d_conv, P90P10_pos_7d_conv,
      
      mean_1h_bat, MAD_1h_bat, P90P10_pos_1h_bat,
      mean_4h_bat, MAD_4h_bat, P90P10_pos_4h_bat,
      mean_1d_bat, MAD_1d_bat, P90P10_pos_1d_bat,
      mean_7d_bat, MAD_7d_bat, P90P10_pos_7d_bat,
      
      mean_1h_sms, MAD_1h_sms, P90P10_pos_1h_sms,
      mean_4h_sms, MAD_4h_sms, P90P10_pos_4h_sms,
      mean_1d_sms, MAD_1d_sms, P90P10_pos_1d_sms,
      mean_7d_sms, MAD_7d_sms, P90P10_pos_7d_sms,
      
      mean_1h_call, MAD_1h_call, P90P10_pos_1h_call,
      mean_4h_call, MAD_4h_call, P90P10_pos_4h_call,
      mean_1d_call, MAD_1d_call, P90P10_pos_1d_call,
      mean_7d_call, MAD_7d_call, P90P10_pos_7d_call
    )
} else if (predictors == "onlyESM") {
  
  data <- data %>%
    select(
      day,
      uid,
      q1:q9,
      g1:g2, q9_lag
    )
}

###### Helper function #######
# This function splits the data into two sets train+val set and test set
# train + val will be splitted later 

SlidingWindow_CV <- function(data, origin, horizon){
  # origin I put 54*3 (162)
  # horizon I put 3
  samplesize <- nrow(data)
  trainindex <- list()
  testindex <- list()
  
  trainindex[[1]] <- 1:origin
  testindex[[1]] <- (origin + 1):(origin+horizon)
  
  counter <- 1
  index <- testindex[[1]][horizon]+1
  
  while(testindex[[counter]][horizon] < (samplesize-horizon)){
    index <- testindex[[counter]][horizon]+1
    trainindex[[counter+1]] <- (index-origin+1):index-1
    testindex[[counter+1]] <- (index):(index+horizon-1)
    counter <- counter + 1
  }
  return(list(trainindex,testindex))
}


#################### START PREDICTIONS ###################

# Error in if (rss/tss == 0) rsq <- NA : 
#   missing value where TRUE/FALSE needed

Overall_results_allparticipants <- list() # Saves the prediction performance
Overall_pred_allparticipants <- list() # Saves all predicted and observed values

PNumbers <- unique(data$uid)

#for(h in 1:length(PNumbers)){
for(h in 10:15){
  ##### Create empty lists
  Overall_results <- list() # Saves the prediction performance
  Overall_pred <- list() # Saves all predicted and observed values
  counter <- 1


  Participant1 = data[data["uid"] == PNumbers[h],  ] # We will loop through each participant

  Features <- Participant1 %>% select(-day, -uid, -q9_lag) %>% colnames() # remove outcome and non needed variables here here

  set.seed(12361488)

  timeSlices <- SlidingWindow_CV(Participant1,54*3,3)   # Split into train+val and test set
  trainSlices <- timeSlices[[1]] # train and validation set
  testSlices <- timeSlices[[2]] # test set

  lengthtest <- length(unlist(testSlices))

  Results_pred <- data.frame(index = unlist(testSlices), pred = rep(NA,lengthtest), true = rep(NA,lengthtest))

  fitControl <- trainControl(method="timeslice", horizon = 3*10, initialWindow = 132, fixedWindow = TRUE) #Splits train+val into train and validation set, horizon predicts 3


  for(i in 1:length(trainSlices)){


        ### Create prediction model
        formula2 <- as.formula(paste("q9_lag ~", paste(colnames(Participant1[,colnames(Participant1) %in% Features]), collapse = "+")))

        # see https://topepo.github.io/caret/
        try <- try(model <- caret::train(formula2, data=Participant1[trainSlices[[i]],],trControl = fitControl, method="rf",metric = "MAE", preProcess = c("knnImpute"),
                                         na.action = na.pass)) # You could add hyperparameters here e.g., mtry: number of randomly selected predictors, you could also add other preProcess options e.g., nzv

        # Next we store the results from the prediction model by making predictions for our test set. The predict() function automatically uses the "best" model (based on leave last three out from the train+val set)
        dataindex <- which(Results_pred$index == testSlices[[i]])
        Results_pred$pred[dataindex] <- NA

        if(any(class(try) != "try-error")){
          try(Results_pred$pred[dataindex] <- predict(model,Participant1[testSlices[[i]],], na.action = na.pass)) #final model automatically used with predict (na.pass, because missing values are getting imputed)
          Results_pred$true[dataindex] <- Participant1$q9_lag[testSlices[[i]]] # Ahhyun: Changed to Q9_lag
        }


        Results_pred$Participant <- PNumbers[h]

      }

      # Performance rf
      rss <- sum((Results_pred$pred - Results_pred$true) ^ 2,na.rm = TRUE)  # residual sum of squares
      tss <- sum((Results_pred$true - mean(Results_pred$true, na.rm = TRUE))^ 2,na.rm = TRUE)  # total sum of squares
      rsq <- 1 - rss/tss
      if(rss/tss == 0) rsq <- NA

      complete <- complete.cases(Results_pred$true,Results_pred$pred)

      mape <- mape(Results_pred$true[complete],Results_pred$pred[complete])
      smape <- smape(Results_pred$true[complete],Results_pred$pred[complete])
      rae <- rae(Results_pred$true[complete],Results_pred$pred[complete])
      mae <- mae(Results_pred$true[complete],Results_pred$pred[complete])
      cor <- cor(Results_pred$true,Results_pred$pred, use = "pairwise.complete.obs")

      Results <- cbind(rsq,cor,mae,mape,rae,smape,PNumbers[h])

      print(Results)

      Overall_results[[counter]] <- Results
      Overall_pred[[counter]] <- Results_pred
      counter <- counter + 1


}

Overall_results_allparticipants[[h]] <- as.data.frame(do.call("rbind",Overall_results))
Overall_pred_allparticipants[[h]] <- as.data.frame(do.call("rbind",Overall_pred)) #Combine lists into one dataframe

##### Save results
Overall_results_allparticipants <- as.data.frame(do.call("rbind",Overall_results_allparticipants)) # Save performance
Overall_pred_allparticipants <- as.data.frame(do.call("rbind",Overall_pred_allparticipants)) # Save predictions

View(Overall_results_allparticipants)
View(Overall_pred_allparticipants)
write.csv(Overall_results_allparticipants,file = paste("Results/OverallResults","_step",step,"_",dataset,"_",lag,"_",outcome,"_",predictors,".csv",sep = ""))
write.csv(Overall_pred_allparticipants,file = paste("Results/OverallPred","_step",step,"_",dataset,"_",lag,"_",outcome, "_",predictors,".csv",sep = ""))



### Ahhyun: This still needs to be adapted!

#######################################################
############### USE MEAN AS PREDICTION ################
#######################################################
# change mood to q9_log

Overall_results_allparticipants_mean <- list() # Stores the prediction performance
Overall_pred_allparticipants_mean <- list() # create a list to store the predicted and observed values

for(h in 1:length(PNumbers)){
  
  counter <- 1
  
  Overall_results_mean <- list()
  Overall_pred_mean <- list()
  
  # Participant1 = data[data["ParticipantNumber"] == PNumbers[h] & data["timescale_beforeESM"] == 6,  ] # Here the level of aggregation is not important (because we only use positive affect). Thus, we could have also used a different level than 6h
  #  ## Ahhyun: we already have the calcuated mean? 
  # 
  # q9_lag <- Participant1[,outcome]
  # Participant1$q9_lag <- q9_lag + 1
  # 
  for(k in c(15,20,30)){
    
    source("BlockedCV.R")
    timeSlices <- SlidingWindow_CV(Participant1,k,step)   
    trainSlices <- timeSlices[[1]]
    testSlices <- timeSlices[[2]]
    
    
    lengthtest <- length(unlist(testSlices))
    
    Results_pred_mean <- data.frame(index = unlist(testSlices),pred = rep(NA,lengthtest),true = rep(NA,lengthtest))
    
    for(i in 1:length(trainSlices)){
      dataindex <- which(Results_pred_mean$index == testSlices[[i]])
      
      Results_pred_mean$pred[dataindex] <- mean(Participant1$q9_lag[trainSlices[[i]]], na.rm = TRUE) #we use the Q9 as a prediction.
      Results_pred_mean$true[dataindex] <- Participant1$q9_lag[testSlices[[i]]]
      Results_pred_mean$w <- k
      Results_pred_mean$Participant <- PNumbers[h]
      
      
    }
    
    # Performance
    rss <- sum((Results_pred_mean$pred - Results_pred_mean$true) ^ 2,na.rm = TRUE)  # residual sum of squares
    tss <- sum((Results_pred_mean$true - mean(Results_pred_mean$true, na.rm = TRUE))^ 2,na.rm = TRUE)  # total sum of squares
    rsq <- 1 - rss/tss
    
    complete <- complete.cases(Results_pred_mean$true,Results_pred_mean$pred)
    
    mape <- mape(Results_pred_mean$true[complete],Results_pred_mean$pred[complete])
    smape <- smape(Results_pred_mean$true[complete],Results_pred_mean$pred[complete])
    rae <- rae(Results_pred_mean$true[complete],Results_pred_mean$pred[complete])
    mae <- mae(Results_pred_mean$true[complete],Results_pred_mean$pred[complete])
    cor <- cor(Results_pred_mean$true,Results_pred_mean$pred, use = "pairwise.complete.obs")
    
    
    Overall_results_mean[[counter]] <- cbind(rsq,cor,mae,mape,rae,smape,k,n,PNumbers[h])
    Overall_pred_mean[[counter]] <- Results_pred_mean
    counter <- counter + 1
  }
  
  Overall_results_allparticipants_mean[[h]] <- as.data.frame(do.call("rbind",Overall_results_mean)) 
  Overall_pred_allparticipants_mean[[h]] <- as.data.frame(do.call("rbind",Overall_pred_mean)) #Combine lists into one dataframe
}

##### Save results

Overall_results_allparticipants_mean <- as.data.frame(do.call("rbind",Overall_results_allparticipants_mean))
Overall_pred_allparticipants_mean <- as.data.frame(do.call("rbind",Overall_pred_allparticipants_mean))
write.csv(Overall_results_allparticipants_mean,file = paste("Results/OverallResultsMEAN","_step",step,"_",outcome,".csv",sep = ""))
write.csv(Overall_pred_allparticipants_mean,file = paste("Results/OverallPredMEAN","_step",step,"_",outcome,".csv",sep = ""))


###########################################################################
############### USE MEAN + Dynamic indices PREDICTION #################
###########################################################################

###########################################################################
############### USE Combined PREDICTION (EMA + Passive) #################
###########################################################################