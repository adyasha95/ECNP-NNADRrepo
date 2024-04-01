
library(mlr)
library(tidyverse)
library(reshape2)
library(e1071)
library(caret)
library(ggplot2)
library(dplyr)
library(future)
library(parallel)

setwd('/volume/mitnvp1_scratch/AK_Testing/ECNP')

load("/volume/projects/JF_ECNP/Analysis/HC_SCZ_manuscript/Models_new/regression_age_unif_dist_age_sel_PID.RData")

data0 <- read.csv("/volume/projects/JF_ECNP/Analysis/HC_SCZ_manuscript/Models_new/classif_noICA_all_atlas_HC_SCZ_check_mri.csv")


data_matched <- data0 %>% 
  select(c("PID","SITE_ID", "AGE", "SEX","Diagnosis","TIV","QC_grade", starts_with("schaefer200"),starts_with("aal3"),starts_with("hammers")))
unique(data_matched$SITE_ID)
unique(data_matched$Diagnosis)
# Create a subset for the classification task (without age and gender)
data_matched <- data_matched %>% mutate_if(is.character, as.numeric)
dim(data_matched)
data_matched <-na.omit(data_matched , cols = c(starts_with("schaefer200")))
dim(data_matched)

dim(filter(data_matched,SITE_ID==1))
dim(filter(data_matched,SITE_ID==8))
dim(filter(data_matched,SITE_ID==6))
dim(filter(data_matched,SITE_ID==9))
TIV_matched <- data_matched$TIV
data_neuroimg <- data_matched  %>% 
  select(c(starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
data_neuroimg<-data.frame((data_neuroimg/TIV_matched)*1000)
min(data_neuroimg)
max(data_neuroimg)

data_neuroimg <- cbind(data_matched$SITE_ID,data_matched$AGE, data_matched$Diagnosis,data_matched$PID, data_neuroimg)
data_neuroimg$AGE <- data_neuroimg$`data_matched$AGE`
data_neuroimg$Diagnosis <- data_neuroimg$`data_matched$Diagnosis`
data_neuroimg$SITE_ID <- data_neuroimg$`data_matched$SITE_ID`
data_neuroimg$PID <- data_neuroimg$`data_matched$PID`
data_neuroimg<-subset(data_neuroimg, SITE_ID %in% c(1,2,3,6,8,9))

unique(data_neuroimg$Diagnosis)
unique(data_neuroimg$SITE_ID)
data_neuroimg$HC_idx <- ifelse(data_neuroimg$Diagnosis == 0, 1, 0)

# dummy code sites
#data_neuroimg$LMU <- ifelse(data_neuroimg$SITE_ID == 1, 1, 0)
#data_neuroimg$Edinburgh <- ifelse(data_neuroimg$SITE_ID == 3, 1, 0)
#data_neuroimg$Oslo <- ifelse(data_neuroimg$SITE_ID == 9, 1, 0)
#data_neuroimg$LMU
data_neuroimg_hc<-filter(data_neuroimg, Diagnosis == 0)
dim(data_neuroimg_hc)
data_neuroimg_regr <- data_neuroimg_hc %>% 
  select(c("PID","SITE_ID","HC_idx","AGE", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
data_neuroimg_regr_unif<-data_neuroimg_regr[data_neuroimg_regr$PID %in% sel_PID, ]
data_neuroimg_regr_rem <- data_neuroimg_regr[!data_neuroimg_regr$PID %in% sel_PID, ]
dim(data_neuroimg_regr_rem)
#unif sample
data_neuroimg_regr<-data_neuroimg_regr_unif
data_neuroimg_regr$PID<-NULL
dim(data_neuroimg_regr)
data_neuroimg_regr <- data_neuroimg_regr %>% 
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
#for remaining HC applicaion
data_neuroimg_regr_rem$PID<-NULL
dim(data_neuroimg_regr_rem)
data_neuroimg_regr_rem <- data_neuroimg_regr_rem %>% 
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))


# Filter the patient group (Diagnosis == 1) from the original data_neuroimg
data_neuroimg_patients <- filter(data_neuroimg, Diagnosis == 1)
dim(data_neuroimg_patients)
# Select the relevant features from the patient group for prediction
data_neuroimg_patients_regr <- data_neuroimg_patients %>%
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
#scz all
data_neuroimg_patients_scz_all<-data_neuroimg_patients_regr %>%
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
#schz MUC
data_neuroimg_patients_scz_MUC<-filter(data_neuroimg_patients_regr, SITE_ID==1)
data_neuroimg_patients_scz_MUC<-data_neuroimg_patients_scz_MUC %>%
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
#scz VER
data_neuroimg_patients_scz_VER<-filter(data_neuroimg_patients_regr, SITE_ID==8)
data_neuroimg_patients_scz_VER<-data_neuroimg_patients_scz_VER %>%
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
#scz OSL
data_neuroimg_patients_scz_OSL<-filter(data_neuroimg_patients_regr, SITE_ID==9)
data_neuroimg_patients_scz_OSL<-data_neuroimg_patients_scz_OSL %>%
  select(c("AGE","HC_idx","SITE_ID", starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))


dim(data_neuroimg_regr)
dim(data_neuroimg_regr_rem)

dim(data_neuroimg_patients_scz_all)
dim(data_neuroimg_patients_scz_MUC)
dim(data_neuroimg_patients_scz_VER)
dim(data_neuroimg_patients_scz_OSL)

makePreprocWrapperGlobMeanDiff <- function(learner, group.column = "TRUE", subgroup.column = "TRUE") {
  trainfun <- function(data, target, args = list(group.column = group.column, subgroup.column = subgroup.column)) {
    # check if control.col is not all 0 (i.e. mean difference should be computed 
    # across the whole group)
    col.names <- colnames(data)
    ids <- 1:nrow(data)
    data$artificial.id <- ids
    if (args$subgroup.column %in% col.names){
      # subset data 
      suppressMessages(hdata.part1.long <- data %>% filter(!!as.symbol(args$subgroup.column) == 1) %>% 
                         dplyr::select(-!!as.symbol(args$subgroup.column)) %>%
                         pivot_longer(
                           cols = -c(!!as.symbol(args$group.column), artificial.id, !!as.symbol(target)), 
                           names_to = "feature", 
                           values_to = "value")%>% 
                         group_by(feature) %>% 
                         mutate(overall.mean = mean(value, na.rm = TRUE)) %>% 
                         group_by(feature, !!as.symbol(args$group.column)) %>% 
                         mutate(group.mean = mean(value, na.rm = TRUE)) %>% 
                         ungroup() %>%
                         mutate(h.value = value - (group.mean - overall.mean)))
      
      suppressMessages(hdata.part1 <- hdata.part1.long %>%
                         dplyr::select(-c(value, overall.mean, group.mean)) %>% 
                         pivot_wider(names_from = feature, values_from = h.value))
      
      suppressMessages(estimates <- hdata.part1.long %>% 
                         dplyr::select(!!as.symbol(args$group.column), feature, overall.mean, group.mean) %>% 
                         distinct())
      
      suppressMessages(hdata.part2 <- data %>% filter(!!as.symbol(args$subgroup.column) == 0) %>%
                         dplyr::select(-!!as.symbol(args$subgroup.column)) %>%
                         pivot_longer(
                           cols = -c(!!as.symbol(args$group.column), artificial.id, !!as.symbol(target)), 
                           names_to = "feature", 
                           values_to = "value") %>% 
                         left_join(estimates, ) %>% 
                         mutate(h.value = value - (group.mean - overall.mean)) %>% 
                         dplyr::select(-c(value, overall.mean, group.mean)) %>% 
                         pivot_wider(names_from = feature, values_from = h.value))
      
      suppressMessages(data <- hdata.part1 %>% full_join(hdata.part2) %>% 
                         arrange(artificial.id) %>% 
                         dplyr::select(-artificial.id, -!!as.symbol(args$group.column)))
      
      
    } else {
      suppressMessages(hdata.long <- data %>% 
                         pivot_longer(
                           cols = -c(!!as.symbol(args$group.column), artificial.id, !!as.symbol(target)), 
                           names_to = "feature", 
                           values_to = "value") %>% 
                         group_by(feature) %>% 
                         mutate(overall.mean = mean(value, na.rm = TRUE)) %>% 
                         group_by(feature, !!as.symbol(args$group.column)) %>% 
                         mutate(group.mean = mean(value, na.rm = TRUE)) %>% 
                         ungroup() %>%
                         mutate(h.value = value - (group.mean - overall.mean)))  
      
      suppressMessages(data <- hdata.long %>% 
                         dplyr::select(-c(value, overall.mean, group.mean)) %>% 
                         pivot_wider(names_from = feature, values_from = h.value) %>% 
                         arrange(artificial.id) %>% 
                         dplyr::select(-artificial.id, -!!as.symbol(args$group.column)))
      
      suppressMessages(estimates <- hdata.long %>% 
                         dplyr::select(!!as.symbol(args$group.column), feature, overall.mean, group.mean) %>% 
                         distinct())
    }
    
    return(list(data = data, control = list(estimates = estimates, args = args)))
  }
  
  predictfun <- function(data, target, args, control) {
    
    #print(args$group.column)
    # Your prediction logic here
    col.names <- colnames(data)
    ids <- 1:nrow(data)
    data$artificial.id <- ids
    
    if(!(setequal(control$estimates[[control$args$group.column]], data[[args$group.column]]))) {
      if (control$args$subgroup.column %in% col.names) {
        # subset data 
        suppressMessages(hdata.pred <- data %>% filter(!!as.symbol(control$args$group.column) %in% setdiff(data[[control$args$group.column]], control$estimates[[control$args$group.column]])) %>%
                           filter(!!as.symbol(control$args$subgroup.column) == 1) %>% 
                           dplyr::select(-!!as.symbol(control$args$subgroup.column)) %>%
                           pivot_longer(
                             cols = -c(!!as.symbol(control$args$group.column), artificial.id), 
                             names_to = "feature", 
                             values_to = "value")%>% 
                           group_by(feature) %>% 
                           #mutate(overall.mean = mean(value, na.rm = TRUE)) %>% 
                           group_by(feature, !!as.symbol(control$args$group.column)) %>% 
                           mutate(group.mean = mean(value, na.rm = TRUE)) %>% 
                           ungroup() %>% dplyr::select(!!as.symbol(control$args$group.column), feature, group.mean))
        
        #print(head(hdata.pred))
      } else {
        suppressMessages(hdata.pred <- data %>% filter(!!as.symbol(control$args$group.column) %in% setdiff(data[[control$args$group.column]], control$estimates[[control$args$group.column]])) %>%
                           pivot_longer(
                             cols = -c(!!as.symbol(control$args$group.column), artificial.id), 
                             names_to = "feature", 
                             values_to = "value") %>% 
                           group_by(feature) %>% 
                           group_by(feature, !!as.symbol(control$args$group.column)) %>% 
                           mutate(group.mean = mean(value, na.rm = TRUE)) %>% 
                           ungroup() %>% dplyr::select(!!as.symbol(control$args$group.column), feature, group.mean))
        
      }
      
      overall.means <- control$estimates %>% select(feature, overall.mean) %>% drop_na() %>% distinct()
      suppressMessages(hdata.pred <- hdata.pred %>% left_join(overall.means, by = "feature"))
      
      suppressMessages(control$estimates <- control$estimates %>% full_join(hdata.pred) %>% distinct())
      
    } 
    
    
    if (control$args$subgroup.column %in% col.names){
      data <- data %>% dplyr::select(- !!as.symbol(control$args$subgroup.column)) 
    }
    
    suppressMessages(hdata <- data %>% 
                       pivot_longer(
                         cols = -c(!!as.symbol(control$args$group.column),artificial.id), 
                         names_to = "feature", 
                         values_to = "value") %>% 
                       left_join(control$estimates, by = c( control$args$group.column, "feature")) %>% 
                       mutate(h.value = value - (group.mean - overall.mean)) %>% 
                       dplyr::select(-!!as.symbol(args$group.column)))
    
    #print(head(hdata))
    suppressMessages(data <- hdata %>% 
                       dplyr::select(-c(value, overall.mean, group.mean)))
    
    suppressMessages(data <- data %>% 
                       pivot_wider(names_from = feature, values_from = h.value) %>% 
                       arrange(artificial.id) %>% 
                       dplyr::select(-artificial.id))
    
    
    return(data)
    
    
    col.names <- colnames(data)
    ids <- 1:nrow(data)
    data$artificial.id <- ids
    if (control$args$subgroup.column %in% col.names){
      data <- data %>% dplyr::select(-!!as.symbol(args$subgroup.column)) 
    }
    
    suppressMessages(hdata <- data %>% 
                       pivot_longer(
                         cols = -c(!!as.symbol(control$args$group.column),artificial.id), 
                         names_to = "feature", 
                         values_to = "value") %>% 
                       left_join(control$estimates, by = c( control$args$group.column, "feature")) %>% 
                       mutate(h.value = value - (group.mean - overall.mean)) %>% 
                       dplyr::select(-!!as.symbol(args$group.column)))
    
    
    suppressMessages(data <- hdata %>% 
                       dplyr::select(-c(value, overall.mean, group.mean)))
    
    suppressMessages(data <- data %>% 
                       pivot_wider(names_from = feature, values_from = h.value) %>% 
                       arrange(artificial.id) %>% 
                       dplyr::select(-artificial.id))
    
    
    return(data)
  }
  makePreprocWrapper(
    learner,
    train = trainfun,
    predict = predictfun,
    # par.set = makeParamSet(
    #   makeLogicalLearnerParam("group.column"),
    #   makeLogicalLearnerParam("subgroup.column")
    # ),
    par.vals = list(group.column = group.column, subgroup.column = subgroup.column)
  )
}
# ---------------------------------- 2: set up machine learning pipeline (nested resampling) -----------------------------------
# --------- 2.1: make a learning task ---------
task <- makeRegrTask(data = data_neuroimg_regr, target = "AGE")
learner <- makeLearner("regr.LiblineaRL2L1SVR")
# Define the tuning strategy
# Define the search space for hyperparameter tuning
param <- makeParamSet(
  makeNumericParam("cost",
                   lower = -6,
                   upper = 4,
                   trafo = function(x) 2^x)
)

ctrl <- makeTuneControlRandom(budget = 50) # budget should at least be 50, better around 250

# Define the inner resampling strategy
inner_resampling <- makeResampleDesc(
  "RepCV",
  folds = 5,
  reps = 5,
  stratify = FALSE
)

learner <-
  makePreprocWrapperCaret(
    learner,
    ppc.zv = TRUE,
    ppc.pca = FALSE,
    ppc.scale = TRUE,
    # ensures variables are scaled, is necessary for meaningful results in PCA
    ppc.center = TRUE
  )
# site correction 
learner <- 
  makePreprocWrapperGlobMeanDiff(
    learner, 
    group.column = "SITE_ID",
    subgroup.column = "HC_idx"
  )
# Define the nested cross-validation
learner <- makeTuneWrapper(
  learner = learner,
  resampling = inner_resampling,
  measures = list(rmse, mae, rsq),
  par.set = param,
  control = ctrl,
  show.info = FALSE
)

# Define the outer resampling strategy
outer_resampling <- makeResampleDesc(
  "RepCV",
  folds = 5,
  reps = 5,
  predict = "both"
)

# Set up future parallel backend
num_cpus <- detectCores()
cl <- makeCluster(num_cpus)
# Set seed for reproducibility
set.seed(1)

# Parallelize each resampling step
res <- mclapply(
  list(task),
  function(task) {
    resample(
      learner, 
      task, 
      outer_resampling, 
      measures = list(rmse, mae, rsq),
      extract = getTuneResult,
      show.info = TRUE,
      models = TRUE
    )
  },
  mc.cores = num_cpus
)

stopCluster(cl)


res<-res[[1]]
save(res,file = "/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_allsites_all_atlas_meanoffset_UniDistr_normHC_withgroup_offset.RData")

# Create a new task for the patient group
task_patients_regr_rem <- makeRegrTask(data = data_neuroimg_regr_rem, target = "AGE")
task_patients_regr_rem
task_patients_scz_all <- makeRegrTask(data = data_neuroimg_patients_scz_all, target = "AGE")
task_patients_scz_all
task_patients_scz_MUC <- makeRegrTask(data = data_neuroimg_patients_scz_MUC, target = "AGE")
task_patients_scz_MUC
task_patients_scz_VER <- makeRegrTask(data = data_neuroimg_patients_scz_VER, target = "AGE")
task_patients_scz_VER
task_patients_scz_OSL <- makeRegrTask(data = data_neuroimg_patients_scz_OSL, target = "AGE")
task_patients_scz_OSL
# Create a new task for the patient group
predictions_patients_regr_rem <- predict(res$model, task = task_patients_regr_rem)
predictions_patients_regr_rem
save(predictions_patients_regr_rem, file = "/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedremHC_withgroup_offset.RData")

predictions_patients_scz <- predict(res$model, task = task_patients_scz_all)
predictions_patients_scz
save(predictions_patients_scz, file = "/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedsczall_withgroup_offset.RData")

predictions_patients_scz_MUC <- predict(res$model, task = task_patients_scz_MUC)
predictions_patients_scz_MUC
save(predictions_patients_scz_MUC, file = "/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedsczMUC_withgroup_offset.RData")

predictions_patients_scz_VER <- predict(res$model, task = task_patients_scz_VER)
predictions_patients_scz_VER
save(predictions_patients_scz_VER, file = "/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedsczVER_withgroup_offset.RData")

predictions_patients_scz_OSL <- predict(res$model, task = data_neuroimg_patients_scz_OSL)
predictions_patients_scz_OSL
save(predictions_patients_scz_OSL, file = "/volume/mitnvp1_scratch/AK_Testing/ECNP/regression_age_all_atlas_meanoffset_UniDistr_appliedsczOSL_withgroup_offset.RData")
