
library(mlr)
library(tidyverse)
library(reshape2)
library(e1071)
library(caret)
library(ggplot2)
library(dplyr)
library(future)
library(parallel)


# ---------------------------------- 1: load & preprocess data  -----------------------------------

# load data
PIDs_matched_matched_pat <- read.csv("/usr/local/vipar/projects/project_17/PIDs_matchedHC_toPAT_manuscript1.csv", header = FALSE,col.names = c("PID"))

#PIDs_matched <- read.csv("/volume/HARMONY/ViPAR_analyses/ECNP_classif_ROI/PIDs_matchedHC_withoutMilan.csv")
PIDs_list_matchedtopat <- PIDs_matched_matched_pat$PID

data0 <- Clinical_MRI
#data0 <- read.csv("/volume/projects/JF_ECNP/Analysis/HC_SCZ_manuscript/Models_new/classif_noICA_all_atlas_HC_SCZ_check_mri.csv")
#matched HC to pat (only HC) data for later application
data_matched_topat <- data0[data0$PID %in% PIDs_list_matchedtopat, ]
dim(data_matched_topat)


unique(data_matched_topat$SITE_ID)


data_matched_topat_LMU<-filter(data_matched_topat, SITE_ID == 1)
data_matched_topat_Verona<-filter(data_matched_topat, SITE_ID == 8)
data_matched_topat_Turku<-filter(data_matched_topat, SITE_ID == 6)
data_matched_topat_Oslo<-filter(data_matched_topat, SITE_ID == 9)
dim(data_matched_topat_LMU)
dim(data_matched_topat_Verona)
dim(data_matched_topat_Turku)
dim(data_matched_topat_Oslo)

# For our analysis, we will need only the volumetric information.
data_matched_topat <- data_matched_topat %>% 
  select(c("PID","SITE_ID", "AGE", "SEX","Diagnosis","TIV","QC_grade", starts_with("schaefer200"),starts_with("aal3"),starts_with("hammers")))
data_pat=data0 %>% 
  select(c("PID","SITE_ID", "AGE", "SEX","Diagnosis","TIV","QC_grade", starts_with("schaefer200"),starts_with("aal3"),starts_with("hammers")))

data_pat<-filter(data_pat, Diagnosis == 1)
dim(data_pat)
data_pat <-na.omit(data_pat, TIV)
dim(data_pat)
unique(data_pat$SITE_ID)

data_pat_LMU<-filter(data_pat, SITE_ID == 1)
data_pat_Verona<-filter(data_pat, SITE_ID == 8)
data_pat_Turku<-filter(data_pat, SITE_ID == 6)
data_pat_Oslo<-filter(data_pat, SITE_ID == 9)
dim(data_pat_LMU)
dim(data_pat_Verona)
dim(data_pat_Turku)
dim(data_pat_Oslo)

data_HC_pat<-rbind(data_matched_topat,data_pat)
data_HC_pat<-subset(data_HC_pat, SITE_ID %in% c(1,6,8,9))

data_HC_pat_LMU<-filter(data_HC_pat, SITE_ID == 1)
data_HC_pat_Verona<-filter(data_HC_pat, SITE_ID == 8)
data_HC_pat_Turku<-filter(data_HC_pat, SITE_ID == 6)
data_HC_pat_Oslo<-filter(data_HC_pat, SITE_ID == 9)
dim(data_HC_pat_LMU)
dim(data_HC_pat_Verona)
dim(data_HC_pat_Turku)
dim(data_HC_pat_Oslo)


dim(data_HC_pat)
unique(data_HC_pat$Diagnosis)
# For our analysis, we will need only the volumetric information.
data_matched <- data_HC_pat %>% 
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

data_neuroimg <- cbind(data_matched$SITE_ID,data_matched$Diagnosis, data_neuroimg)
data_neuroimg$Diagnosis <- data_neuroimg$`data_matched$Diagnosis`
data_neuroimg$SITE_ID <- data_neuroimg$`data_matched$SITE_ID`
dim(data_neuroimg)
unique(data_neuroimg$Diagnosis)
unique(data_neuroimg$SITE_ID)

##HC remaining
data_remaining <- data0[!(data0$PID %in% PIDs_list_matchedtopat), ]
data_remaining<-filter(data_remaining, Diagnosis == 0)
data_remaining <- data_remaining %>% 
  select(c("PID","SITE_ID", "AGE", "SEX","Diagnosis","TIV","QC_grade", starts_with("schaefer200"),starts_with("aal3"),starts_with("hammers")))
unique(data_remaining$SITE_ID)
unique(data_remaining$Diagnosis)
data_remaining <- data_remaining %>% mutate_if(is.character, as.numeric)
dim(data_remaining)
data_remaining <-na.omit(data_remaining , cols = c(starts_with("schaefer200")))
dim(data_remaining)

TIV_remaining <- data_remaining$TIV
data_neuroimg_rem <- data_remaining  %>% 
  select(c(starts_with("schaefer200_gmv"),starts_with("aal3"),starts_with("hammers_gmv")))
data_neuroimg_rem<-data.frame((data_neuroimg_rem/TIV_remaining)*1000)
min(data_neuroimg_rem)
max(data_neuroimg_rem)

data_neuroimg_rem <- cbind(data_remaining$SITE_ID,data_remaining$Diagnosis, data_neuroimg_rem)
data_neuroimg_rem$Diagnosis <- data_neuroimg_rem$`data_remaining$Diagnosis`
data_neuroimg_rem$SITE_ID <- data_neuroimg_rem$`data_remaining$SITE_ID`
dim(data_neuroimg_rem)
unique(data_neuroimg_rem$Diagnosis)
unique(data_neuroimg_rem$SITE_ID)

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
                       mutate(h.value = value - (overall.mean - group.mean)) %>% 
                       dplyr::select(-!!as.symbol(args$group.column)))
    
    
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
                       mutate(h.value = value - (overall.mean - group.mean)) %>% 
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
  
data_neuroimg$Diagnosis <- as.factor(data_neuroimg$Diagnosis)

data_neuroimg <- subset(data_neuroimg, SITE_ID %in% c(1,8,9)) #only data for munich, oslo and verona
dim(data_neuroimg)
data_neuroimg$HC_idx <- ifelse(data_neuroimg$Diagnosis == 0, 1, 0)
dim(data_neuroimg)
data_neuroimg_rem$HC_idx <- ifelse(data_neuroimg_rem$Diagnosis == 0, 1, 0)
dim(data_neuroimg_rem)
############## Separating for model application
data_neuroimg_rest_munich_oslo=subset(data_neuroimg, SITE_ID %in% c(1,9))
data_neuroimg_verona=filter(data_neuroimg, SITE_ID==8)
data_neuroimg_rest_munich_verona=subset(data_neuroimg, SITE_ID %in% c(1,8))
data_neuroimg_oslo=filter(data_neuroimg, SITE_ID==9)
data_neuroimg_rest_oslo_verona=subset(data_neuroimg, SITE_ID %in% c(9,8))
data_neuroimg_munich=filter(data_neuroimg, SITE_ID==1)
data_neuroimg_rem_HC_int<- subset(data_neuroimg_rem, SITE_ID %in% c(1,8,9))
dim(data_neuroimg_rem_HC_int)
data_neuroimg_rem_HC_ext<- subset(data_neuroimg_rem, SITE_ID %in% c(5,3,2,6))
dim(data_neuroimg_rem_HC_ext)


data_neuroimg_classif <- data_neuroimg %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))

##model trained on munich and oslo and applied to verona
data_neuroimg_classif_rest_munich_oslo <- data_neuroimg_rest_munich_oslo %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_rest_munich_oslo)
data_neuroimg_classif_verona <- data_neuroimg_verona %>% 
  select(c("SITE_ID","HC_idx","Diagnosis",starts_with("aal3")))
dim(data_neuroimg_classif_verona)
##model trained on munich and verona and applied to oslo
data_neuroimg_classif_rest_munich_verona <- data_neuroimg_rest_munich_verona %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_rest_munich_verona)
data_neuroimg_classif_oslo <- data_neuroimg_oslo %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_oslo)

data_neuroimg_classif_rest_oslo_verona <- data_neuroimg_rest_oslo_verona %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_rest_oslo_verona)
data_neuroimg_classif_munich <- data_neuroimg_munich %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_munich)

data_neuroimg_classif_rem_HC_int <- data_neuroimg_rem_HC_int %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_rem_HC_int)
data_neuroimg_classif_rem_HC_ext <- data_neuroimg_rem_HC_ext %>% 
  select(c("SITE_ID","Diagnosis","HC_idx",starts_with("aal3")))
dim(data_neuroimg_classif_rem_HC_ext)


# ---------------------------------- 2: set up machine learning pipeline (nested resampling) -----------------------------------
# --------- 2.1: make a learning task ---------
task_neuroimg_classif <-
  makeClassifTask(
    "MUC_neuroimg_task",
    data = data_neuroimg_classif,
    target = "Diagnosis",
    positive = 1
  )
task_neuroimg_classif_rest_munich_verona <-
  makeClassifTask(
    "MUC_VER_neuroimg_task",
    data = data_neuroimg_classif_rest_munich_verona,
    target = "Diagnosis",
    positive = 1
  )
task_neuroimg_classif_rest_munich_oslo <-
  makeClassifTask(
    "MUC_OSL_neuroimg_task",
    data = data_neuroimg_classif_rest_munich_oslo,
    target = "Diagnosis",
    positive = 1
  )
task_neuroimg_classif_rest_oslo_verona <-
  makeClassifTask(
    "OSL_VER_neuroimg_task",
    data = data_neuroimg_classif_rest_oslo_verona,
    target = "Diagnosis",
    positive = 1
  )


task_neuroimg_classif_verona <-
  makeClassifTask(
    "VER_neuroimg_task",
    data = data_neuroimg_classif_verona,
    target = "Diagnosis",
    positive = 1
  )

task_neuroimg_classif_oslo <-
  makeClassifTask(
    "OSL_neuroimg_task",
    data = data_neuroimg_classif_oslo,
    target = "Diagnosis",
    positive = 1
  )

task_neuroimg_classif_munich <-
  makeClassifTask(
    "MUC_neuroimg_task",
    data = data_neuroimg_classif_munich,
    target = "Diagnosis",
    positive = 1
  )

task_neuroimg_classif_HCrem_int <-
  makeClassifTask(
    "HCrem_int_neuroimg_task",
    data = data_neuroimg_classif_rem_HC_int,
    target = "Diagnosis",
    positive = 0
  )

task_neuroimg_classif_HCrem_ext <-
  makeClassifTask(
    "HCrem_ext_neuroimg_task",
    data = data_neuroimg_classif_rem_HC_ext,
    target = "Diagnosis",
    positive = 0
  )

# let's check out the ML task we've created
task_neuroimg_classif
task_neuroimg_classif_rest_munich_oslo
task_neuroimg_classif_rest_munich_verona
task_neuroimg_classif_rest_oslo_verona
task_neuroimg_classif_verona
task_neuroimg_classif_oslo
task_neuroimg_classif_munich
task_neuroimg_classif_HCrem_int
task_neuroimg_classif_HCrem_ext

# --------- 2.2: construct a learner ---------
lrn_classif <- makeLearner("classif.LiblineaRL1L2SVC")

# define hyperparameters of SVM + search space we are going to cover

ps_classif <- makeParamSet(
  makeNumericParam(
    "cost",
    lower = -8,
    upper = 8,
    trafo = function(x) 2^x
  )
)

# define the optimization algorithm (aka tuning method)
ctrl_classif <-
  makeTuneControlRandom(budget = 150) # budget should at least be 50, better around 250

# define an evaluation method, (i.e., a resampling strategy)
inner_classif <-
  makeResampleDesc("RepCV",
                   folds = 5,
                   reps = 5,
                   stratify = TRUE)

# --------- 2.3: feature preprocessing ---------
# remove 0 variance feats, imputation, scaling, PCA 
# Define the scaling range
scaling_range <- c(-1, 1)

lrn_classif <-
  makePreprocWrapperCaret(
    lrn_classif,
    ppc.zv = TRUE,
    ppc.pca = FALSE,
    # ensures variables are scaled, is necessary for meaningful results in PCA
    ppc.scale = TRUE,
    ppc.center = TRUE
    )


# site correction 
lrn_classif <- 
  makePreprocWrapperGlobMeanDiff(
    lrn_classif, 
    group.column = "SITE_ID",
    subgroup.column = "HC_idx"
  )

# --------- 2.4: initialize inner resampling scheme ---------
lrn_classif <- makeTuneWrapper(
  lrn_classif,
  resampling = inner_classif,
  par.set = ps_classif,
  measures = list(bac, tpr, tnr),
  # we use BAC as our perfomance measure, and TPR/TNR are returned in addition
  control = ctrl_classif,
  show.info = FALSE
)

# --------- 2.5: initialize outer resampling scheme ---------
outer_classif <-
  makeResampleDesc(
    "RepCV",
    folds = 5,
    reps = 5,
    predict = "both",
    # save results from both train/test folds
    stratify = TRUE
  )
# run the analysis
# run the analysis
# Set the number of CPUs you want to use
num_cpus <- detectCores()
cl <- makeCluster(num_cpus)
# Set seed for reproducibility
set.seed(1)

# Parallelize each resampling step
resampling_neuroimg <- mclapply(
  list(task_neuroimg_classif),
  function(task) {
    resample(
      lrn_classif,
      task,
      resampling = outer_classif,
      measures = list(bac, tpr, tnr),
      extract = getTuneResult,
      show.info = TRUE,
      models = TRUE
    )
  },
  mc.cores = num_cpus
)
resampling_neuroimg<-resampling_neuroimg[[1]]

  
save(resampling_neuroimg,file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_MUC_OSL_VER.RData")


set.seed(1)

resampling_neuroimg_rest_munich_oslo <- mclapply(
  list(task_neuroimg_classif_rest_munich_oslo),
  function(task) {
    resample(
      lrn_classif,
      task,
      resampling = outer_classif,
      measures = list(bac, tpr, tnr),
      extract = getTuneResult,
      show.info = TRUE,
      models = TRUE
    )
  },
  mc.cores = num_cpus
)
resampling_neuroimg_rest_munich_oslo<-resampling_neuroimg_rest_munich_oslo[[1]]

save(resampling_neuroimg_rest_munich_oslo,file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_MUC_OSL.RData")

set.seed(1)
resampling_neuroimg_rest_munich_verona <- mclapply(
  list(task_neuroimg_classif_rest_munich_verona),
  function(task) {
    resample(
      lrn_classif,
      task,
      resampling = outer_classif,
      measures = list(bac, tpr, tnr),
      extract = getTuneResult,
      show.info = TRUE,
      models = TRUE
    )
  },
  mc.cores = num_cpus
)
resampling_neuroimg_rest_munich_verona<-resampling_neuroimg_rest_munich_verona[[1]]

save(resampling_neuroimg_rest_munich_verona,file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_MUC_VER.RData")

set.seed(1)

resampling_neuroimg_rest_oslo_verona <- mclapply(
  list(task_neuroimg_classif_rest_oslo_verona),
  function(task) {
    resample(
      lrn_classif,
      task,
      resampling = outer_classif,
      measures = list(bac, tpr, tnr),
      extract = getTuneResult,
      show.info = TRUE,
      models = TRUE
    )
  },
  mc.cores = num_cpus
)
resampling_neuroimg_rest_oslo_verona<-resampling_neuroimg_rest_oslo_verona[[1]]
stopCluster(cl)

save(resampling_neuroimg_rest_oslo_verona,file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_OSL_VER.RData")

##model application to verona
model_munich_oslo<-resampling_neuroimg_rest_munich_oslo$model
preds_verona=predict(model_munich_oslo,task=task_neuroimg_classif_verona)
save(preds_verona, file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_appliedonVerona.RData")
# Calculate the overall balanced accuracy
preds=preds_verona
overall_bacc <- 0
overall_sen <- 0
overall_spe <- 0
balanced_acc_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_sen_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_spe_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector

for (i in 1:length(preds)) {
  predicted_labels <- preds[[i]]$data$response
  true_labels <- preds[[i]]$data$truth
  conf_matrix <- confusionMatrix(predicted_labels, true_labels)
  
  # Calculate sensitivity
  TP <- conf_matrix$table[2, 2]
  FN <- conf_matrix$table[2, 1]
  sensitivity <- TP / (TP + FN)
  
  # Calculate specificity
  TN <- conf_matrix$table[1, 1]
  FP <- conf_matrix$table[1, 2]
  specificity <- TN / (TN + FP)
  
  # Calculate balanced accuracy as the mean of sensitivity and specificity
  balanced_acc <- (sensitivity + specificity) / 2
  overall_sen<-overall_sen+sensitivity
  overall_spe<-overall_spe+specificity
  overall_bacc <- overall_bacc + balanced_acc
  balanced_acc_all[i]<-balanced_acc
  overall_sen_all[i]<-sensitivity
  overall_spe_all[i]<-specificity
}
balanced_acc_all
overall_sen_all
overall_spe_all
overall_bacc <- overall_bacc / length(preds)
overall_sen <- overall_sen / length(preds)
overall_spe <- overall_spe / length(preds)
overall_bacc
overall_sen
overall_spe
mean(balanced_acc_all)
mean(overall_sen_all)
mean(overall_spe_all)

################model application to oslo
model_munich_verona<-resampling_neuroimg_rest_munich_verona$model
preds_oslo=predict(model_munich_verona,task=task_neuroimg_classif_oslo)
save(preds_oslo, file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_appliedonOslo.RData")

preds=preds_oslo
overall_bacc <- 0
overall_sen <- 0
overall_spe <- 0
balanced_acc_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_sen_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_spe_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector

for (i in 1:length(preds)) {
  predicted_labels <- preds[[i]]$data$response
  true_labels <- preds[[i]]$data$truth
  conf_matrix <- confusionMatrix(predicted_labels, true_labels)
  
  # Calculate sensitivity
  TP <- conf_matrix$table[2, 2]
  FN <- conf_matrix$table[2, 1]
  sensitivity <- TP / (TP + FN)
  
  # Calculate specificity
  TN <- conf_matrix$table[1, 1]
  FP <- conf_matrix$table[1, 2]
  specificity <- TN / (TN + FP)
  
  # Calculate balanced accuracy as the mean of sensitivity and specificity
  balanced_acc <- (sensitivity + specificity) / 2
  overall_sen<-overall_sen+sensitivity
  overall_spe<-overall_spe+specificity
  overall_bacc <- overall_bacc + balanced_acc
  balanced_acc_all[i]<-balanced_acc
  overall_sen_all[i]<-sensitivity
  overall_spe_all[i]<-specificity
}
balanced_acc_all
overall_sen_all
overall_spe_all
overall_bacc <- overall_bacc / length(preds)
overall_sen <- overall_sen / length(preds)
overall_spe <- overall_spe / length(preds)
overall_bacc
overall_sen
overall_spe
mean(balanced_acc_all)
mean(overall_sen_all)
mean(overall_spe_all)



#####application to munich
model_oslo_verona<-resampling_neuroimg_rest_oslo_verona$models
preds_munich=predict(model_oslo_verona,task=task_neuroimg_classif_munich)
save(preds_munich, file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_appliedonMunich.RData")

preds=preds_munich
overall_bacc <-0
overall_sen <- 0
overall_spe <- 0
balanced_acc_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_sen_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_spe_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector

for (i in 1:length(preds)) {
  predicted_labels <- preds[[i]]$data$response
  true_labels <- preds[[i]]$data$truth
  conf_matrix <- confusionMatrix(predicted_labels, true_labels)
  
  # Calculate sensitivity
  TP <- conf_matrix$table[2, 2]
  FN <- conf_matrix$table[2, 1]
  sensitivity <- TP / (TP + FN)
  
  # Calculate specificity
  TN <- conf_matrix$table[1, 1]
  FP <- conf_matrix$table[1, 2]
  specificity <- TN / (TN + FP)
  
  # Calculate balanced accuracy as the mean of sensitivity and specificity
  balanced_acc <- (sensitivity + specificity) / 2
  overall_sen<-overall_sen+sensitivity
  overall_spe<-overall_spe+specificity
  overall_bacc <- overall_bacc + balanced_acc
  balanced_acc_all[i]<-balanced_acc
  overall_sen_all[i]<-sensitivity
  overall_spe_all[i]<-specificity
}
balanced_acc_all
overall_sen_all
overall_spe_all
overall_bacc <- overall_bacc / length(preds)
overall_sen <- overall_sen / length(preds)
overall_spe <- overall_spe / length(preds)
overall_bacc
overall_sen
overall_spe
mean(balanced_acc_all)
mean(overall_sen_all)
mean(overall_spe_all)

####model application to remaining HC
model_munich_verona_osl<-resampling_neuroimg$model
preds_HC_int=predict(model_munich_verona_osl,task=task_neuroimg_classif_HCrem_int)
save(preds_HC_int, file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_appliedonHC_int.RData")

preds=preds_HC_int
overall_bacc <- 0
overall_sen <- 0
overall_spe <- 0
balanced_acc_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_sen_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_spe_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector

for (i in 1:length(preds)) {
  predicted_labels <- preds[[i]]$data$response
  true_labels <- preds[[i]]$data$truth
  conf_matrix <- confusionMatrix(predicted_labels, true_labels)
  
  # Calculate sensitivity
  TP <- conf_matrix$table[2, 2]
  FN <- conf_matrix$table[2, 1]
  sensitivity <- TP / (TP + FN)
  
  # Calculate specificity
  TN <- conf_matrix$table[1, 1]
  FP <- conf_matrix$table[1, 2]
  specificity <- TN / (TN + FP)
  
  # Calculate balanced accuracy as the mean of sensitivity and specificity
  balanced_acc <- (sensitivity + specificity) / 2
  overall_sen<-overall_sen+sensitivity
  overall_spe<-overall_spe+specificity
  overall_bacc <- overall_bacc + balanced_acc
  balanced_acc_all[i]<-balanced_acc
  overall_sen_all[i]<-sensitivity
  overall_spe_all[i]<-specificity
}
balanced_acc_all
overall_sen_all
overall_spe_all
overall_bacc <- overall_bacc / length(preds)
overall_sen <- overall_sen / length(preds)
overall_spe <- overall_spe / length(preds)
overall_bacc
overall_sen
overall_spe
mean(balanced_acc_all)
mean(overall_sen_all)
mean(overall_spe_all)


####model application to remaining HC
model_munich_verona_osl<-resampling_neuroimg$model
preds_HC_ext=predict(model_munich_verona_osl,task=task_neuroimg_classif_HCrem_ext)
save(preds_HC_ext, file = "/usr/local/vipar/projects/project_17/ECNP_classif_HCSCZ_aal3_globaloffset_appliedonHC_ext.RData")

preds=preds_HC_int
overall_bacc <- 0
overall_sen <- 0
overall_spe <- 0
balanced_acc_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_sen_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector
overall_spe_all <- vector("numeric", length = length(preds)) # Initialize as a numeric vector

for (i in 1:length(preds)) {
  predicted_labels <- preds[[i]]$data$response
  true_labels <- preds[[i]]$data$truth
  conf_matrix <- confusionMatrix(predicted_labels, true_labels)
  
  # Calculate sensitivity
  TP <- conf_matrix$table[2, 2]
  FN <- conf_matrix$table[2, 1]
  sensitivity <- TP / (TP + FN)
  
  # Calculate specificity
  TN <- conf_matrix$table[1, 1]
  FP <- conf_matrix$table[1, 2]
  specificity <- TN / (TN + FP)
  
  # Calculate balanced accuracy as the mean of sensitivity and specificity
  balanced_acc <- (sensitivity + specificity) / 2
  overall_sen<-overall_sen+sensitivity
  overall_spe<-overall_spe+specificity
  overall_bacc <- overall_bacc + balanced_acc
  balanced_acc_all[i]<-balanced_acc
  overall_sen_all[i]<-sensitivity
  overall_spe_all[i]<-specificity
}
balanced_acc_all
overall_sen_all
overall_spe_all
overall_bacc <- overall_bacc / length(preds)
overall_sen <- overall_sen / length(preds)
overall_spe <- overall_spe / length(preds)
overall_bacc
overall_sen
overall_spe
mean(balanced_acc_all)
mean(overall_sen_all)
mean(overall_spe_all)



