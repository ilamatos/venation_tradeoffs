###################################################
# R-code LEAF VENATION FORM-FUNCTION TRADE-OFFS   #
###################################################

##### GRADIENT BOOSTED MODELS ################################

# Load packages -------------------
library(h2o)
library(tidyverse)
library(data.table)
library(iml)
library(ggpubr)
library(hrbrthemes)

# GBM Kleaf_max ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "Kleaf_max"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 1) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_Kleaf_max<-top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_major") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "VD_minor")
int_all<- rbind(int1$results, int2$results, int3$results)
# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction
inter3<-Partial$new(predictor.gbm, c("VD_medium", "VD_minor")) 

inter3<- inter3 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("GBM outputs/Fig5c_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM P50 ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_imput.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "P50"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 4) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_P50<-top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_major") 
int2<- Interaction$new(predictor.gbm, feature = "MST_minor")
int3<- Interaction$new(predictor.gbm, feature = "VD_minor")
int_all<- rbind(int1$results, int2$results, int3$results)
# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

h2o.shutdown() # shut down cluster


# GBM P88 ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_imput.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "P88"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 6) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_P88<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot
ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_minor") 
int2<- Interaction$new(predictor.gbm, feature = "ER_major")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

h2o.shutdown() # shut down cluster

# GBM SWP_M ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "SWP_M"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 10) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_SWP_M<-top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot
ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_major") 
int2<- Interaction$new(predictor.gbm, feature = "MST_major")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")
int_all<- rbind(int1$results, int2$results, int3$results)
# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

h2o.shutdown() # shut down cluster

# GBM SWP_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "SWP_L"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 12) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_SWP_L<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_medium") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_major")
int_all<- rbind(int1$results, int2$results, int3$results)
# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

h2o.shutdown() # shut down cluster

# GBM SWS_M ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "SWS_M"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 15) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_SWS_M<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot
ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_major") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "MST_minor")
int_all<- rbind(int1$results, int2$results, int3$results)
# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction
inter6<-Partial$new(predictor.gbm, c("MST_medium", "MST_major")) 

inter6<- inter6 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter6

ggsave(inter6,file=paste0("GBM outputs/Fig5f_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWS_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "SWS_L"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 19) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_SWS_L<-top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot
ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_major") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "MST_minor")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction
inter7<-Partial$new(predictor.gbm, c("MST_medium", "ER_major")) 

inter7<- inter7 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter7

ggsave(inter7,file=paste0("GBM outputs/Fig5g_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM Phe ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "Phe"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 21) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_Phe<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_medium") 
int2<- Interaction$new(predictor.gbm, feature = "VD_minor")
int3<- Interaction$new(predictor.gbm, feature = "VD_major")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction
inter8<-Partial$new(predictor.gbm, c("VD_major", "MST_medium")) 

inter8<- inter8 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter8

ggsave(inter8,file=paste0("GBM outputs/Fig5h_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM Dkleaf_mean ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "Dkleaf_mean"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 24) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mpd_Dkleaf_mean<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_medium") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_medium")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction
inter2<-Partial$new(predictor.gbm, c("VD_medium", "MST_medium")) 

inter2<- inter2 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2

ggsave(inter2,file=paste0("GBM outputs/Fig5b_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

inter4<-Partial$new(predictor.gbm, c("ER_major", "MST_medium")) 

inter4<- inter4 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter4

ggsave(inter4,file=paste0("GBM outputs/Fig5d_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM LMA ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "LMA"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 31) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_LMA<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification",
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_medium") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "VD_minor")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction

inter1<-Partial$new(predictor.gbm, c("ER_medium", "VD_minor")) 

inter1<- inter1 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter1

ggsave(inter1,file=paste0("GBM outputs/Fig5a",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM e_W ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "e_W"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 32) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_e_W<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_medium") 
int2<- Interaction$new(predictor.gbm, feature = "ER_major")
int3<- Interaction$new(predictor.gbm, feature = "VD_minor")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

h2o.shutdown() # shut down cluster

# GBM e_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "e_L"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 39) # set seed for reproducibility


# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_e_L<- top_mod

# Goodness of fit and R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_major") 
int2<- Interaction$new(predictor.gbm, feature = "VD_minor")
int3<- Interaction$new(predictor.gbm, feature = "MST_major")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

# coplot pairwise interaction
inter5<-Partial$new(predictor.gbm, c("ER_minor", "ER_major")) 

inter5<- inter5 %>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter5

ggsave(inter5,file=paste0("GBM outputs/Fig5e_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM ISI ------------------------

h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "ISI"
target<-dependent

# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]

# Fit  GBM using h2o.automl to do the grid and algorithm search 
v1 <- h2o.automl(y = dependent, 
                 x = independent, 
                 training_frame = v_train, 
                 validation_frame = v_test,
                 include_algos = c("GBM"),  # GBM model
                 max_runtime_secs = 30, # maximum running time
                 nfolds = 3, # number of K-folds cross-validations
                 seed = 40) # set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod
top_mod_ISI<- top_mod

# R-squared for training model
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
vimp<- as.data.frame(h2o.varimp(top_mod))
# save variable importance values
write_csv(vimp,paste0("GBM outputs/gbm_",target,"_vimp.csv")) 

# PDP and ICE Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_major")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_minor")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_major")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_minor")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_major")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_minor")

g_pdps<-ggarrange(p1+labs(title = NULL), 
                  p2+labs(title = NULL), 
                  p3+labs(title = NULL), 
                  p4+labs(title = NULL), 
                  p5+labs(title = NULL), 
                  p6+labs(title = NULL), 
                  p7+labs(title = NULL), 
                  p8+labs(title = NULL), 
                  p9+labs(title = NULL),
                  ncol = 3, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  align='hv',
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                           "(f)", "(g)", "(h)", "(i)"));g_pdps
# save PDP + ICE plots
ggsave(g_pdps,file=paste0("GBM outputs/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("GBM outputs/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)

# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
response <- as.data.frame(venation2%>% select(target)%>%na.omit)[,1]
# Create custom predict function that returns the predicted values as a
#    vector (probability of purchasing in our example)
pred <- function(model, newdata)  {
  results <- as.data.frame(h2o.predict(model, as.h2o(newdata)))
  return(results)
}
# create predictor object to pass to explainer functions
predictor.gbm <- Predictor$new(
  model = top_mod, 
  data = features, 
  y = response, 
  predict.fun = pred,
  class = "classification"
)

# Calculate H-statistic
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("GBM outputs/gbm_",target,"_interaction.csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_minor") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "VD_minor")
int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("GBM outputs/gbm_",target,"_interpair.csv"))

h2o.shutdown() # shut down cluster

# Plot Heatmaps [Figure 3] --------------------

# list of functional traits
res_list<- c("P50", "P88", "ISI", "SWP_M", "SWP_L",
             "SWS_M", "SWS_L", "Phe", "Dkleaf_mean","Kleaf_max", 
             "e_W", "e_L","LMA")

# list files with variable importance values
files<- list.files(path = "GBM outputs/",pattern = "*_vimp.csv", full.names = T )

# bind files with variable importance values
response<- list()
for (i in 1:length(res_list)){
  response[[i]]<- rep(res_list[i], times=10) 
}

dt<-lapply(files, read_csv) %>% 
  bind_rows()%>%
  mutate(response = unlist(response))
glimpse(dt)

# Plot heatmap importance
dt$response<- factor(dt$response, levels = res_list)
heat_vimp<-ggplot(dt, aes(response, variable, fill= scaled_importance)) + 
  geom_tile()+
  scale_fill_viridis(option = "D") +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom");heat_vimp

# list files with variable interaction values
files<- list.files(path = "GBM outputs/",pattern = "*_interaction.csv", full.names = T )
files

# bind files with variable interaction values
response<- list()
for (i in 1:length(res_list)){
  response[[i]]<- rep(res_list[i], times=10) 
}

dt<-lapply(files, read_csv) %>% 
  bind_rows()%>%
  mutate(response = unlist(response))
glimpse(dt)

# Plot heatmap pairwise interaction
dt$response<- factor(dt$response, levels = res_list)
heat_inter<-ggplot(dt, aes(response, .feature, fill= .interaction)) + 
  geom_tile()+
  scale_fill_viridis(option = "D") +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom");heat_inter

# plot total variance explained 

# get R-Squared values
r2_train[1] <- h2o.r2(top_mod_P50)
r2_train[2]<- h2o.r2(top_mod_P88)
r2_train[3]<- h2o.r2(top_mod_ISI)
r2_train[4]<- h2o.r2(top_mod_SWP_M)
r2_train[5]<- h2o.r2(top_mod_SWP_L)
r2_train[6]<- h2o.r2(top_mod_SWS_M)
r2_train[7]<- h2o.r2(top_mod_SWS_L)
r2_train[8]<- h2o.r2(top_mod_Phe)
r2_train[9]<- h2o.r2(top_mpd_Dkleaf_mean)
r2_train[10]<- h2o.r2(top_mod_Kleaf_max)
r2_train[11]<- h2o.r2(top_mod_e_W)
r2_train[12]<- h2o.r2(top_mod_e_L)
r2_train[13]<- h2o.r2(top_mod_LMA)

dt_var<-data.frame(res_list,r2_train)

var_exp<-ggplot(data=dt_var%>%mutate (r2_train = r2_train*100), aes(x=reorder(res_list, -r2_train), y=r2_train)) +
  geom_bar(stat="identity", fill = "gray50")+
  ylim(0,100)+
  ylab ("% Variance explained")+
  geom_hline(yintercept=50, linetype="dashed", color = "black", size =1)+
  theme_classic()+
  theme (axis.title.x=element_blank(),
         axis.text = element_text(size = 14),
         axis.title = element_text(size = 14),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)); var_exp

# create Figure 3
g_heat<-ggarrange(var_exp, heat_vimp, heat_inter,
                  ncol = 1, nrow = 3, common.legend = TRUE,
                  legend='bottom',
                  labels=c("(a)", "(b)", "(c)"));g_heat


# Test resilience across venation types -----------

dt<-read_csv("tables/TableS1_spp_list.csv")
glimpse(dt)
kruskal.test(Dkleaf_mean ~ venation_type, data = dt)
