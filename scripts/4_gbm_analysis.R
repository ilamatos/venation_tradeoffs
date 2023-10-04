###########################################################
# R-code LEAF VENATION ARCHITECTURE-FUNCTION TRADE-OFFS   #
###########################################################

##### GRADIENT BOOSTED MODELS ################################

# Load packages -------------------
library(h2o)
library(tidyverse)
library(data.table)
library(iml)
library(ggpubr)
library(hrbrthemes)

# GBM SCALED RMIN ###############################################
# GBM Kleaf_max ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 4) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)

# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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

Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_large") 
int2<- Interaction$new(predictor.gbm, feature = "VD_small")
int3<- Interaction$new(predictor.gbm, feature = "ER_small")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)
# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependece co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_medium", "ER_small")) 
inter2<-Partial$new(predictor.gbm, c("VD_small", "VD_large")) 
inter3<-Partial$new(predictor.gbm, c("VD_small", "ER_small"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


h2o.shutdown() # shut down cluster

# GBM P50 ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_imput_trans.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 6) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_large", "VD_medium")) 
inter2<-Partial$new(predictor.gbm, c("MST_large", "MST_medium")) 
inter3<-Partial$new(predictor.gbm, c("VD_medium", "MST_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM P88 ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_imput_trans.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "VD_large")
int3<- Interaction$new(predictor.gbm, feature = "ER_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_small", "VD_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_small", "ER_large")) 
inter3<-Partial$new(predictor.gbm, c("VD_medium", "VD_large"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


h2o.shutdown() # shut down cluster


# GBM SWP_M ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 16)# set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_medium") 
int2<- Interaction$new(predictor.gbm, feature = "VD_large")
int3<- Interaction$new(predictor.gbm, feature = "VD_small")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("ER_medium", "VD_large")) 
inter2<-Partial$new(predictor.gbm, c("ER_large", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("ER_small", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWP_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 14)# set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_small") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_medium", "ER_small")) 
inter2<-Partial$new(predictor.gbm, c("VD_large", "ER_small")) 
inter3<-Partial$new(predictor.gbm, c("VD_medium", "MST_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWS_M ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_large") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_small")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_large", "ER_small")) 
inter2<-Partial$new(predictor.gbm, c("MST_medium", "MST_large")) 
inter3<-Partial$new(predictor.gbm, c("MST_small", "MST_large"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWS_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_large") 
int2<- Interaction$new(predictor.gbm, feature = "ER_small")
int3<- Interaction$new(predictor.gbm, feature = "MST_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_large", "MST_medium")) 
inter2<-Partial$new(predictor.gbm, c("ER_small", "VD_large")) 
inter3<-Partial$new(predictor.gbm, c("ER_small", "MST_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM Phe ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 25) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_medium") 
int2<- Interaction$new(predictor.gbm, feature = "VD_large")
int3<- Interaction$new(predictor.gbm, feature = "VD_small")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_large", "ER_medium")) 
inter2<-Partial$new(predictor.gbm, c("MST_medium", "VD_large")) 
inter3<-Partial$new(predictor.gbm, c("ER_large", "VD_small"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM Dkleaf_mean ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "ER_medium")
int3<- Interaction$new(predictor.gbm, feature = "MST_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_small", "MST_medium")) 
inter2<-Partial$new(predictor.gbm, c("VD_medium", "MST_medium")) 
inter3<-Partial$new(predictor.gbm, c("ER_medium", "MST_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM LMA ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 39) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_large", "ER_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_large", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("MST_small", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM e_W ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 32) #  set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_medium") 
int2<- Interaction$new(predictor.gbm, feature = "MST_large")
int3<- Interaction$new(predictor.gbm, feature = "VD_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_medium", "MST_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_large", "VD_large")) 
inter3<-Partial$new(predictor.gbm, c("MST_large", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM e_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_large") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_large", "VD_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_large", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("MST_small", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM ISI ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod


# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_scaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_scaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_scaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_scaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_scaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_small") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_medium", "ER_medium")) 
inter2<-Partial$new(predictor.gbm, c("MST_medium", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("ER_medium", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_scaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_scaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_scaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_scaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_scaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_scaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM UNSCALED RMIN ###############################################
# GBM Kleaf_max ------------------------

h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 23) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)

# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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

Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)
# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_medium") 
int2<- Interaction$new(predictor.gbm, feature = "VD_medium")
int3<- Interaction$new(predictor.gbm, feature = "MST_large")

int_all<- rbind(int1$results, int2$results, int3$results)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_scaled/gmb_interpair_",target,".csv"))

# Partial dependece co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_medium", "ER_medium")) 
inter2<-Partial$new(predictor.gbm, c("ER_medium", "MST_large")) 
inter3<-Partial$new(predictor.gbm, c("VD_large", "MST_large"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


h2o.shutdown() # shut down cluster

# GBM P50 ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_imput_trans_log.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 7) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "ER_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_small", "MST_medium")) 
inter2<-Partial$new(predictor.gbm, c("VD_medium", "MST_medium")) 
inter3<-Partial$new(predictor.gbm, c("VD_medium", "ER_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


h2o.shutdown() # shut down cluster

# GBM P88 ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_imput_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "VD_large")
int3<- Interaction$new(predictor.gbm, feature = "ER_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_small", "VD_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_small", "ER_large")) 
inter3<-Partial$new(predictor.gbm, c("VD_medium", "VD_large"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWP_M ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 16) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_medium") 
int2<- Interaction$new(predictor.gbm, feature = "ER_small")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_medium", "MST_medium")) 
inter2<-Partial$new(predictor.gbm, c("VD_large", "ER_small")) 
inter3<-Partial$new(predictor.gbm, c("VD_small", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWP_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 16) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "ER_small") 
int2<- Interaction$new(predictor.gbm, feature = "VD_small")
int3<- Interaction$new(predictor.gbm, feature = "MST_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_small", "ER_small")) 
inter2<-Partial$new(predictor.gbm, c("VD_large", "MST_medium")) 
inter3<-Partial$new(predictor.gbm, c("MST_large", "ER_small"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWS_M ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 18) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_large") 
int2<- Interaction$new(predictor.gbm, feature = "VD_small")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_medium", "VD_large")) 
inter2<-Partial$new(predictor.gbm, c("ER_medium", "VD_large")) 
inter3<-Partial$new(predictor.gbm, c("ER_medium", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM SWS_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 24) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_medium") 
int2<- Interaction$new(predictor.gbm, feature = "MST_large")
int3<- Interaction$new(predictor.gbm, feature = "VD_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_small", "VD_medium")) 
inter2<-Partial$new(predictor.gbm, c("ER_medium", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("MST_large", "VD_large"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM Phe ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 25) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_large") 
int2<- Interaction$new(predictor.gbm, feature = "ER_large")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_medium", "ER_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_small", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("VD_large", "ER_large"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM Dkleaf_mean ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 27) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "MST_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_large", "MST_medium")) 
inter2<-Partial$new(predictor.gbm, c("ER_large", "VD_small")) 
inter3<-Partial$new(predictor.gbm, c("MST_medium", "VD_small"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM LMA ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 40) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_scaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_small") 
int2<- Interaction$new(predictor.gbm, feature = "MST_medium")
int3<- Interaction$new(predictor.gbm, feature = "MST_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_large", "MST_large")) 
inter2<-Partial$new(predictor.gbm, c("ER_large", "MST_medium")) 
inter3<-Partial$new(predictor.gbm, c("MST_large", "VD_small"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM e_W ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 33) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_medium") 
int2<- Interaction$new(predictor.gbm, feature = "MST_small")
int3<- Interaction$new(predictor.gbm, feature = "VD_large")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("VD_small", "VD_large")) 
inter2<-Partial$new(predictor.gbm, c("MST_small", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("VD_small", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)
h2o.shutdown() # shut down cluster

# GBM e_L ------------------------
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"

venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 40) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "MST_small") 
int2<- Interaction$new(predictor.gbm, feature = "ER_small")
int3<- Interaction$new(predictor.gbm, feature = "VD_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("ER_small", "MST_small")) 
inter2<-Partial$new(predictor.gbm, c("ER_large", "MST_small")) 
inter3<-Partial$new(predictor.gbm, c("ER_large", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster

# GBM ISI ------------------------

h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 42) # set seed for reproducibility

# See all fitted models
v1 %>% summary
h2o.get_leaderboard(v1)

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod


# Goodness of fit and R-squared for training model
h2o.performance(top_mod, newdata = v_train) # train GOF
r2_train <- h2o.r2(top_mod)
r2_train # R-squared train

# Variable importance
h2o.varimp(top_mod)
vimp<- as.data.frame(h2o.varimp(top_mod))
glimpse(vimp)
# save variable importance values
write_csv(vimp,paste0("figures/GBM_outputs_unscaled/gbm_",target,"_vimp.csv")) 
# plot variables importance
h2o.varimp_plot(top_mod) 
# variable importance over all models
va_plot <- h2o.varimp_heatmap(v1)
va_plot

png(filename= paste0("figures/GBM_outputs_unscaled/va_plot_",target,".png"), width=650, height=650)
va_plot
dev.off()

# PDP and ICE Plots
# Partial Dependence Plots (PDP)
# Individual Conditional Expectation (ICE) Plots
p1<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_large")
p2<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_medium")
p3<-h2o.ice_plot(top_mod, newdata = v_train, column = "VD_small")
p4<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_large")
p5<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_medium")
p6<-h2o.ice_plot(top_mod, newdata = v_train, column = "ER_small")
p7<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_large")
p8<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_medium")
p9<-h2o.ice_plot(top_mod, newdata = v_train, column = "MST_small")

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
ggsave(g_pdps,file=paste0("figures/GBM_outputs_unscaled/PDP_ICE_",target,".png"),width=24,height=30, units = "cm", dpi = 300)

# SHAP (SHapley Additive exPlanation) values 
shap<-as.data.frame(h2o.predict_contributions(top_mod, v_train)) # compute SHAP values
# save SHAP values
write_csv(shap, paste0("figures/GBM_outputs_unscaled/gbm_",target,"_shap.csv"))
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot

ggsave(shap_plot,file=paste0("figures/GBM_outputs_unscaled/SHAP_",target,".png"),width=12,height=12, units = "cm", dpi = 300)


# Pairwise interactions strength using iml
# Create a data frame with just the predictor variables
features <- as.data.frame(venation) %>% dplyr::select(all_of(independent), dependent)%>%na.omit%>% select(-dependent)
# Create a vector with the actual responses
venation2 <- venation %>% 
  as.data.table() %>% 
  mutate(pred = as.data.table(h2o.predict(top_mod, newdata=venation))$predict)
plot(x= as.data.frame(venation2%>% select(target))[,1], y = venation2$pred) # predicted x observed
abline(0,1)
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
#Interaction$new(predictor.gbm) %>% plot()  
interact<- Interaction$new(predictor.gbm)
# Save H values
write_csv(interact$results, paste0("figures/GBM_outputs_unscaled/gmb_interaction_",target,".csv"))

head(arrange(interact$results, desc(.interaction)), n = 6)

# Measure the two-way interactions of the 3 first interactive variables (3 first variables with higher H-statistic) with the target variable
int1<-Interaction$new(predictor.gbm, feature = "VD_medium") 
int2<- Interaction$new(predictor.gbm, feature = "MST_small")
int3<- Interaction$new(predictor.gbm, feature = "MST_medium")

int_all<- rbind(int1$results, int2$results, int3$results)
glimpse(int_all)

# Save H values of pairwise interactions
write_csv(int_all, paste0("figures/GBM_outputs_unscaled/gmb_interpair_",target,".csv"))

# Partial dependence co-plots for the 3 strongest pairwise interactions
# what are the 3 strongest pairwise interactions?
head(arrange(int_all, desc(.interaction)), n = 6)

inter1<-Partial$new(predictor.gbm, c("MST_small", "VD_medium")) 
inter2<-Partial$new(predictor.gbm, c("ER_medium", "VD_medium")) 
inter3<-Partial$new(predictor.gbm, c("MST_large", "VD_medium"))

# Save PDP values
write_csv(inter1$results, paste0("figures/GBM_outputs_unscaled/inter1_",target,".csv"))
write_csv(inter2$results, paste0("figures/GBM_outputs_unscaled/inter2_",target,".csv"))
write_csv(inter3$results, paste0("figures/GBM_outputs_unscaled/inter3_",target,".csv"))


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

ggsave(inter1,file=paste0("figures/GBM_outputs_unscaled/inter1_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter2<- inter2%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter2
ggsave(inter2,file=paste0("figures/GBM_outputs_unscaled/inter2_",target,".png"),width=12,height=8, units = "cm", dpi = 300)


inter3<- inter3%>% plot()+
  scale_fill_viridis_c(option='D')+
  theme_bw()+
  labs(size =2)+ 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); inter3

ggsave(inter3,file=paste0("figures/GBM_outputs_unscaled/inter3_",target,".png"),width=12,height=8, units = "cm", dpi = 300)

h2o.shutdown() # shut down cluster


# Plot total variance explained [Figure 4a]-----------------

dt_var<- read_csv("data/variance_explained.csv")
glimpse(dt_var)

var_exp<-ggplot(data=dt_var%>%mutate (r2_train = r2_train*100), aes(x=reorder(`Response variable`, -r2_train), y=r2_train, col = status, fill = `Response variable` )) +
  geom_bar(stat="identity", position="dodge")+
  ylim(0,100)+
  ylab ("% Variance explained")+
  theme_classic()+
  scale_color_manual(values = c("scaled" = "black", "unscaled" = "black"))+
  scale_fill_manual(values = c("P50"= "#0072B2",
                               "P88"= "#0072B2",
                               "ISI"="#0072B2",
                               "Dkleaf_mean"= "#E69F00",
                               "e_W" = "#CC79A7",
                               "e_L" = "#CC79A7",
                               "SWS_M" = "#009E73",
                               "SWS_L" = "#009E73",
                               "SWP_M" = "#009E73",
                               "SWP_L" = "#009E73",
                               "Phe" = "#009E73",
                               "Kleaf_max" = "#D55E00",
                               "LMA" =  "#56B4E9"))+
  theme (axis.title.x=element_blank(),
         axis.text = element_text(size = 14),
         axis.title = element_text(size = 14),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)); var_exp

ggsave(var_exp,file="figures/Figure 4a.png",width=28,height=18, units = "cm", dpi = 300)

# Heatmaps --------------------
library(tidyverse)
library(viridis)
library(hrbrthemes)

res_list<- c("P50", "P88", "ISI", "SWP_M", "SWP_L",
             "SWS_M", "SWS_L", "Phe", "Dkleaf_mean","Kleaf_max", 
             "e_W", "e_L","LMA")

# Heatmaps importance [Figures 4b-c] ---------------------------------------

# Plot heatmap importance - scaled rmin --------------------
dt<-read_csv("data/predictor_importance.csv")
glimpse(dt)
dt$Response<- factor(dt$Response, levels = res_list)
heat_vimp_scaled<-dt%>% filter(approach == "scaled")%>%
  ggplot(aes(Response, Predictor, fill= scaled_importance)) + 
  geom_tile()+
  scale_fill_viridis(option = "D", limits = c(0, 1), oob = scales::squish) +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom");heat_vimp_scaled

# save heatmap importance
ggsave(heat_vimp_scaled,file="figures/Figure 4b.png",width=12,height=16, units = "cm", dpi = 300)

# Plot heatmap importance - unscaled rmin ---------------------

heat_vimp_unscaled<-dt%>% filter(approach == "unscaled")%>%
  ggplot(aes(Response, Predictor, fill= scaled_importance)) + 
  geom_tile()+
  scale_fill_viridis(option = "D", limits = c(0, 1), oob = scales::squish) +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom");heat_vimp_unscaled

# save heatmap importance
ggsave(heat_vimp_unscaled,file="figures/Figure 4c.png",width=12,height=16, units = "cm", dpi = 300)


# Heatmaps interactions [Figures 6a-b] --------------------------------------
# Plot heatmap importance - scaled rmin -------------
dt<-read_csv("data/predictor_interaction.csv")
glimpse(dt)
dt$Response<- factor(dt$Response, levels = res_list)
heat_inter_scaled<-dt%>%filter(approach == "scaled")%>%
  ggplot(aes(Response, predictor, fill= interaction)) + 
  geom_tile()+
  scale_fill_viridis(option = "D", limits = c(0,0.5), oob = scales::squish) +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom");heat_inter_scaled

# save heatmap importance
ggsave(heat_inter_scaled,file="figures/Figure 6a.png",width=12,height=16, units = "cm", dpi = 300)

# Plot heatmap importance - unscaled rmin -------------
heat_inter_unscaled<-dt%>%filter(approach == "unscaled")%>%
  ggplot(aes(Response, predictor, fill= interaction)) + 
  geom_tile()+
  scale_fill_viridis(option = "D", limits = c(0,0.5), oob = scales::squish) +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom");heat_inter_unscaled

# save heatmap importance
ggsave(heat_inter_unscaled,file="figures/Figure 6b.png",width=12,height=16, units = "cm", dpi = 300)

# Coplots strongest pairwise interactions [Figures 6c-d]------------------

# ISI - MST_medium x ER_medium --------------------------

data<-read_csv("data/data_for_gbm_trans.csv")
glimpse(data)
setDT(data)

pdat <- expand_grid(MST_medium = seq( min(data$MST_medium,na.rm=T), max(data$MST_medium,na.rm=T), length.out=100),
                    ER_medium = seq( min(data$ER_medium,na.rm=T), max(data$ER_medium,na.rm=T), length.out=100),
                    VD_large = mean(data$VD_large,na.rm=T),
                    ER_small = mean(data$ER_small,na.rm=T),
                    VD_small = mean(data$VD_small,na.rm=T),
                    VD_medium = mean(data$VD_medium,na.rm=T),
                    ER_large =mean(data$ER_large,na.rm=T),
                    MST_small = mean(data$MST_small,na.rm=T),
                    MST_large = mean(data$MST_large,na.rm=T),
                    clade2 = mean(data$clade2,na.rm=T))%>% 
  drop_na()

h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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

pred <- predict(top_mod, newdata=as.h2o(pdat))
pred2<-as.vector(pred$predict)
pdat$pred<-pred2


pair1<-ggplot(data=pdat, aes(MST_medium, ER_medium,fill=pred))+
  geom_tile()+
  scale_fill_viridis_c(option='D', name = "ISI")+
  coord_cartesian(expand=F)+
  labs(title="H-statistic = 0.32",
       x ="MST medium", y = "ER medium")+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14)); pair1
ggsave(pair1, file = "figures/Figure 6c.png",width=14,height=12, units = "cm", dpi = 300)




# SWP_midrib - VD_medium x MST_medium --------------------------

data<-read_csv("data/data_for_gbm_trans_log.csv")
glimpse(data)
setDT(data)

pdat <- expand_grid(MST_medium = seq( min(data$MST_medium,na.rm=T), max(data$MST_medium,na.rm=T), length.out=100),
                    VD_medium = seq( min(data$VD_medium,na.rm=T), max(data$VD_medium,na.rm=T), length.out=100),
                    VD_large = mean(data$VD_large,na.rm=T),
                    ER_small = mean(data$ER_small,na.rm=T),
                    VD_small = mean(data$VD_small,na.rm=T),
                    ER_medium = mean(data$ER_medium,na.rm=T),
                    ER_large =mean(data$ER_large,na.rm=T),
                    MST_small = mean(data$MST_small,na.rm=T),
                    MST_large = mean(data$MST_large,na.rm=T),
                    clade2 = mean(data$clade2,na.rm=T))%>% 
  drop_na()

h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans_log.csv"
venation <- h2o.uploadFile(path = venation_path)

# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_small", "VD_medium", "VD_large",
                 "ER_small", "ER_medium", "ER_large",
                 "MST_small", "MST_medium", "MST_large")
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
                 seed = 16)# set seed for reproducibility

# Get top model 
top_mod <- h2o.get_best_model(v1)
top_mod

pred <- predict(top_mod, newdata=as.h2o(pdat))
pred2<-as.vector(pred$predict)
pdat$pred<-pred2


pair2<-ggplot(data=pdat, aes(VD_medium, MST_medium,fill=pred))+
  geom_tile()+
  scale_fill_viridis_c(option='D', name = "SWP midrib")+
  coord_cartesian(expand=F)+
  labs(title="H-statistic = 0.43",
       x ="VD medium", y = "MST medium")+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14)); pair2
ggsave(pair1, file = "figures/Figure 6d.png",width=14,height=12, units = "cm", dpi = 300)



# Test resilience across venation types -----------
library(tidyverse)
dt<-read_csv("tables/TableS1_spp_list.csv")
glimpse(dt)
kruskal.test(Dkleaf_mean ~ venation_type, data = dt)
kruskal.test(Dkleaf_M ~ venation_type, data = dt)
kruskal.test(Dkleaf_L ~ venation_type, data = dt)
