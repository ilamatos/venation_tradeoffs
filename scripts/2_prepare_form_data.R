###########################################################
# R-code LEAF VENATION ARCHITECTURE-FUNCTION TRADE-OFFS   #
###########################################################

#### PREPARE ARCHITECTURE (FORM) DATA ####################################
# Load packages ----------
library (tidyverse)
library(magrittr)
library(ggpubr)
library(BHPMF)

# Import datasets --------

# Leaf venation architecture dataset
dt_form <- read_csv("data/UCBG_venation_form_orig.csv")
glimpse(dt_form)

# Taxonomic information for 122 species
spp_list<- read_csv("data/UCBG_122spp_list.csv")
glimpse(spp_list)

# Leaf venation trait dataset original (non-imputed) values
trait<-read_csv("data/UCBG_venation_function.csv")
glimpse(trait)

# Leaf venation trait dataset imputed values
imput_final<- read_csv("data/UCBG_venation_function_imput.csv")
glimpse(imput_final)

# Figure out how many scales we have sufficient data at --------
# [Figure S1]
data_at_scale <- dt_form %>% 
  mutate(bin = cut(dt_form$widthThreshold,breaks=100)) %>%
  group_by(bin) %>% 
  summarize(num_codes=length(unique(spp_code))) %>%
  mutate(bin_label=as.numeric(gsub("\\(","",sapply(strsplit(as.character(bin),","),head,1)))) %>%
  filter(bin_label > 0)

g_data_available <- ggplot(data_at_scale,aes(x=bin_label,y=num_codes)) + 
  geom_line(size=1) + 
  theme_bw() +
  geom_vline(xintercept = 0,color='red') +
  geom_vline(xintercept = 0.5,color='red') +
  xlab(expression(paste("r"["min"], " (mm)"))) +
  ylab("Log (Number of samples)")+
  scale_y_log10();g_data_available

# save Figure S1
ggsave(g_data_available, file='figures/Figure S1.png',width=8,height=5,dpi = 300)



# Filter out rmin out of the range -----------

dt_form_filter<- dt_form%>%
  filter(widthThreshold < 0.5)%>% # filtering out rmin > 0.5 mm
  filter (widthThreshold > 0.01) # filtering out rmin < 0.01 mm
glimpse(dt_form_filter)

# Plot leaf venation architecture traits across vein sizes ---------------------

# filter out species without venation architecture data
spp_list_form<- spp_list%>%
  filter (spp_code != "Auc_jap",
          spp_code != "Nym")
glimpse(spp_list_form)

# add clade to the venation architecture dataset
dt_form_join <- dt_form_filter  %>% 
  left_join(spp_list_form, by = "spp_code") # add species taxonomy
glimpse(dt_form_join)

# reorder clade levels 
dt_form_join$clade <- factor(dt_form_join$clade, levels = c("ferns", "basal angiosperms","monocots","basal eudicots", "rosids", "asterids"))
levels(as.factor(dt_form_join$clade))

# plot vein density (VD)
g_vd_clade<- ggplot(dt_form_join, 
                    aes(x=widthThreshold, y=vein_density, 
                        color = clade, group = spp_code )) +
  geom_path(linewidth =2, alpha =.5)+
  theme_bw() +
  scale_x_sqrt()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 40),
        legend.position = "none",
        strip.text.x = element_text(size = 35),
        panel.spacing = unit(2, "lines"))+
  ylab (label = expression("Vein density (mm mm"^-2*")"))+
  xlab(label = expression ("r"["min"]*" (mm)"))+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#0072B2",# dark blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+ # grey
  facet_wrap(~clade); g_vd_clade

# save Figure S4
ggsave(g_vd_clade, file='figures/Figure S4.png', width = 50, height = 40, units = "cm", dpi = 300, limitsize = FALSE)

# Plot Minimum spanning tree ratio (MST)

g_mst_clade<- ggplot(dt_form_join, aes(x=widthThreshold, y=MSTRatio, color = clade, group = spp_code )) +
  geom_path(linewidth =2, alpha =.5)+
  theme_bw() +
  scale_x_sqrt()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 40),
        legend.position = "none",
        strip.text.x = element_text(size = 35),
        panel.spacing = unit(2, "lines"))+
  ylab (label = "Minimum spanning tree ratio")+
  xlab(label = expression ("r"["min"]*" (mm)"))+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#0072B2",# dark blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+ # grey
  facet_wrap(~clade); g_mst_clade
# save Figure S5
ggsave(g_mst_clade, file='figures/Figure S5.png', width = 50, height = 40, units = "cm", dpi = 300, limitsize = FALSE)

# Plot loop elongation ratio (ER)

g_er_clade<- ggplot(dt_form_join, aes(x=widthThreshold, y=Elongation, color = clade, group = spp_code )) +
  geom_line(linewidth =2, alpha =.5)+
  theme_bw() +
  scale_x_sqrt()+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 40),
        legend.position = "none",
        strip.text.x = element_text(size = 35),
        panel.spacing = unit(2, "lines"))+
  ylab (label = "Loop elongation ratio")+
  xlab(label = expression ("r"["min"]*" (mm)"))+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#0072B2",# dark blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+ # grey
  facet_wrap(~clade, scale = "free"); g_er_clade

# save Figure S6
ggsave(g_er_clade, file='figures/Figure S6.png', width = 50, height = 40, units = "cm", dpi = 300, limitsize = FALSE)

# Divide veins into e bins (small, medium, large veins) ------------------
# 1st approach - scaled rmin -----------------------
# Scale rmin values ------------------

names <- unique(dt_form_join$spp_code)
dt_std <- list()
# standardize rmin to vary between 0-1 within each species
for (i in 1:length(names)){
  dt<- dt_form_join%>%
    filter(spp_code == names[[i]])%>%
    mutate(rminstd = widthThreshold/max(widthThreshold))
  dt_std[[i]]<-dt
}
dt_form_std <- bind_rows(dt_std)
glimpse(dt_form_std)

# Divide rmin values into 3 bins (small, medium, large veins) --------------

# divide data into 3 bins
width_bins <- c(0.3,0.6,0.9) 
dt_form_std_bin3 <- dt_form_std %>%
  mutate(rmin_binned = case_when(rminstd <= 0.3 ~ 0.3, # small veins
                                 rminstd > 0.3 & rminstd <= 0.6 ~ 0.6, # medium veins
                                 rminstd > 0.6 ~ 0.9)) # large veins

glimpse(dt_form_std_bin3)

# calculate median values of venation architecture traits for each bin

dt_form_median <- dt_form_std_bin3 %>%
  mutate (ER = Elongation -1)%>% # Elongation ratio - start it at zero for easier fitting later
  group_by(spp_code, rmin_binned) %>%
  summarize(VD.median=median(vein_density, na.rm = T), MST.median=median(MSTRatio, na.rm = T),ER.median=median(ER, na.rm = T),CR.median=median(Circularity, na.rm = T)) %>%
  ungroup

glimpse(dt_form_median)
unique(dt_form_median$spp_code) # checking the total number of species

# insert 0 anywhere the bin exists but the data don't
dt_form_median <- dt_form_median %>% replace(is.na(.), 0)%>%
  left_join(spp_list_form, by="spp_code") # add taxonomy info
unique(dt_form_median$spp_code) # checking the total number of species
glimpse(dt_form_median)

# write venation architecture data 3 bins (scaled rmin)
write_csv(dt_form_median,"data/UCBG_form_binned3_std_median.csv")

# Join leaf venation architecture and function datasets ----------------------------

# re-organize venation data so that each species has a single row
dt_form_final <- dt_form_median %>%
  mutate(rmin_binned = str_c("rmin_binned", rmin_binned))%>%
  pivot_wider(names_from = rmin_binned, values_from = c("VD.median", "MST.median", "ER.median", "CR.median"), names_sep="")%>%
  rename(VD_small = VD.medianrmin_binned0.3,
         VD_medium = VD.medianrmin_binned0.6,
         VD_large = VD.medianrmin_binned0.9,
         MST_small = MST.medianrmin_binned0.3,
         MST_medium = MST.medianrmin_binned0.6,
         MST_large = MST.medianrmin_binned0.9,
         ER_small = ER.medianrmin_binned0.3,
         ER_medium = ER.medianrmin_binned0.6,
         ER_large = ER.medianrmin_binned0.9,
         CR_small = CR.medianrmin_binned0.3,
         CR_medium = CR.medianrmin_binned0.6,
         CR_large = CR.medianrmin_binned0.9)

dt_form_final <- dt_form_final %>% replace(is.na(.), 0)

glimpse(dt_form_final)
unique(dt_form_final$spp_code) # checking number of spp

# 2nd approach - unscaled rmin ---------------------
# Log10-transform rmin ----------------
dt_form_log<-dt_form_join%>% mutate(rminlog = log10(widthThreshold) + 2.5 )
# Divide rmin values into 3 bins (small, medium, large veins) -------------------
width_bins <- c(0.5,1.0,1.5) # dividing data into 3 bins
dt_form_log_bin3 <- dt_form_log %>%
  mutate(rmin_binned = case_when(rminlog <= 1 ~ 0.5, # small veins
                                 rminlog > 1 & rminlog <= 1.5 ~ 1, # medium veins
                                 rminlog > 1.5 ~ 1.5)) # large veins

glimpse(dt_form_log_bin3)

# # write venation architecture data 3 bins (unscaled rmin)
write_csv(dt_form_log_bin3, "UCBG_form_binned3_log_median.csv") 

# calculate median values of venation form traits for each bin

dt_form_median_log <- dt_form_log_bin3 %>%
  mutate (ER = Elongation -1)%>% # Elongation ratio - start it at zero for easier fitting later
  group_by(spp_code, rmin_binned) %>%
  summarize(VD.median=median(vein_density, na.rm = T), MST.median=median(MSTRatio, na.rm = T),ER.median=median(ER, na.rm = T),CR.median=median(Circularity, na.rm = T)) %>%
  ungroup

glimpse(dt_form_median_log)
unique(dt_form_median_log$spp_code) # checking the total number of species

# insert 0 anywhere the bin exists but the data don't
dt_form_median_log <- dt_form_median_log %>% replace(is.na(.), 0)%>%
  left_join(spp_list_form, by="spp_code") # add taxonomy info
unique(dt_form_median_log$spp_code) # checking the total number of species
glimpse(dt_form_median_log)

write_csv(dt_form_median_log,"UCBG_form_binned3_log_median.csv")


# Join leaf venation architecture and function datasets -----------------------
# re-organize venation data so that each species has a single row
dt_form_final_log <- dt_form_median_log %>%
  mutate(rmin_binned = str_c("rmin_binned", rmin_binned))%>%
  pivot_wider(names_from = rmin_binned, values_from = c("VD.median", "MST.median", "ER.median", "CR.median"), names_sep="")%>%
  rename(VD_small = VD.medianrmin_binned0.5,
         VD_medium = VD.medianrmin_binned1,
         VD_large = VD.medianrmin_binned1.5,
         MST_small = MST.medianrmin_binned0.5,
         MST_medium = MST.medianrmin_binned1,
         MST_large = MST.medianrmin_binned1.5,
         ER_small = ER.medianrmin_binned0.5,
         ER_medium = ER.medianrmin_binned1,
         ER_large = ER.medianrmin_binned1.5,
         CR_small = CR.medianrmin_binned0.5,
         CR_medium = CR.medianrmin_binned1,
         CR_large = CR.medianrmin_binned1.5)

dt_form_final_log <- dt_form_final_log %>% replace(is.na(.), 0)

glimpse(dt_form_final_log)
# Prepare data for GBM analysis ---------------------------------------

# non-imputed functional trait data
dt_all <- as.data.frame(left_join(dt_form_final, trait, by = "spp_code"))%>%
  mutate(clade = clade, Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>%
  dplyr::select(spp_code, spp_name_WP, family_WP, order_apg4, clade,
                VD_small, VD_medium, VD_large,
                MST_small, MST_medium, MST_large,
                ER_small, ER_medium, ER_large,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI)
glimpse(dt_all)

data_for_gbm_trans <- dt_all%>%
  # coding clade as numeric to make models interpretation easier
  mutate(clade2 = case_when(clade == "ferns" ~ 1, 
                            clade == "basal angiosperms" ~ 2,
                            clade == "monocots" ~ 3,
                            clade == "basal eudicots" ~ 4,
                            clade == "rosids" ~ 5,
                            clade == "asterids" ~ 6,
                            TRUE ~ 0))%>% #sqrt transform VD and ER traits
mutate(VD_large = sqrt(VD_large), VD_medium= sqrt(VD_medium), VD_small = sqrt(VD_small))%>%
mutate(ER_large = sqrt(ER_large), ER_medium= sqrt(ER_medium), ER_small = sqrt(ER_small))


dt_all_log <- as.data.frame(left_join(dt_form_final_log, trait, by = "spp_code"))%>%
  mutate(Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>%
  dplyr::select(spp_code, spp_name_WP, family_WP, order_apg4, clade,
                VD_small, VD_medium, VD_large,
                MST_small, MST_medium, MST_large,
                ER_small, ER_medium, ER_large,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI)
glimpse(dt_all_log)

data_for_gbm_trans_log <- dt_all_log %>%
  mutate(clade2 = case_when(clade == "ferns" ~ 1, # coding clade as numeric to make models interpretation easier
                            clade == "basal angiosperms" ~ 2,
                            clade == "monocots" ~ 3,
                            clade == "basal eudicots" ~ 4,
                            clade == "rosids" ~ 5,
                            clade == "asterids" ~ 6,
                            TRUE ~ 0))%>%
  #sqrt transform VD and ER traits
  mutate(VD_large = sqrt(VD_large), VD_medium= sqrt(VD_medium), VD_small = sqrt(VD_small))%>%
  mutate(ER_large = sqrt(ER_large), ER_medium= sqrt(ER_medium), ER_small = sqrt(ER_small))


# write data for GBM analysis
# non-imputed trait values
write.csv(data_for_gbm_trans,"data/data_for_gbm_trans.csv")
write_csv(data_for_gbm_trans_log, "data/data_for_gbm_trans_log.csv" )


# imputed functional trait data
dt_all_imputed <- as.data.frame(left_join(dt_form_final, imput_final, by = "spp_code"))%>%
  mutate(clade = clade, Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>%
  dplyr::select(spp_code, spp_name_WP, family_WP, order_apg4, clade,
                VD_small, VD_medium, VD_large,
                MST_small, MST_medium, MST_large,
                ER_small, ER_medium, ER_large,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI)

glimpse(dt_all_imputed)

dt_all_imputed_log <- as.data.frame(left_join(dt_form_final_log, imput_final, by = "spp_code"))%>%
  mutate( Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>%
  dplyr::select(spp_code, spp_name_WP, family_WP, order_apg4, clade,
                VD_small, VD_medium, VD_large,
                MST_small, MST_medium, MST_large,
                ER_small, ER_medium, ER_large,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI)

glimpse(dt_all_imputed_log)

data_for_gbm_imput_trans <- dt_all_imputed%>%
  mutate(clade2 = case_when(clade == "ferns" ~ 1, # coding clade as numeric to make models interpretation easier
                            clade == "basal angiosperms" ~ 2,
                            clade == "monocots" ~ 3,
                            clade == "basal eudicots" ~ 4,
                            clade == "rosids" ~ 5,
                            clade == "asterids" ~ 6,
                            TRUE ~ 0))%>% #sqrt transform VD and ER traits
mutate(VD_large = sqrt(VD_large), VD_medium= sqrt(VD_medium), VD_small = sqrt(VD_small))%>%
  mutate(ER_large = sqrt(ER_large), ER_medium= sqrt(ER_medium), ER_small = sqrt(ER_small))

glimpse(data_for_gbm_imput_trans)

data_for_gbm_imput_trans_log <- data_for_gbm_imput_trans%>%
  mutate(clade2 = case_when(clade == "ferns" ~ 1, # coding clade as numeric to make models interpretation easier
                            clade == "basal angiosperms" ~ 2,
                            clade == "monocots" ~ 3,
                            clade == "basal eudicots" ~ 4,
                            clade == "rosids" ~ 5,
                            clade == "asterids" ~ 6,
                            TRUE ~ 0))%>%
  #sqrt transform VD and ER traits
  mutate(VD_large = sqrt(VD_large), VD_medium= sqrt(VD_medium), VD_small = sqrt(VD_small))%>%
  mutate(ER_large = sqrt(ER_large), ER_medium= sqrt(ER_medium), ER_small = sqrt(ER_small))

glimpse(data_for_gbm_imput_trans_log)

# write data for GBM analysis
# imputed trait values
write_csv(data_for_gbm_imput_trans, "data/data_for_gbm_imput_trans.csv" )
write_csv(data_for_gbm_imput_trans_log, "data/data_for_gbm_imput_trans_log.csv" )

# Histograms rmin [Figure S2]--------------------------------
glimpse(dt_form_join)

hist1<-ggplot(dt_form_join, aes(x=widthThreshold)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 30)+
  geom_density(alpha=.2, fill="grey50")+
  theme_classic()+
  xlab(expression ("Original r"["min"]*" (mm)"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)); hist1

hist2<-ggplot(dt_form_std, aes(x=rminstd)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 30)+
  geom_density(alpha=.2, fill="grey50")+
  theme_classic()+
  xlab(expression ("Scaled r"["min"]*" (mm)"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))+
  geom_vline(xintercept=0.3,  linetype="dashed", size=1)+
  geom_vline(xintercept=0.6,  linetype="dashed", size=1)+
  annotate("text", x = 0.15, y = 6.5, label = "small veins", size = 5)+
  annotate("text", x = 0.45, y = 6.5, label = "medium veins", size = 5)+
  annotate("text", x = 0.75, y = 6.5, label = "large veins", size = 5); hist2


hist3<-ggplot(dt_form_join, aes(x=log10(widthThreshold)+2.5)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 30)+
  geom_density(alpha=.2, fill="grey50")+
  theme_classic()+
  xlab(expression ("Log10 (+2.5) r"["min"]*" (mm)"))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))+
  geom_vline(xintercept=1,  linetype="dashed", size=1)+
  geom_vline(xintercept=1.5,  linetype="dashed", size=1)+
  annotate("text", x = 0.75, y = 2.5, label = "small veins", size = 5)+
  annotate("text", x = 1.25, y = 2.5, label = "medium veins", size = 5)+
  annotate("text", x = 1.75, y = 2.5, label = "large veins", size = 5); hist3

g_hist<-ggarrange(hist1, hist2, hist3, nrow =3, ncol =1,
                  labels=c("(a)", "(b)", "(c)"));g_hist

ggsave(g_hist,file='Figure S2.png',width=18,height=32, units = "cm", dpi = 300)

# Prepare data for PCA ---------------------------------

# filter out rmin < 0.01
dt_form2 <- dt_form%>% filter(widthThreshold >0.01)
glimpse(dt_form2)

# define maximum and minimum rmin values
width_low = 0.01 # minimum rmin 
width_high = 0.5 # maximum rmin

# bin the axis-wise data into 50 bins
width_bins <- seq(width_low, width_high,length.out=50)
dt_form2$rmin_binned <- cut(dt_form2$widthThreshold,breaks=width_bins,labels=FALSE,include.lowest=TRUE)
dt_form2$rmin_binned <- width_bins[dt_form2$rmin_binned]

glimpse(dt_form2)
all_scales = expand.grid(rmin_binned=width_bins)

# get bin-median values 
dt_median <- dt_form2 %>% 
  group_by(spp_code, rmin_binned) %>%
  summarize(VD.median=median(vein_density), MST.median=median(MSTRatio),ER.median=median(Elongation)-1,CR.median=median(Circularity)) %>%
  ungroup
glimpse(dt_median)
unique(dt_median$spp_code) # checking the total number of species

# add species taxonomic info
dt_binned<- dt_median%>%
  left_join(spp_list_form, by="spp_code")%>%
  rename(species=spp_name_WP, family = family_WP, order = order_apg4)
glimpse(dt_binned)
unique(dt_binned$spp_code) # cheking number of species


# join form (50 bins) and function (imputed values) datasets for PCA
data_for_pca <-as.data.frame(left_join(dt_binned, imput_final, by = "spp_code"))%>%
  mutate(Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>% # calculate Dkleaf_mean
  dplyr::select(spp_code, species, family,order, clade,
                VD.median,
                MST.median,
                ER.median,
                CR.median,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI, rmin_binned)%>%
  mutate(VD.median=sqrt(VD.median), ER.median=sqrt(ER.median))%>% # sqrt-transform data for more normality
  na.omit

glimpse(data_for_pca)
sum(is.na(data_for_pca)) # no NAs in the PCA data

# save data for PCA
write_csv(data_for_pca, "data/data_for_pca.csv")
# END