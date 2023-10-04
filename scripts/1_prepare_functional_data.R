###########################################################
# R-code LEAF VENATION ARCHITECTURE-FUNCTION TRADE-OFFS   #
###########################################################

#### PREPARE FUNCTIONAL DATA ####################################
# Load packages ----------
library (tidyverse)
library(magrittr)
library(ggpubr)
library(BHPMF)

# Import datasets --------

# Leaf venation function dataset
fun_data<- read_csv("data/UCBG_venation_function.csv") 
glimpse(fun_data)

# Taxonomic information for 122 species
spp_list<- read_csv("data/UCBG_122spp_list.csv")
glimpse(spp_list)

# Impute missing data using BHPMF ------- 

# BHPMF - Schrodt et al 2015 Glob. Ecol. Biogeogr

# set working directory
setwd("data/bays_imput") 

# specify a temporary output directory
tmp.dir <- dirname("data/bays_imput/tmp")

# create species list with phylogenetic info (species name, genus and family)
phylo <- spp_list %>% mutate(species_id = seq(from= 1, to = 122, by =1),
                             species = spp_name_WP,
                             genus = genus_WP,
                             family = family_WP)%>%
  dplyr::select(species_id, species, genus, family)

glimpse(phylo)

# organize functional trait dataset
trait <- fun_data %>% dplyr::select (spp_code, LMA, SWP_M, SWP_L,
                                     SWS_M, SWS_L, e_W, e_L, Dkleaf_M,
                                     Dkleaf_L, Kleaf_max, P50, P88, Phe, ISI)
glimpse(trait)

# check if number of rows is the same for both datasets
nrow(phylo) == nrow(trait)

# function to z-transform data
scale2 <- function(x, na.rm = T) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# log10 and z-transform trait data
trait_trans<- trait %>%dplyr::select(-spp_code)%>% # exclude spp_code column
  mutate (Dkleaf_M = Dkleaf_M + 100,
          Dkleaf_L = Dkleaf_L + 100)%>% # adding 100 to DKleaf (because we cannot take log10 from neg values)
  mutate_all (.,log10)%>% # take the log10
  mutate_all(.,scale2) # z-transform variables

glimpse(trait_trans)

# gap fill missing data using BHPMF - repeat imputation 50 times
for (i in 1:50){
  capture.output(GapFilling(as.matrix(trait_trans), as.data.frame(phylo),
                            mean.gap.filled.output.path = paste0("data/bays_imput/means/imput_means", i,".csv"),
                            std.gap.filled.output.path=paste0("data/bays_imput/sd/imput_sd", i,".csv"),
                            rmse.plot.test.data = T), file=paste0("data/bays_imput/RMSE/RMSE", i,".txt"))
}

# imput_means.csv = contains the imputed mean values
# imput_sd.csv = contains the standard deviation for the imputed values
# RMSE = contains the average root mean square error for each imputation

# Back transform imputed values ---------------------

# take the log10 mean of the original (non-imputed) trait data
tmean<- trait%>% dplyr::select(-spp_code)%>% #remove spp_code column
  mutate (Dkleaf_M = Dkleaf_M + 100,
          Dkleaf_L = Dkleaf_L + 100)%>% # adding 100 to DKleaf (because we cannot take log10 from neg values)mutate_all(.,log10)%>%
  mutate_all(.,log10)%>%
  summarise_all(., mean, na.rm = T) # take the log10 mean
glimpse(tmean)

# take the log10 sd of the original (non-imputed) trait data
tsd<- trait%>% dplyr::select(-spp_code)%>% #remove spp_code column
  mutate (Dkleaf_M = Dkleaf_M + 100,
          Dkleaf_L = Dkleaf_L + 100)%>% # adding 100 to DKleaf (because we cannot take log10 from neg values)mutate_all(.,log10)%>%
  mutate_all(.,log10)%>%
  summarise_all(., sd, na.rm = T) # take the log10 mean
glimpse(tsd)

# read each of the 50 files with imputed mean data
files = list.files("data/bays_imput/means", pattern = "*.csv", full.names = TRUE) 
length(files) # length should be 50

# back transform the imputed mean values for each trait
back_list<- list()
for (i in 1:50){
  back_list[[i]]<- read_delim(files[i], delim = "\t")%>%
    mutate(N = i,
           LMA = 10^((LMA*tsd$LMA)+tmean$LMA), 
           SWP_M = 10^((SWP_M*tsd$SWP_M)+tmean$SWP_M),
           SWP_L = 10^((SWP_L*tsd$SWP_L)+tmean$SWP_L),
           SWS_M = 10^((SWS_M*tsd$SWS_M)+tmean$SWS_M),
           SWS_L = 10^((SWS_L*tsd$SWS_L)+tmean$SWS_L),
           e_W = 10^((e_W*tsd$e_W)+tmean$e_W),
           e_L = 10^((e_L*tsd$e_L)+tmean$e_L),
           Dkleaf_M = 10^((Dkleaf_M*tsd$Dkleaf_M)+tmean$Dkleaf_M),
           Dkleaf_L = 10^((Dkleaf_L*tsd$Dkleaf_L)+tmean$Dkleaf_L),
           Kleaf_max = 10^((Kleaf_max*tsd$Kleaf_max)+tmean$Kleaf_max),
           P50 = 10^((P50*tsd$P50)+tmean$P50),
           P88 = 10^((P88*tsd$P88)+tmean$P88),
           Phe = 10^((Phe*tsd$Phe)+tmean$Phe),
           ISI = 10^((ISI*tsd$ISI)+tmean$ISI))
}

# list of imputed and back-transformed values
glimpse(back_list) 

# join imputed and back transformed results 
imput_means<-back_list%>% 
  bind_rows() 
glimpse(imput_means) 

# calculate the range of trait values in the original dataset
tmax<- trait%>% dplyr::select(-spp_code)%>% #remove spp_code column
  summarise_all(., max, na.rm = T) # maximum value for each trait
glimpse(tmax)

tmin<- trait%>% dplyr::select(-spp_code)%>% #remove spp_code column
  summarise_all(., min, na.rm = T) # minimum value for each trait
glimpse(tmin)

# filter out imputed values "out of range"
t_mean <- list()
t_sd <- list()
for (i in 1:length(tmax)){
  t<-sapply(back_list,"[[",i) # select trait across lists 
  t_filter<-t%>% # filtering out trait values out of the range
    replace(.,. >as.double(tmax[,i]), NA)%>% # filter out values over the range
    replace(.,. <as.double(tmin[,i]), NA) # filter out values below the range
  t_mean[[i]]<-as.data.frame(apply(t_filter, 1, mean, na.rm = T)) # calculate trait mean of imputed values
  t_sd[[i]]<-as.data.frame(apply(t_filter, 1, sd, na.rm = T)) # calculate trait sd of imputed values
}

# create datadet woth imputed, back-transformed, and filtered mean trait values
imput_means<-t_mean%>% bind_cols()%>%
  magrittr::set_colnames(c("LMA_imput", "SWP_M_imput", "SWP_L_imput", "SWS_M_imput", "SWS_L_imput",
                           "e_W_imput", "e_L_imput", "Dkleaf_M_imput", "Dkleaf_L_imput", "Kleaf_max_imput", "P50_imput", 
                           "P88_imput", "Phe_imput", "ISI_imput"))%>%
  mutate(Dkleaf_M_imput = Dkleaf_M_imput - 100,
         Dkleaf_L_imput = Dkleaf_L_imput - 100)
glimpse(imput_means)

# combine original and imputed mean trait values into a single data frame
comp_mean<-cbind(trait, imput_means)
glimpse(comp_mean)

# Validate imputation -------------------------------------
# plot original vs. imputed trait mean values to validate imputation
# [Appendix S1 - Figure 2]

g_kmax<-ggplot(comp_mean, aes(x=Kleaf_max, y = Kleaf_max_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,40)+
  ylim(0,40)+
  ylab(label = expression(atop("Kleaf"["max"]*" (mmol m"^-2*"s"^-1*"MPa"^-1*")",
                               "imputed")))+
  xlab(label = expression(atop("Kleaf"["max"]*" (mmol m"^-2*"s"^-1*"MPa"^-1*")",
                               "original")))+
  stat_regline_equation(label.y = 40, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 38, aes(label = ..rr.label..), size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_kmax

g_p50<-ggplot(comp_mean, aes(x=P50, y = P50_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,3.5)+
  ylim(0,3.5)+
  ylab(label = "P50 (-MPa) imputed")+
  xlab(label = "P50 (-MPa) original")+
  stat_regline_equation(label.y = 3.5, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 3.2, aes(label = ..rr.label..), size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_p50

g_p88<-ggplot(comp_mean, aes(x=P88, y = P88_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,8)+
  ylim(0,8)+
  ylab(label = "P88 (-MPa) imputed")+
  xlab(label = "P88 (-MPa) original")+
  stat_regline_equation(label.y = 8, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 7.5, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_p88

g_isi<-ggplot(comp_mean, aes(x=P88, y = P88_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,10)+
  ylim(0,10)+
  ylab(label = "ISI imputed")+
  xlab(label = "ISI original")+
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 9.2, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_isi

g_swpm<-ggplot(comp_mean, aes(x=SWP_M, y = SWP_M_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,1500)+
  ylim(0,1500)+
  ylab(label = expression(atop("SWP"["midrib"]*" (kJ m"^-2*"m"^-1*")",
                               "imputed")))+
  xlab(label = expression(atop("SWP"["midrib"]*" (kJ m"^-2*"m"^-1*")",
                               "original")))+
  stat_regline_equation(label.y = 1500, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 1400, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_swpm

g_swpl<-ggplot(comp_mean, aes(x=SWP_L, y = SWP_L_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,800)+
  ylim(0,800)+
  ylab(label = expression(atop("SWP"["lamina"]*" (kJ m"^-2*"m"^-1*")",
                               "imputed")))+
  xlab(label = expression(atop("SWP"["lamina"]*" (kJ m"^-2*"m"^-1*")",
                               "original")))+
  stat_regline_equation(label.y = 800, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 750, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_swpl

g_swsm<-ggplot(comp_mean, aes(x=SWS_M, y = SWS_M_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,10000)+
  ylim(0,10000)+
  ylab(label = expression(atop("SWS"["midrib"]*" (J m"^-2*")",
                               "imputed")))+
  xlab(label = expression(atop("SWS"["midrib"]*" (J m"^-2*")",
                               "original")))+
  stat_regline_equation(label.y = 10000, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 9500, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_swsm

g_swsl<-ggplot(comp_mean, aes(x=SWS_L, y = SWS_L_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,50000)+
  ylim(0,50000)+
  ylab(label = expression(atop("SWS"["lamina"]*" (J m"^-2*")",
                               "imputed")))+
  xlab(label = expression(atop("SWS"["lamina"]*" (J m"^-2*")",
                               "original")))+
  stat_regline_equation(label.y = 50000, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 45000, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_swsl

comp_mean2<- comp_mean%>%
  mutate(Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2,
         Dkleaf_mean_imput = (Dkleaf_M_imput + Dkleaf_L_imput)/2)

g_dkleafm<-ggplot(comp_mean2, aes(x=Dkleaf_mean, y = Dkleaf_mean_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(-100,400)+
  ylim(-100,400)+
  ylab(label = expression(atop("Dkleaf"["mean"]*"(%)",
                               "imputed")))+
  xlab(label = expression(atop("Dkleaf"["mean"]*"(%)",
                               "original")))+
  stat_regline_equation(label.y = 400, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 350, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_dkleafm

g_ew<-ggplot(comp_mean, aes(x=e_W, y = e_W_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,600)+
  ylim(0,600)+
  ylab(label = expression(atop("e"["whole"]*" (MN m"^-2*")",
                               "imputed")))+
  xlab(label = expression(atop("e"["whole"]*" (MN m"^-2*")",
                               "original")))+
  stat_regline_equation(label.y = 600, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 550, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_ew

g_el<-ggplot(comp_mean, aes(x=e_L, y = e_L_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,600)+
  ylim(0,600)+
  ylab(label = expression(atop("e"["lamina"]*" (MN m"^-2*")",
                               "imputed")))+
  xlab(label = expression(atop("e"["lamina"]*" (MN m"^-2*")",
                               "original")))+
  stat_regline_equation(label.y = 600, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 550, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_el

g_lma<-ggplot(comp_mean, aes(x=LMA, y = LMA_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,300)+
  ylim(0,300)+
  ylab(label = expression(atop("LMA (m"^-2*")",
                               "imputed")))+
  xlab(label = expression(atop("LMA (m"^-2*")",
                               "original")))+
  stat_regline_equation(label.y = 300, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 250, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_lma

g_phe<-ggplot(comp_mean, aes(x=Phe, y = Phe_imput))+
  geom_point(color = "gray20", alpha =.3, size =5)+
  geom_smooth(method = "lm", color ="red")+
  xlim(0,10)+
  ylim(0,10)+
  ylab(label = expression(atop("Phenol content (g g"^-1*")",
                               "imputed")))+
  xlab(label = expression(atop("Phenol content (g g"^-1*")",
                               "original")))+
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 3) +
  stat_regline_equation(label.y = 9, aes(label = ..rr.label..),  size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10));g_phe

# combine all figures
g_imput<-ggarrange(g_kmax, g_p50, g_p88,g_isi,
                   g_swpm, g_swpl, g_swsm, g_swsl,
                   g_dkleafm, g_ew, g_el, g_lma, g_phe,
                   common.legend=TRUE,legend='bottom',
                   align='hv',
                   nrow=5,ncol=3,
                   labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                            "(f)", "(g)", "(h)", "(i)", "(j)", 
                            "(k)", "(l)", "(m)"));g_imput

# save Appendix S1 - Figure 2
ggsave(g_imput,file='figures/AppendixS1_Fig2.png',width=24,height=36, units = "cm")

# replace NAs in the original trait data with the imputed values
imput_final<-comp_mean%>%
  mutate(SWP_M = coalesce(SWP_M,SWP_M_imput),
         SWS_M = coalesce(SWS_M,SWS_M_imput),
         SWS_L = coalesce(SWS_L,SWS_L_imput),
         e_W = coalesce(e_W,e_W_imput),
         e_L = coalesce(e_L,e_L_imput),
         Dkleaf_M = coalesce(Dkleaf_M,Dkleaf_M_imput),
         Dkleaf_L = coalesce(Dkleaf_L,Dkleaf_L_imput),
         P50 = coalesce(P50,P50_imput),
         P88 = coalesce(P88,P88_imput),
         Phe = coalesce(Phe,Phe_imput))%>%
  dplyr::select(spp_code,LMA, SWP_M, SWP_L, SWS_M, SWS_L,
                e_W, e_L, Dkleaf_M, Dkleaf_L, Kleaf_max,
                P50, P88, Phe, ISI)

glimpse(imput_final)
sum(is.na(imput_final)) # no NAs in the final imputed trait dataset

# save imputed functional trait dataset
write_csv(imput_final, "data/UCBG_venation_function_imput.csv")


# Plot functional traits across clades ---------------------------
# [Figure S3]

# add clade column to functional datasets 
trait$clade <- spp_list$clade # original trait values
imput_final$clade <- spp_list$clade # imputed trait values

# reorder clade levels 
trait$clade <- factor(trait$clade, levels = c("ferns", "basal angiosperms","monocots","basal eudicots", "rosids", "asterids"))
levels(as.factor(trait$clade))

imput_final$clade <- factor(imput_final$clade, levels = c("ferns", "basal angiosperms","monocots","basal eudicots", "rosids", "asterids"))
levels(as.factor(imput_final$clade))

# Violin plots + kruskal wallis test + pairwise wilcox test

# Kleaf_max 
kruskal.test(Kleaf_max ~ clade, data = trait)

pairwise.wilcox.test(trait$Kleaf_max, trait$clade,
                     p.adjust.method = "BH")

g_kmax_box<-ggplot(trait, aes(x=clade, y=Kleaf_max, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = Kleaf_max), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab((label = expression(atop("Kleaf"["max"], 
                                "(mmol m"^-2*"s"^-1*"MPa"^-1*")"))))+
  annotate("text", label = c("b", "ab", "ab", "ab", "ab", "a"), x = c(1,2,3,4,5,6), y =35, size = 5);g_kmax_box

# P50
kruskal.test(P50 ~ clade, data = trait)

g_p50_box<-ggplot(trait, aes(x=clade, y=P50, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = P50), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = "P50 (-MPa)");g_p50_box

# P88
kruskal.test(P88 ~ clade, data = trait)

g_p88_box<-ggplot(trait, aes(x=clade, y=P88, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = P88), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = "P88 (-MPa)");g_p88_box

# ISI
kruskal.test(ISI ~ clade, data = trait)

pairwise.wilcox.test(trait$ISI, trait$clade,
                     p.adjust.method = "BH")

g_isi_box<-ggplot(trait, aes(x=clade, y=ISI, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = ISI), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = "Implosion safety index")+
  annotate("text", label = c("bc", "b", "c", "ab", "ab", "a"), x = c(1,2,3,4,5,6), y =0.85, size = 5);g_isi_box

# SWP_M
kruskal.test(SWP_M ~ clade, data = trait)

g_swpm_box<-ggplot(trait, aes(x=clade, y=SWP_M, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = SWP_M), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("SWP"["midrib"]*" (kJ m"^-2*"m"^-1*")"));g_swpm_box

# SWP_L
kruskal.test(SWP_L ~ clade, data = trait)

pairwise.wilcox.test(trait$SWP_L, trait$clade,
                     p.adjust.method = "BH")

g_swpl_box<-ggplot(trait, aes(x=clade, y=SWP_L, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = SWP_L), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("SWP"["lamina"]*" (kJ m"^-2*"m"^-1*")"))+
  annotate("text", label = c("ab", "ab", "a", "a", "b", "ab"), x = c(1,2,3,4,5,6), y =800, size = 5);g_swpl_box

# SWS_M
kruskal.test(SWS_M ~ clade, data = trait)

g_swsm_box<-ggplot(trait, aes(x=clade, y=SWS_M, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = SWS_M), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("SWS"["midrib"]*" (J m"^-2*")"));g_swsm_box

# SWS_L
kruskal.test(SWS_L ~ clade, data = trait)

pairwise.wilcox.test(trait$SWS_L, trait$clade,
                     p.adjust.method = "BH")

g_swsl_box<-ggplot(trait, aes(x=clade, y=SWS_L, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = SWS_L), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("SWS"["lamina"]*" (J m"^-2*")"))+
  annotate("text", label = c("a", "b", "b", "b", "b", "b"), x = c(1,2,3,4,5,6), y =50000, size = 5);g_swsl_box

# DKleaf_mean
kruskal.test(Dkleaf_mean ~ clade, data = trait%>% mutate (Dkleaf_mean = (Dkleaf_M+ Dkleaf_L)/2))

g_dkleaf_box<-trait%>% mutate (Dkleaf_mean = (Dkleaf_M+ Dkleaf_L)/2)%>% # calculate Dkleaf_mean
  ggplot(., aes(x=clade, y=Dkleaf_mean, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final%>% mutate (Dkleaf_mean = (Dkleaf_M+ Dkleaf_L)/2),  # calculate Dkleaf_mean
              aes (x = clade, y =Dkleaf_mean), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("Dkleaf"["mean"]*" (%)"));g_dkleaf_box

# e_W
kruskal.test(e_W ~clade, data = trait)

pairwise.wilcox.test(trait$e_W, trait$clade,
                     p.adjust.method = "BH")

g_ew_box<-ggplot(trait, aes(x=clade, y=e_W, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = e_W), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("e"["whole"]*" (MN m"^-2*")"))+
  annotate("text", label = c("ab", "bc", "a", "bc", "c", "bc"), x = c(1,2,3,4,5,6), y =520, size = 5);g_ew_box

# LMA
kruskal.test(LMA ~ clade, data = trait)

g_lma_box<-ggplot(trait, aes(x=clade, y=LMA, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = LMA), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("LMA (m"^-2*")"));g_lma_box

# Phe
kruskal.test(Phe ~ clade, data = trait)

pairwise.wilcox.test(trait$Phe, trait$clade,
                     p.adjust.method = "BH")
g_phe_box<-ggplot(trait, aes(x=clade, y=Phe, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = Phe), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("Phenol content (g g"^-1*")"))+
  annotate("text", label = c("a", "a", "b", "a", "a", "a"), x = c(1,2,3,4,5,6), y =8, size = 5);g_phe_box

# e_L
kruskal.test(e_L ~ clade, data = trait)

g_el_box<-ggplot(trait, aes(x=clade, y=e_L, fill=clade)) + 
  geom_violin()+
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" =  "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999"))+
  geom_jitter(data = imput_final, aes (x = clade, y = e_L), inherit.aes = F, shape = 1, size=1, alpha=0.3, width = 0.2, stroke = 1.2)+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none",
        axis.text.x=element_blank())+
  ylab(label = expression("e"["lamina"]*" (MN m"^-2*")"));g_el_box
#annotate("text", label = c("a", "a", "b", "a", "a", "a"), x = c(1,2,3,4,5,6), y =8, size = 5)

# make Figure S3

g_fun_violin<-ggarrange(g_p50_box, g_p88_box,g_isi_box,
                        g_swpm_box, g_swpl_box, g_swsm_box, g_swsl_box, g_phe_box,
                        g_dkleaf_box,g_kmax_box, g_ew_box, g_el_box,g_lma_box, 
                        common.legend=TRUE,legend='bottom',
                        align='hv',
                        nrow=5,ncol=3,
                        labels=c("(a)", "(b)", "(c)", "(d)", "(e)",
                                 "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)"));g_fun_violin

ggsave(g_fun_violin,file='figures/Figure S3.png',width=24,height=40, units = "cm")
# END
