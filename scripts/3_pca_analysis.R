###########################################################
# R-code LEAF VENATION ARCHITECTURE-FUNCTION TRADE-OFFS   #
###########################################################

#### PRINCIPAL COMPONENT ANALYSIS (PCA) ##################

# Load packages --------------
library(tidyverse)
library(vegan)
library(ggrepel)
library(viridis)
library(ggpubr)

# Import datasets ---------
data_for_pca<- read_csv("data/data_for_pca.csv")
glimpse(data_for_pca)

# Run PCA -----------------
pca_all <- prcomp(data_for_pca%>%select_if(.,is.numeric)%>%
                    select(-rmin_binned, -Dkleaf_M, -Dkleaf_L, -CR.median), # removing variables not for pca
                  center=TRUE, # center data
                  scale=TRUE) # z-transform data
pca_all

# save PCA loadings [Table S3]
write_csv(as.data.frame(pca_all$rotation), "tables/TableS3_PCA_loadings.csv")

# Broken stick method 
# how many PCs to retain [Figure S7]
png(filename="figures/Figure S7.png", width=1200, height=650)
screeplot(pca_all, main = "Screeplot of Wolf Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")
dev.off()

# extract scores
pca_trajectories <- data.frame(clade=as.factor(data_for_pca$clade),
                               code=as.factor(data_for_pca$spp_code), 
                               rmin=data_for_pca$rmin_binned, 
                               pca_all$x[,c(1,2,3)])

# save PCA scores [Table S4]
write_csv(as.data.frame(pca_trajectories), "tables/TableS4_PCA_scores.csv")

# variance explained
varexp = 100*pca_all$sdev^2/sum(pca_all$sdev^2)
round(varexp, 2)

# get axes
pca_rotations = as.data.frame(pca_all$rotation[,1:2])
pca_rotations$var = row.names(pca_rotations)
pca_rotations$x0 = 0
pca_rotations$y0 = 0
row.names(pca_rotations) = NULL

scale_factor = 10

# Plot PCA results -----------------------------

# overall PCA space - by clade PC1 x PC2 [Figure 3]
g_pca_clade1 <- ggplot(pca_trajectories %>% mutate(clade.factor=clade),
                       aes(x=PC1, y=PC2,col=clade.factor,group=clade.factor)) + 
  geom_point(size = 1.5, alpha = 0.1) + 
  stat_ellipse(size =0.5, alpha = 0.8) + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC1,yend=scale_factor*PC2),inherit.aes = FALSE,size=0.5,color='grey20') +
  geom_text_repel(data=pca_rotations,aes(x=scale_factor*PC1,y=scale_factor*PC2, label = var),inherit.aes = FALSE,size=4,color='black',
                  box.padding = 0.5, max.overlaps = Inf, nudge_x = -0.1, direction = "x") +  
  theme_bw() +
  xlab(sprintf("PC1 (%.2f%%)",varexp[1])) +
  ylab(sprintf("PC2 (%.2f%%)",varexp[2])) +
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" = "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999")) +
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#0072B2",#dark blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"), name = "Clades",
                     breaks=c('ferns', 'basal angiosperms', 'monocots', 'basal eudicots', 'rosids', 'asterids'))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position="right");g_pca_clade1


# overall PCA space - by rmin PC1 x PC2
g_pca_rmin1 <- ggplot(pca_trajectories %>% mutate(rmin.factor=rmin),
                      aes(x=PC1, y=PC2,col=rmin.factor,group=rmin.factor)) + 
  geom_point(size = 1.5, alpha = 0.1) + 
  stat_ellipse(size =0.5, alpha = 0.8) + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC1,yend=scale_factor*PC2),inherit.aes = FALSE,size=0.5,color='grey20') +
  geom_text_repel(data=pca_rotations,aes(x=scale_factor*PC1,y=scale_factor*PC2,label=gsub("\\.median","",var)),inherit.aes = FALSE,size=4,color='black',
                  box.padding = 0.5, max.overlaps = Inf, nudge_x = -0.1, direction = "x") +    
  theme_bw() +
  xlab(sprintf("PC1 (%.2f%%)",varexp[1])) +
  ylab(sprintf("PC2 (%.2f%%)",varexp[2])) +
  scale_fill_viridis(name=expression(paste("r"["min"], " (mm)"))) +
  scale_color_viridis(name=expression(paste("r"["min"], " (mm)")))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        legend.key.size = unit(1, 'cm'),
        legend.position="bottom");g_pca_rmin1 

g_pca1<-ggarrange(g_pca_clade1, g_pca_rmin1,
                  common.legend=F,legend='right',
                  align='hv',
                  nrow=2,ncol=1,
                  labels=c("(a)", "(b)"));g_pca1
# save Figure 3
ggsave(g_pca1,file='figures/Figure 3.png',width=18,height=32, units = "cm", dpi = 300)


# overall PCA space - by clade PC1 x PC3 [Figure S8]----------------

# get axes
pca_rotations = as.data.frame(pca_all$rotation[,1:3])
pca_rotations$var = row.names(pca_rotations)
pca_rotations$x0 = 0
pca_rotations$y0 = 0
row.names(pca_rotations) = NULL

scale_factor = 10

g_pca_clade2 <- ggplot(pca_trajectories %>% mutate(clade.factor=clade),
                       aes(x=PC1, y=PC3,col=clade.factor,group=clade.factor)) + 
  geom_point(size = 1.5, alpha = 0.1) + 
  stat_ellipse(size =0.5, alpha = 0.8) + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC1,yend=scale_factor*PC3),inherit.aes = FALSE,size=0.5,color='grey20') +
  geom_text(data=pca_rotations,aes(x=scale_factor*PC1,y=scale_factor*PC3, label=gsub("\\.median","",var)),inherit.aes = FALSE,size=4,color='black') +  
  theme_bw() +
  xlab(sprintf("PC1 (%.2f%%)",varexp[1])) +
  ylab(sprintf("PC3 (%.2f%%)",varexp[3])) +
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" = "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999")) +
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#0072B2",#dark blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"), name = "Clades",
                     breaks=c('ferns', 'basal angiosperms', 'monocots', 'basal eudicots', 'rosids', 'asterids'))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position="right");g_pca_clade2


# overall PCA space - by rmin PC1 x PC2
g_pca_rmin2 <- ggplot(pca_trajectories %>% mutate(rmin.factor=rmin),
                      aes(x=PC1, y=PC3,col=rmin.factor,group=rmin.factor)) + 
  geom_point(size = 1.5, alpha = 0.1) + 
  stat_ellipse(size =0.5, alpha = 0.8) + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC1,yend=scale_factor*PC3),inherit.aes = FALSE,size=0.5,color='grey20') +
  geom_text(data=pca_rotations,aes(x=scale_factor*PC1,y=scale_factor*PC3,label=gsub("\\.median","",var)),inherit.aes = FALSE,size=4,color='black') +    
  theme_bw() +
  xlab(sprintf("PC1 (%.2f%%)",varexp[1])) +
  ylab(sprintf("PC3 (%.2f%%)",varexp[3])) +
  scale_fill_viridis(name=expression(paste("r"["min"], " (mm)"))) +
  scale_color_viridis(name=expression(paste("r"["min"], " (mm)")))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        # legend.key.size = unit(1, 'cm'),
        legend.position="bottom");g_pca_rmin2 


# overall PCA space - by clade PC2 x PC3 ----------------

# get axes
pca_rotations = as.data.frame(pca_all$rotation[,2:3])
pca_rotations$var = row.names(pca_rotations)
pca_rotations$x0 = 0
pca_rotations$y0 = 0
row.names(pca_rotations) = NULL

scale_factor = 10

g_pca_clade3 <- ggplot(pca_trajectories %>% mutate(clade.factor=clade),
                       aes(x=PC2, y=PC3,col=clade.factor,group=clade.factor)) + 
  geom_point(size = 1.5, alpha = 0.1) + 
  stat_ellipse(size =0.5, alpha = 0.8) + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC2,yend=scale_factor*PC3),inherit.aes = FALSE,size=0.5,color='grey20') +
  geom_text(data=pca_rotations,aes(x=scale_factor*PC2,y=scale_factor*PC3, label=gsub("\\.median","",var)),inherit.aes = FALSE,size=4,color='black') +  
  theme_bw() +
  xlab(sprintf("PC2 (%.2f%%)",varexp[2])) +
  ylab(sprintf("PC3 (%.2f%%)",varexp[3])) +
  scale_fill_manual(values = c("asterids"="#D55E00",# dark orange
                               "rosids"="#CC79A7", # pink
                               "basal angiosperms" = "#0072B2", #dark blue
                               "monocots" = "#009E73",# green
                               "basal eudicots" = "#F0E442",#yellow
                               "ferns" = "#999999")) +
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#0072B2",#dark blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"), name = "Clades",
                     breaks=c('ferns', 'basal angiosperms', 'monocots', 'basal eudicots', 'rosids', 'asterids'))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position="right");g_pca_clade3


# overall PCA space - by rmin PC2 x PC3
g_pca_rmin3 <- ggplot(pca_trajectories %>% mutate(rmin.factor=rmin),
                      aes(x=PC2, y=PC3,col=rmin.factor,group=rmin.factor)) + 
  geom_point(size = 1.5, alpha = 0.1) + 
  stat_ellipse(size =0.5, alpha = 0.8) + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC2,yend=scale_factor*PC3),inherit.aes = FALSE,size=0.5,color='grey20') +
  geom_text(data=pca_rotations,aes(x=scale_factor*PC2,y=scale_factor*PC3,label=gsub("\\.median","",var)),inherit.aes = FALSE,size=4,color='black') +    
  theme_bw() +
  xlab(sprintf("PC2 (%.2f%%)",varexp[1])) +
  ylab(sprintf("PC3 (%.2f%%)",varexp[2])) +
  scale_fill_viridis(name=expression(paste("r"["min"], " (mm)"))) +
  scale_color_viridis(name=expression(paste("r"["min"], " (mm)")))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        # legend.key.size = unit(1, 'cm'),
        legend.position="bottom");g_pca_rmin3

# arrange Figure S8
g_pca_all<-ggarrange(g_pca_clade2, g_pca_rmin2,
                     g_pca_clade3, g_pca_rmin3,
                     common.legend=F,legend='bottom',
                     align='hv',
                     nrow=2,ncol=2,
                     labels=c("(a)", "(b)", "(c)", "(d)"));g_pca_all
# save Figure S8
ggsave(g_pca_all,file='figures/Figure S8.png',width=32,height=32, units = "cm", dpi = 300)
# END