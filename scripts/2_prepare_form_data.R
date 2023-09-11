###################################################
# R-code LEAF VENATION FORM-FUNCTION TRADE-OFFS   #
###################################################

#### PREPARE FORM DATA ####################################
# Load packages ----------
library (tidyverse)
library(magrittr)
library(ggpubr)
library(BHPMF)

# Import datasets --------

# Taxonomic information for 122 species
spp_list<- read_csv("data/UCBG_122spp_list.csv")
glimpse(spp_list)

# Leaf venation trait dataset original (non-imputed) values
trait<-read_csv("data/UCBG_venation_function.csv")
glimpse(trait)

# Leaf venation trait dataset imputed values
imput_final<- read_csv("data/UCBG_venation_function_imput.csv")
glimpse(imput_final)

# Join leaf venation form (HLD) datasets --------
# read in files - Hierarchical loop decomposition (HLD) 
setwd("data/HLD/") 

spp_list_form<-spp_list%>%filter(spp_code != "Auc_jap", # removing species without venation form data
                                 spp_code != "Nym")%>%
  dplyr::select(spp_code)
glimpse(spp_list_form$spp_code)

form<- list()
data_form<-list()
for (i in 1:nrow(spp_list_form)){
  file_name <- unlist(list.files(path=paste0("data/HLD/",spp_list_form$spp_code[i])))
  form[[i]]<- str_subset(file_name, "_HLD.csv")
  data_form[[i]]<-read_csv(paste0("data/HLD",spp_list_form$spp_code[i], "/",form[[i]]))
}
# combine all files into one
dt_form<-bind_rows(data_form)
glimpse(dt_form)
unique(dt_form$name) # checking if all species were imported

# adding column with spp_codes
dt2<-dt_form %>%
  mutate(spp_code = case_when(name == "UCBG_Abemex_931304_B01-1_crop_HLD" ~ 'Abe_mex',
                              name == "UCBG_Aexpun_730753_B01-2_crop_HLD" ~ 'Aex_pun',
                              name == "Copy of UCBG-Alecor200230360-B01-2_20200114_2246_01_01_01__stitch_crop_HLD" ~ 'Ale_cor',
                              name == "UCBG-Alslig830808-B01-6_crop_recrop_CNN_1" ~ 'Als_lig',
                              name == "UCBG-Anecal620082-B01-3_20200217_0148_01_01_01__stitch_crop_CNN_1" ~ 'Ane_cal',
                              name == "UCBG-Argspi20030221-B01-5_20200216_1632_01_01_01__stitch_crop_HLD" ~ 'Arg_spi',
                              name == "UCBG-Aribae20120405-B01-3_crop_HLD" ~ 'Ari_bae',
                              name == "UCBG_Asitri822026_822026_B01_1_crop_HLD" ~ 'Asi_tri',
                              name == "Copy of UCreim-sc-AspObl-1099_CNN_1" ~ 'Asp_obl',
                              name == "UCBG-AstBan-8_20200312_1033_01_01_01__crop_crop_recrop_CNN_1" ~ 'Ast_ban',
                              name == "UCBG_Barpur000000_B01_5_newcrop_recrop_CNN_1" ~ 'Bar_pur',
                              name == "UCBG-Begfor20050211-B01-1_20200201_0139_01_01_01__stitch_crop_HLD" ~ 'Beg_for',
                              name == "UCBG-Beitar20010025-B01-1_20200201_0134_01_01_01__stitch_crop_HLD" ~ 'Bei_tar',
                              name == "UCBG-Buxhar790114-B01-4_20200131_2119_01_01_01__stitch_crop_HLD" ~ 'Bux_har',
                              name == "UCBG_Calnut_891820_B01_newcrop_crop_CNN_1" ~ 'Cal_nut',
                              name == "UCBG-Calpan-6_20200131_2205_01_01_01__stitch_crop_HLD" ~ 'Cal_pan',
                              name == "UCBG-Camsin-1_20200131_2326_01_01_01__stitch_crop_HLD" ~ 'Cam_sin',
                              name == "UCBG-Cancan611429-B01-1_20200216_1637_01_01_01__stitch_crop_HLD" ~ 'Can_can',
                              name == "UCBG_Carpal970479_B01_1_crop_HLD" ~ 'Car_pal',
                              name == "UCBG-CarPen-3-redo2_20200315_0134_01_01_01__stitch_newcrop_CNN_1" ~ 'Car_pen',
                              name == "UCBG-Catedu780415-B01-5_20200131_2138_01_01_01__stitch_crop_HLD" ~ 'Cat_edu',
                              name == "UCBG-Celaus790592-B01-2_crop_HLD" ~ 'Cel_aus',
                              name == "UCBG-CelTen-1_20200312_1112_01_01_01__stitch_crop" ~ 'Cel_ten',
                              name == "UCBG-Chinit760477-B01-1_20200201_0152_01_01_01__stitch_crop_HLD" ~ 'Chi_nit',
                              name == "UCBG-Chlspi-2_20200215_1432_01_01_01__stitch_crop" ~ 'Chl_spi',
                              name == "UCBG-Cisant20010620-B01-1_20200131_2253_01_01_01__stitch_crop_HLD" ~ 'Cis_ant',
                              name == "UCBG-Cisant20010620-B01-1_20200131_2253_01_01_01__stitch_crop_HLD_part2" ~ 'Cis_ant',
                              name == "UCBG-Cispop20041097-B01-3_crop" ~ 'Cis_pop',
                              name == "UCBG-Citmon20140599-B01-2_20200217_0135_01_01_01__stitch_crop_HLD" ~ 'Cit_mon',
                              name == "UCBG_Cluala891280_B01_3_partial_crop_corner_CNN_1" ~ 'Clu_ala',
                              name == "UCBG-Cnetri20041196-B01-3_crop_HLD" ~ 'Cne_tri',
                              name == "UCBG-Comery730581-B01-7_20200131_2222_01_01_01__stitch_crop_HLD" ~ 'Com_ery',
                              name == "UCBG_Comumb941236-B01-2_crop_HLD" ~ 'Com_umb',
                              name == "UCBG-Corbud505941025-B01-4_20200131_2239_01_01_01__stitch_crop_HLD" ~ 'Cor_bud',
                              name == "UCBG-Corlae770128-B01-1_20200201_0118_01_01_01__stitch_crop_HLD" ~ 'Cor_lae',
                              name == "UCBG_Costap_20090527_B01-2_crop_CNN_1" ~ 'Cos_tap',
                              name == "UCBG-Cripat541123-B01-5_crop_HLD" ~ 'Cri_pat',
                              name == "UCBGredo-CupNud-4063" ~ 'Cup_nud',
                              name == "UCBG-Cusspi730574-B01-5_partial_crop_CNN_1" ~ 'Cus_spi',
                              name == "UCBG_Cyaful20120101_B01_1-1_crop_HLD" ~ 'Cya_ful',
                              name == "UCBG-Cycbip931315-B01-2_20200217_0237_01_01_01__stitch_crop_HLD" ~ 'Cyc_bip',
                              name == "UCBG-Daphim20020188-B01_20200217_0226_01_01_01__stitch_crop" ~ 'Dap_him',
                              name == "UCreim-sc-DelFis-1109_CNN_3" ~ 'Del_fis',
                              name == "UCBG_Desspi830792_B01_5_crop_HLD" ~ 'Des_spi',
                              name == "UCBG_Dicant941056_B01_2_crop_CNN_1" ~ 'Dic_ant',
                              name == "UCBG-Dicfeb-6_20200131_2335_01_01_01__stitch_crop_CNN_1" ~ 'Dic_feb',
                              name == "UCBG_Dicsqu20100344_B01_3.5_crop_HLD" ~ 'Dic_squ',
                              name == "UCBG-Diolyc491820-B01-7_20200201_0111_01_01_01__stitch_CNN_1" ~ 'Dio_lyc',
                              name == "UCBG-Diopol-x_20200215_1437_01_01_01__stitch_crop_HLD" ~ 'Dio_pol',
                              name == "UCBG_Dipesc560655_B01_4_crop_HLD" ~ 'Dip_esc',
                              name == "UCBG_Dirocc870079_B01-4_crop_HLD" ~ 'Dir_occ',
                              name == "UCBG_Dooasp541237_B01_4_crop_HLD" ~ 'Doo_asp',
                              name == "UCBG-Dryrem980390-B01-3_crop_CNN_1" ~ 'Dry_rem',
                              name == "UCBG_Ehracu_980746_B01-2_crop" ~ 'Ehr_acu',
                              name == "UCBG_Equtel750357_B01_3repeat_crop_recrop_CNN_1" ~ 'Equ_tel',
                              name == "UCBG-Escher040493-B01-1_20200217_0155_01_01_01__stitch_crop_HLD" ~ 'Esc_her',
                              name == "UCBG-redo-EucCri-10029_CNN_2" ~ 'Euc_cri',
                              name == "UCBG-Eurjap610733-B01-5_20200128_1052_01_01_01__stitch_crop_HLD" ~ 'Eur_jap',
                              name == "UCredo-sc-EusJap-4105_CNN_3" ~ 'Eus_jap',
                              name == "UCBG_Frason831525_B01_2_crop_HLD" ~ 'Fra_son',
                              name == "UCBG-Gomkeu20010213-B01-3_20200216_1640_01_01_01__stitch_crop" ~ 'Gom_keu',
                              name == "UCBG-Helchi20020209-B01-9_20200131_2226_01_01_01__stitch_crop_HLD" ~ 'Hel_chi',
                              name == "UCredo-sc-HibSca-3096_CNN_3" ~ 'Hib_sca',
                              name == "UCredo-IllLan-6_CNN_3" ~ 'Ill_lan',
                              name == "UCBG-Iteili860414-B01-5_20200131_2157_01_01_01__stitch_crop_HLD" ~ 'Ite_ili',
                              name == "UCBG-Kniuva340175-B01-x_20200131_2212_01_01_01__stitch_newcrop_recrop_CNN_1" ~ 'Kni_uva',
                              name == "UCBG_Lonjap000000_B01_2_crop" ~ 'Lon_jap',
                              name == "UCBG-LuzNiv-2-robot_20200305_1339_01_01_01__stitch_newcrop_recrop_CNN_1" ~ 'Luz_niv',
                              name == "UCBG-Magche-2_20200131_2302_01_01_01__stitch_crop_HLD" ~ 'Mag_che',
                              name == "UCBG-Malpal851139-B01-1_20200217_0142_01_01_01__stitch_crop" ~ 'Mal_pal',
                              name == "UCredo-sc-MelMaj-10095" ~ 'Mel_maj',
                              name == "UCBG-Metdav930209-B01_crop_HLD" ~ 'Met_dav',
                              name == "UCBG_Micmul20120510_B01_3_crop_HLD" ~ 'Mic_mul',
                              name == "UCBG-Micpla500577-B01-3_20200217_0129_01_01_01__stitch_crop_HLD" ~ 'Mic_pla',
                              name == "UCBG_Moncar750110_B01_2_crop_HLD" ~ 'Mon_car',
                              name == "UCBG-Myrrub650176-B01_20200114_1551_01_01_01__stitch_crop_HLD" ~ 'Myr_rub',
                              name == "UCBG-Nandom860621-B01-3_20200115_0447_01_01_01__stitch_crop_HLD" ~ 'Nan_dom',
                              name == "UCBG-Notfus900611-B01_20200114_1445_01_01_01__stitch_crop_HLD" ~ 'Not_fus',
                              name == "UCBG-Oleeur551110-B01-8_20200128_1103_01_01_01__stitch_crop_HLD" ~ 'Ole_eur',
                              name == "UCBG_Onosen800321_B01_2 _crop_HLD" ~ 'Ono_sen',
                              name == "UCBG_Osmber_900317_B01-3_crop_recrop_CNN_1" ~ 'Osm_ber',
                              name == "UCBG-Paebro-3_20200312_1055_01_01_01__stitch_crop" ~ 'Pae_bro',
                              name == "UCBG-Partor20030187-B01-2_20200114_2304_01_01_01__stitch_newcrop_recrop_CNN_1" ~ 'Par_tor',
                              name == "UCBG-PerAcu630851_B01_crop_CNN_1" ~ 'Per_acu',
                              name == "UCBG_Peubol361375_B01_4_crop_HLD" ~ 'Peu_bol',
                              name == "UCBG-Philat7901174-B01-1_20200216_1645_01_01_01__stitch_crop_HLD" ~ 'Phi_lat',
                              name == "UCreim-sc-PhiVar-2098_CNN_3" ~ 'Phi_var',
                              name == "UCBGredo-PhyPur-2_CNN_1" ~ 'Phy_pur',
                              name == "UCBG-Pitvir810591-B01-12_20200128_1044_01_01_01__stitch_crop_HLD" ~ 'Pit_vir',
                              name == "UCBG_Plaori701760216_B01_x_crop_HLD" ~ 'Pla_ori',
                              name == "UCreim.sc-PolCor094_CNN_2" ~ 'Pol_cor',
                              name == "UCBG-Polcri20030330-B01-3_20200131_2243_01_01_01__stitch_crop_HLD" ~ 'Pol_cri',
                              name == "UCBG-PolVir-1_20200312_1107_01_01_01__stitch_crop_HLD" ~ 'Pol_vir',
                              name == "UCBG_Prohoo651083_B01_1_crop_HLD" ~ 'Pro_hoo',
                              name == "UCreim-sc-PteEsc-4108_leaflet_CNN_1" ~ 'Pte_esc',
                              name == "UCBG_Puyalp560083_B01_1_newcrop_CNN_1" ~ 'Puy_alp',
                              name == "UCBG_Quisap361419_B01_2_crop_HLD" ~ 'Qui_sap',
                              name == "UCBG_Ranlae910152_B01_4_crop_HLD" ~ 'Ran_lae',
                              name == "UCredo-ResLut-3" ~ 'Res_lut',
                              name == "UCBG_Rhapir671218_B01_1_crop_HLD" ~ 'Rha_pir',
                              name == "UCredo-sc-RibVib-6103_CNN_1" ~ 'Rib_vib',
                              name == "UCBG_Romcou501640-B01_crop_HLD" ~ 'Rom_cou',
                              name == "UCBG_Rubodo760803_B01_2_crop_HLD" ~ 'Rub_odo',
                              name == "UCreim-sc-RubTin-6111" ~ 'Rub_tin',
                              name == "UCBGredo-SagLat-4062_CNN_2" ~ 'Sag_lat',
                              name == "UCBG_Salsco770174_B01_1_crop_HLD" ~ 'Sal_sco',
                              name == "UCBG-Samnig-1_20200215_1426_01_01_01__stitch_crop_HLD" ~ 'Sam_nig',
                              name == "UCBG_Saucer780467_B01_1_crop_HLD" ~ 'Sau_cer',
                              name == "UCBG_Saumad900440_B01_2_crop_HLD" ~ 'Sau_mad',
                              name == "UCBG-SchPol-2-2400dpi026_CNN_1" ~ 'Sch_pol',
                              name == "UCBG_Simchi610369_B01_3_crop_HLD" ~ 'Sim_chi',
                              name == "UCBG-SmiHis-1_20200312_1119_01_01_01__stitch_crop_HLD" ~ 'Smi_his',
                              name == "UCBG_Sollig830264_B01_1_crop_HLD" ~ 'Sol_lig',
                              name == "UCBG_Stachi80116_B01_1(1)_crop_HLD" ~ 'Sta_chi',
                              name == "UCBG_Taslan600052_B01_7_crop_HLD" ~ 'Tas_lan',
                              name == "UCBG-Tecspe20020957-B01-3_20200114_2226_01_01_01__stitch_crop_HLD" ~ 'Tec_spe',
                              name == "UCBG_Todbar590561_B01_1_crop_HLD" ~ 'Tod_bar',
                              name == "UCBG_Trilan820779_B01_x_crop_recrop_CNN_1" ~ 'Tri_lan',
                              name == "UCBG-TroAra-6_20200312_1026_01_01_01__stitch_crop_CNN_3" ~ 'Tro_ara',
                              name == "UCBG_Urecar762135_B01_5_crop_HLD" ~ 'Ure_car',
                              name == "UCBG-Viogla530299-B01-1_20200216_1628_01_01_01__stitch_crop_HLD" ~ 'Vio_gla',
                              TRUE ~ 'NA'))


glimpse(dt2)
unique(dt2$spp_code)# checking if all species were correctly renamed

# write joined leaf venation form dataset
write_csv(dt2, "data/UCBG_venation_form_orig.csv")

dt_form <- dt2

# Figure out how many scales we have sufficient data at --------
# [Figure S2]
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

# save Figure S2
ggsave(g_data_available, file='figures/FigureS2.png',width=8,height=5,dpi = 300)



# Filter out rmin out of the range -----------

dt_form_filter<- dt_form%>%
  filter(widthThreshold < 0.5)%>% # filtering out rmin > 0.5 mm
  filter (widthThreshold > 0.01) # filtering out rmin < 0.01 mm
glimpse(dt_form_filter)

# Plot leaf venation form traits across vein sizes ---------------------

# filter out species without venation form data
spp_list_form<- spp_list%>%
  filter (spp_code != "Auc_jap",
          spp_code != "Nym")
glimpse(spp_list_form)

# add clade to the venation form dataset
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
ggsave(g_vd_clade, file='figures/FigureS4.png', width = 50, height = 40, units = "cm", dpi = 300, limitsize = FALSE)

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
ggsave(g_mst_clade, file='figures/FigureS5.png', width = 50, height = 40, units = "cm", dpi = 300, limitsize = FALSE)

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
ggsave(g_er_clade, file='figures/FigureS6.png', width = 50, height = 40, units = "cm", dpi = 300, limitsize = FALSE)

# Prepare data for GBM analysis ------------------
# Standardize rmin values ------------------

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

# divide rmin values into 3 bins (minor, medium, major veins) --------------

# divide data into 3 bins
width_bins <- c(0.3,0.6,0.9) 
dt_form_std_bin3 <- dt_form_std %>%
  mutate(rmin_binned = case_when(rminstd <= 0.3 ~ 0.3, # minor veins
                                 rminstd > 0.3 & rminstd <= 0.6 ~ 0.6, # medium veins
                                 rminstd > 0.6 ~ 0.9)) # major veins

glimpse(dt_form_std_bin3)

# calculate median values of venation form traits for each bin

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

# write venation form data 3 bins
write_csv(dt_form_median,"data/UCBG_form_binned3_std_median.csv")

# Join leaf venation form and function datasets ----------------------------

# re-organize venation data so that each species has a single row
dt_form_final <- dt_form_median %>%
  mutate(rmin_binned = str_c("rmin_binned", rmin_binned))%>%
  pivot_wider(names_from = rmin_binned, values_from = c("VD.median", "MST.median", "ER.median", "CR.median"), names_sep="")%>%
  rename(VD_minor = VD.medianrmin_binned0.3,
         VD_medium = VD.medianrmin_binned0.6,
         VD_major = VD.medianrmin_binned0.9,
         MST_minor = MST.medianrmin_binned0.3,
         MST_medium = MST.medianrmin_binned0.6,
         MST_major = MST.medianrmin_binned0.9,
         ER_minor = ER.medianrmin_binned0.3,
         ER_medium = ER.medianrmin_binned0.6,
         ER_major = ER.medianrmin_binned0.9,
         CR_minor = CR.medianrmin_binned0.3,
         CR_medium = CR.medianrmin_binned0.6,
         CR_major = CR.medianrmin_binned0.9)

dt_form_final <- dt_form_final %>% replace(is.na(.), 0)

glimpse(dt_form_final)
unique(dt_form_final$spp_code) # checking number of spp

# merge form and function datasets for GBM analysis

# non-imputed functional trait data
dt_all <- as.data.frame(left_join(dt_form_final, trait, by = "spp_code"))%>%
  mutate(clade = clade, Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>%
  dplyr::select(spp_code, spp_name_WP, family_WP, order_apg4, clade,
                VD_minor, VD_medium, VD_major,
                MST_minor, MST_medium, MST_major,
                ER_minor, ER_medium, ER_major,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI)
glimpse(dt_all)

data_for_gbm <- dt_all%>%
  # coding clade as numeric to make models interpretation easier
  mutate(clade2 = case_when(clade == "ferns" ~ 1, 
                            clade == "basal angiosperms" ~ 2,
                            clade == "monocots" ~ 3,
                            clade == "basal eudicots" ~ 4,
                            clade == "rosids" ~ 5,
                            clade == "asterids" ~ 6,
                            TRUE ~ 0))
# write data for GBM analysis
# non-imputed trait values
write.csv(dt_all,"data/data_for_gbm.csv")

# imputed functional trait data
dt_all_imputed <- as.data.frame(left_join(dt_form_final, imput_final, by = "spp_code"))%>%
  mutate(clade = clade, Dkleaf_mean = (Dkleaf_M + Dkleaf_L)/2)%>%
  dplyr::select(spp_code, spp_name_WP, family_WP, order_apg4, clade,
                VD_minor, VD_medium, VD_major,
                MST_minor, MST_medium, MST_major,
                ER_minor, ER_medium, ER_major,
                LMA, SWP_M, SWP_L, SWS_M, SWS_L, e_W, e_L, Dkleaf_M, Dkleaf_L, Dkleaf_mean,
                Phe, Kleaf_max, P50, P88, ISI)

glimpse(dt_all_imputed)

data_for_gbm_imput <- dt_all_imputed%>%
  mutate(clade2 = case_when(clade == "ferns" ~ 1, # coding clade as numeric to make models interpretation easier
                            clade == "basal angiosperms" ~ 2,
                            clade == "monocots" ~ 3,
                            clade == "basal eudicots" ~ 4,
                            clade == "rosids" ~ 5,
                            clade == "asterids" ~ 6,
                            TRUE ~ 0))

glimpse(data_for_gbm_imput)

# write data for GBM analysis
# imputed trait values
write_csv(data_for_gbm_imput, "data/data_for_gbm_imput.csv" )

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