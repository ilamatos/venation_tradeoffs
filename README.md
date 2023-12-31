<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Logo.png" alt="Logo" width="500" height="500">
  </a>

<h3 align="center">Leaf venation form versus function trade-offs</h3>

  <p align="center">
   Data and Rcode to reproduce analysis of the manuscript entitled "Leaf venation network architecture coordinates functional trade-offs across vein spatial scales"
    <br />
    <a href="https://github.com/ilamatos/xylem_implosion_safety"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/ilamatos/xylem_implosion_safety">View Demo</a>
    ·
    <a href="https://github.com/ilamatos/xylem_implosion_safety/issues">Report Bug</a>
    ·
    <a href="https://github.com/ilamatos/xylem_implosion_safety/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
    <li>
      <a href="#about-the-project">About the project </a>
      </ul>
   <li>
      <a href="#statistical-analysis">Statistical Analysis</a>
    </ul>
   <li>
      <a href="#leaf-venation-form-traits">Leaf traits</a>
    </ul>
    <li>
      <a href="#getting-started">Getting Started</a>
      </ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    <li>
      <a href="#PCA">Principal component analysis (PCA)</a>
      </ul>
    <li>
      <a href="#GBM">Gradient boosting regression models (GBM)</a>
      </ul>
    </li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#references">References</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project
Leaf venation networks may contribute to different functional axes: 

* []()(1) Damage resistance to drought: leaf ability to avoid water flow interruption due to xylem conduit implosion or due to the formation of air bubbles inside the conduits;
* []()(2) Damage resistance to herbivory: leaf ability to avoid water flow disruption caused by herbivores cutting across the veins; 
* []()(3) Damage resilience to drought and herbivory: leaf capacity to maintain flow after drought- or herbivory-mediated damage has occurred;
* []()(4) Flow efficiency: how efficiently water flows through the leaf;
* []()(5) Mechanical support: leaf capacity to remain upright in space.

All those functions must be traded off against the construction cost of the leaf. Due to biophysical and physiological constraints it might be impossible to construct a leaf network that simultaneously maximizes all those functions, while also minimizing costs (Blonder et al., 2020, 2018). Thus, plants may have found different solutions (i.e. different network forms) to maximize some functions to the detriment of others, depending on the selective forces under which they have evolved.

Features not related to the venation network (e.g. mesophyll conductance, chemical defenses) also influence leaf functions/costs, and may covary, cancel out, or reinforce the venation form versus function trade-offs. Quantifying the contribution of network architecture features to different leaf functions is therefore essential to understand the biophysical constraints that drove the evolution of plants with diverse network architectures. 

In this study, we obtained leaf venation architecture and functional traits from a phylogenetic diverse set of 122 plant species (in 112 families, including ferns, basal angiosperms, monocots, basal eudicots, asterids, and rosids), which allowed a deeper investigation of how venation architecture versus function trade-offs vary across taxa. Our functional dataset included ten traits (see below) and provided a more complete description of the different functional/cost axes. Our architectural multiscale dataset was extracted from fully segmented leaf images, including all veins between 10 and 500 𝜇m, and supported an accurate representation of how key venation features (see below), such as vein density, vein ramification (branching vs. looping) and shape of loops (circular vs. elongated), vary across the entire venation network. By combining those two datasets we were able to (1) describe the leaf function versus architecture trade-offs, and how they vary across plant clades and vein spatial scales; (2) provide an upper bound estimate of the contribution of venation architecture to each functional axis; (3) identify which venation architecture traits at which scale (from a continuum from minor to major veins) predict each leaf function; and (4) determine how venation architecture traits across spatial scales interact (via independence and/or coordination) to regulate leaf function.

Our dataset also allowed a robust assessment of the following specific hypotheses (see Fig.1) that have previously been tested with smaller datasets or not at all: 

* []()(H1a) Resistance to drought (Fig.1a) should be higher in networks with lower density of lager veins , as larger veins seem to be more prone to embolisms (i.e. formation and propagation of air bubbles inside the xylem conduits; Brodribb et al., 2016) and implosion (i.e. the inward collapse of a xylem conduit cell walls) under drought (Blackman et al., 2010). (H1b) Tree-like networks (i.e. branching networks) at small scales should be more resistant, because highly connected veins (i.e. looping networks) have more entry points for air bubbles to form and for embolisms to spread (Mrad et al., 2018). (H1c) Networks with more  circular loops in small veins should be more resistant to drought, because shorter conduits have less connections, thereby a lower likelihood of being in contact with an already embolized neighboring conduit from which embolisms could spread (Loepfe et al., 2007).
* []()(H2a) Resistance to herbivory (Fig.1b) should be higher in networks with higher vein density (Sack and Scoffoni, 2013), particularly of larger veins (Sack et al., 2008). This is because veins are usually harder to break than other non-lignified leaf tissues (Kawai and Okada, 2016; Kitajima and Poorter, 2010), and larger veins are usually more sclerified, making them stronger than smaller veins (Roth-Nebelsick et al., 2001; Choong, 1996). (H2b) We also expected the presence of more circular loops at all scales to increase damage resistance, as they can provide more ways to stop the propagation of mechanical fractures, especially at leaf edges (Fiorin et al., 2016; Niklas, 1999). (H2c) High secondary chemistry investment (e.g. high phenol content) may offset investment in networks with high vein density and more looping.
* []()(H3a) Resilience (Fig. 1c) to both drought and herbivory should be higher in networks with higher density of either small (Duarte et al., 2023) or large veins (Sack and Holbrook, 2006), (H3b) more looping veins at all scales (Katifori et al., 2010; Roth-Nebelsick et al., 2001), (H3c) and/or in networks with more than one primary vein emanating from the leaf base (i.e. palmate and parallel venation; Sack et al., 2008), as those features can provide redundant pathways that enable continued flow after damage (Sack and Scoffoni, 2013). However, very high redundancy, especially in larger veins, might actually decrease resilience, by facilitating the spread of embolisms (Mrad et al., 2021; Loepfe et al., 2007) and diseases (Chatelet et al., 2006).
* []()(H4a) Flow efficiency (Fig.1d) should be higher in venation networks with a higher density of branching small veins, as small veins contribute relatively more to leaf hydraulic conductance (Kawai and Okada, 2016; McKown et al., 2010; Sack and Holbrook, 2006) and tree-like networks distribute resources more efficiently than networks with loops, in the absence of damage (Fiorin et al., 2016; Katifori et al., 2010).
* []()(H5a) Mechanical support (Fig. 1e) should be higher in networks with higher vein density, particularly large veins (Sack and Scoffoni, 2013) , as they disproportionately increase leaf stiffness (Onoda et al., 2015; Choong et al., 1996, 1992; Lucas and Pereira, 1990). (H5b) The presence of looping small veins (e.g. the small transverse veins in grass leaves), may also increase stiffness by providing reinforcing cross-linkages against bending forces (Niklas, 1999).
* []()(H6a) Construction cost (Fig.1f) should be higher in networks with higher vein density (Kawai and Okada, 2016), and this effect should be strongest in networks with more large veins (John et al., 2017; Sack and Scoffoni, 2013; Niinemets et al., 2007) due to the disproportionate scaling of cost with vein size (i.e. cost is proportional to vein radius2).

Alternatively, (H7) architectural versus functional relationships could be weak, indicating that not all aspects of the network architecture have immediate functional linkages (Blonder et al., 2020), and that non-venation features play a greater role in determining leaf function.

<!-- FIGURE 1 -->
<br />
<div align="left">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Figure 2.png" alt="Logo" width="800" height="900">
  </a>

<h3 align="left">Figure 1</h3>
 (a) Radar chart for hypothesized combinations of leaf venation architecture traits (vein density - VD, minimum spanning tree ratio - MST, and loop elongation ratio - ER) at three vein spatial scales (small, medium, large) that would evolve if each leaf functional axis were independently optimized. (b) Damage resistance to drought should be higher in networks with lower density of large veins (low VDlarge), branching large and medium veins (low MSTlarge and MSTmedium) and and more circular loops in small veins (low ERsmall); (c) Damage resistance to herbivory should be higher in networks with higher large vein density (high VDlarge) and more circular loops at all scales (low ER); (d) Damage resilience to drought and herbivory should be higher in networks with higher density of large and small veins (high VDsmall and VDlarge), palmate venation (more than one midrib), and more loops (low MST) at all scales; (e) Flow efficiency should be higher in networks with  higher density (high VDsmall) of branching (high MSTsmall) small veins; (f) Mechanical support should be higher in networks with higher large vein density (high VDlarge) and more loops in small veins (low MSTsmall); and (g) Construction cost should higher in networks with higher density (high VDlarge and VDmedium) of large and medium veins. 


<p align="left">(<a href="#readme-top">back to top</a>)</p>

<!--LEAF TRAITS -->
## Leaf venation functional traits
- Flow efficiency:
 ```sh
- Kleaf max (mmol m-2 s-1 MPa-1): Maximum leaf hydraulic conductance. Describes how much water flows across the leaf in response to a water potential gradient between the leaf and the surrounding atmosphere. Higher Kleaf max indicates a higher flow efficiency.                          
 ```
- Damage resistance - drought:
 ```sh
- P50 (-MPa): Leaf water potential inducing 50% loss of Kleaf max. Lower (more negative) P50 indicates a higher damage resistance to drought.
- P88 (-MPa): Leaf water potential inducing 88% loss of Kleaf max. Lower (more negative) P88 indicates a higher damage resistance to drought.
- ISI (dimensionless): Xylem conduit implosion safety index. Measures the xylem cell walls resistance to implosion during drought. The thicker the double cell wall relative to its maximum diameter (higher ISI), the greater the resistance to drought.                                          
```
- Damage resistance - herbivory:
 ```sh
 - SWP (kJ m-2 m-1): Specific work to punch. Measures the absolute amount of work done to punch a leaf through its midrib (SWP midrib) or lamina (SWP lamina) per unit of leaf width. Higher SWP indicates a higher mechanical resistance against herbivores.
 - SWS (J m-2): Specific work to shear. Measures the absolute amount of work done to shear the leaf midrib (SWS midrib) or lamina (SWS lamina) per unit of leaf width. Higher SWS indicates a higher mechanical resistance against herbivores.
- Phe (g g-1): total phenol content in dried leaves quantified using the Folin–Ciocalteu (F–C) assay. Phe can be used as a partial proxy (as there are other potential leaf chemical defense compounds besides phenolic ones) of chemical defense against herbivores. Higher Phe values indicate higher damage resistance.                                                                                                                   
```
- Damage resilience - drought and herbivory:
 ```sh
 - ∆Kleafmean (%): Average change in Kleaf after simulated herbivory in the leaf midrib and lamina. The percentage decrease (∆Kleaf <0) or increase (∆Kleaf >0) in the leaf hydraulic conductance 48 hours after the leaf midrib and lamina has been physically damaged. Higher ∆Kleaf mean indicates a higher resilience.                       
 ```
- Mechanical support:
 ```sh
-  ε (MN m-2): Leaf flexural modulus of elasticity. Measures the whole  leaf  (ε whole) or leaf lamina (ε lamina) resistance to deformation from bending forces. Higher ε (stiffer leaves) indicates a higher mechanical support.                                                      
 ```
- Construction cost:
 ```sh
- LMA (kg m-2): Leaf mass per area. Describes the amount of resources invested on the construction of each unit of leaf area (Pérez-Harguindeguy et al., 2016). Higher LMA indicates a higher construction cost.                                                                                               
 ```    
## Leaf venation form traits

<!-- FIGURE 2 -->
<br />
<div align="left">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Table2.png" alt="Figure2" width="1000" height="700">
  </a>

<!-- STATISTICAL ANALYSIS -->
## Statistical analysis

To describe how the leaf function versus architecture trade-offs vary across clades and vein spatial scales, we carried out a PCA with all leaf traits (architecture and functional traits) to evaluate whether species in different clades and/or veins at different spatial scales occupy different portions of the architecture-function space. Prior to the PCA we imputed missing values in the functional dataset using a Bayesian hierarchical probabilistic matrix factorization (BHPMF; Schrodt et al., 2015). BHPMF imputes values based on both the phylogenetic tree (taxonomic hierarchy) and the correlation structure within the matrix of species trait values. For more details of BHPMF imputation and validation, see Appendix S1. The multiscale venation statistics (VD, ER, MST) were binned into 50 rmin bins spanning 0.01 mm (rmin < 0.01 mm, veins too small to be distinguishable in our cleared leaf images) to 0.5 mm (rmin > 0.5 mm, too few veins sampled at this size, see Fig. S1). After the imputation and binning, all trait values were centered and scaled (z-transformed) to improve comparability among them and to reduce bias towards traits with higher variance. VD and ER values were also square-root-transformed to improve normality. Then, we ran the PCA and used the broken-stick method to determine how many principal components to retain. We visualized the retained principal component scores using 95% confidence ellipses at each clade and also at each vein spatial scale (i.e. at each  rmin bin).

To quantify the direct and independent contribution of venation architecture traits to each functional axis, we fitted GBM models. GBM is a machine-learning ensemble method that can effectively capture complex non-linear interactions between predictor variables (Natekin and Knoll, 2013). We fitted our models using one functional trait (P50, P88, ISI, SWPmidrib, SWPlamina, SWSmidrib, SWSlamina, ∆Kleafmean, Kleafmax, 𝜀whole, 𝜀lamina, LMA, Phe) at a time as the response variable, and clade (ferns, basal angiosperms, monocots, basal eudicots, rosids, asterids) plus the three venation architecture traits (VD, MST, ER) as the predictor variables. To make an interpretable assessment of the contribution of venation architecture traits at each spatial scale, we binned VD, MST, and ER at three scales.

Because the range of vein sizes vary across leaves, we used two complementary approaches (scaled and unscaled  rmin) to classify veins into size categories (Fig. S2), allowing us to investigate how both relative (scaled  rmin) and absolute veins sizes (unscaled rmin ) influence the architectural-functional trade-offs. In the first approach (scaled  rmin), we scaled rmin within each species (by dividing each rmin value by the maximum rmin for each species), so that rmin across all species varied from 0 to 1. After this standardization, we classified rmin values into small (0 < rmin ≤ 0.3), medium  (0.3 < rmin ≤ 0.6), and large veins (rmin > 0.6), and calculated the median VD, MST and ER for each bin. As the venation architecture traits were extracted from whole-leaf networks we considered any absences of observed veins at a certain scale as true absences (Blonder et al., 2020), and assigned a median value of zero in those cases. In the second approach (unscaled rmin ), we log-transformed  rmin values and added 2.5 to the log-transformed values to obtain positive values of rmin. Next, we classified veins into small (0.5 < rmin  ≤ 1.0), medium (1.0 < rmin  ≤ 1.5), and large (rmin  > 1.5) bins, and then calculated the median VD, MST and ER for each bin.

To fit each GBM, we split data 80%/20% between training and test sets. Then, we used the h2o.automl as implemented in the h2o R-package to perform a hyperparameter search over the GBM parameters in order to obtain the best model for each functional trait. To prevent model overfitting, hyperparameter tuning was done with a maximum running time of 30 seconds and a 3-fold cross validation. Model performance was assessed using the Root Mean-Square Error (RMSE). From each best model, we obtained the total variance explained, as an estimate of the contribution of venation architecture traits to each functional axis. Due to the lower sample size (N = 32), GBM models for P50 and P88 were fitted using the BHPMF-imputed data, for all other response variables the original (non-imputed) functional trait dataset was used.

To identify which venation architecture traits at which scale predict the different leaf functions, we obtained the influence value of each predictor variable in each best GBM model. Variable influence was obtained using both permutation variable importance and magnitude of variable attributions (SHapley Additive exPlanations -SHAP values, Štrumbelj and Kononenko, 2014). Variable importance, ranging from 0 (lowest importance) to 1 (highest importance), was determined using the h2o.varimp function, which measures the increase in the prediction error of the model (RMSE) after variable values are permuted. SHAP values were determined using h2o.predict_contributions function, which measures the impact of every variable on the prediction by the model for each specific instance of the data, so that predictor variables with large absolute SHAP values are more important to determine the response variable.

To determine how venation architecture traits across spatial scales interact to predict each leaf function, we measure the strength of pairwise interactions between predictor variables using the H-statistic (Friedman and Popescu, 2008), as implemented with the function Interaction$new from the iml (Interpretable Machine Learning) R-package. H-statistics measure how much of the variation of the predicted outcome depends on a given pairwise interaction between predictor variables, and vary from 0 (no interaction) to 1 (100% of variance is due to interactions). 

To test for differences in the leaf functional traits across plant clades, and also to test our specific hypothesis (H3c) that resilience (𝛥Kleaf) varies across venation types (parallel, palmate, pinnate), we used Kruskal-Wallis tests followed by pairwise Wilcox tests with Benjamini and Hochberg (1995) p-value adjustment method.

All analyses were carried out using the R version 4.3.1.



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

You will need R version 4.3.1 (or greater) and the following R-packages installed and loaded in your computer to run the Rcode to reproduce the analysis of this project

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/ilamatos/venation_tradeoffs.git
   ```
2. Install the necessary R-packages
   ```sh
   install.packages(c("tidyverse", "magrittr", "ggpubr", "vegan", "ggrepel", "viridis", "data.table", "iml"))
   ```
3. Some packages may need to be installed from the source
   
```sh
  # Install BHPMF R-package
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/BHPMF/BHPMF_1.0.tar.gz"
  install.packages(packageurl, repos=NULL, type="source
  
  # Install h2o R-package 
  install.packages("C:/Users/ilain/Downloads/h2o-3.42.0.3/h2o-3.42.0.3/R/h2o_3.42.0.3.tar.gz",
                   repos = NULL, type = "source")
   ```
4. Run the R-scripts following the order below
 * []()"1_prepare_functional_data.R"
 * []()"2_prepare_form_data.R"
 * []()"3_pca_analysis.R"
 * []()"4_gbm_analysis.R"


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- PCA  -->
## Principal component analysis (PCA)

Import dataset and run the PCA
```sh
# Import datasets ---------
data_for_pca<- read_csv("data/data_for_pca.csv")
glimpse(data_for_pca)

# Run PCA -----------------
pca_all <- prcomp(data_for_pca%>%select_if(.,is.numeric)%>%
                    select(-rmin_binned, -Dkleaf_M, -Dkleaf_L, -CR.median), # removing variables not for pca
                  center=TRUE, # center data
                  scale=TRUE) # z-transform data
```
Use the Broken-stick method to determine the number of principal components (PCs) to be retained

```sh
screeplot(pca_all, main = "Screeplot of Wolf Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")
```
<!-- FIGURE 2 -->
<br />
<div align="left">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Figure S7.png" alt="Broken_stick" width="1000" height="700">
  </a>

<h3 align="left">Figure 2</h3>
Eigenvalues (grey bars) for each principal component (PC) with null model values generated by the broken stick model (red broken line). The point at which the eigenvalues cross the broken stick model distribution is considered to indicate the maximum number of components to retain.

In this case, we retain the 3 first PCs.

```sh
# extract PCA scores
pca_trajectories <- data.frame(clade=as.factor(data_for_pca$clade),
                               code=as.factor(data_for_pca$spp_code), 
                               rmin=data_for_pca$rmin_binned, 
                               pca_all$x[,c(1,2,3)])
# get axes
pca_rotations = as.data.frame(pca_all$rotation[,1:2])
pca_rotations$var = row.names(pca_rotations)
pca_rotations$x0 = 0
pca_rotations$y0 = 0
row.names(pca_rotations) = NULL
scale_factor = 10
```
Visualize the retained principal component scores using 95% confidence ellipses at each clade and also at each vein spatial scale (i.e. at each  rmin bin)

```sh
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
```
<!-- FIGURE 3 -->
<br />
<div align="left">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Figure 3.png" alt="Broken_stick" width="700" height="900">
  </a>

<h3 align="left">Figure 3</h3>
First (PC1) and second (PC2) principal components of leaf venation form (VD, MST, ER, colored in red) and functional traits (Kleaf max, P50, P88, ISI, SWPmidrib, SWPlamina, SWSmidrib, SWSlamina, 𝛥Kleafmean, 𝜀whole, 𝜀lamina, LMA, Phe, colored in black) across 50 bins of vein diameter sizes (r min). Note that leaf form traits vary across r min sizes, while functional traits do not. 95% confidence ellipses enclose the data (a) at each plant phylogenetic clade (ferns, basal angiosperms, monocots, basal eudicots, rosids, asterids); and (b) at each vein diameter class (r min). Parenthetical values indicate the percentage variance explained by the first (PC1) and second (PC2) principal component axes. 

<!-- GBM  -->
## Gradient Boosting Regression Models (GBM)
Example of GBM analysis, for the leaf functional trait Kleaf max.

Start the h2o cluster and import the leaf trait dataset.
```sh
h2o.init() # start cluster

# Import dataset
venation_path <- "data/data_for_gbm_trans.csv"
venation <- h2o.uploadFile(path = venation_path)
```
Then, define the predictor and response variables.

```sh
# Select predictor (independent) and response (dependent) variables
independent <- c("clade2", "VD_minor", "VD_medium", "VD_major",
                 "ER_minor", "ER_medium", "ER_major",
                 "MST_minor", "MST_medium", "MST_major")
dependent <- "Kleaf_max"
```
Split the dataset into training and test. In this case we are using a 80:20 ratio.

```sh
# Split dataset giving the training dataset 80% of the data
v_dat <- h2o.splitFrame(venation, ratios = 0.8, seed =1)
v_train <- v_dat[[1]]
v_test <- v_dat[[2]]
```
Fit the GBM model using h2o.automl to do the grid and algorithm search, and get the best model.

```sh
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
```
We can visualize the SHAP summary plot.

```sh
# SHAP plot
shap_plot<-h2o.shap_summary_plot(top_mod, newdata = v_train)+
  labs(title = target)+
  theme(legend.position="none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14));shap_plot
```

<!-- FIGURE 5 -->
<br />
<div align="left">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Figure 5.png" alt="Broken_stick" width="700" height="900">
  </a>

<h3 align="left">Figure 5</h3>
SHapley Additive exPlanations (SHAP) summary plots showing how vein density (VD), loop elongation ratio (ER), and minimum spanning tree ratio (MST) at three spatial scales (large, medium, and small) plus clade (predictor variables) influence 13 leaf functional traits (response variables): (a) P50; (b) P88; (c) ISI; (d) SWPmidrib; (e) SWPlamina; (f) SWSmidrib; (g) SWSlamina; (h) Phe; (i) 𝛥Kleafmean; (j) Kleafmax; (k) 𝜀whole; (l) 𝜀lamina; (m) LMA. SHAP measures the impact of predictor variables on the response variable taking into account the interactions between predictors. Each point in the SHAP summary plot represents a row from the original dataset.The y-axis indicates the predictor variable names: VDsmall, VDmedium, VDlarge, MSTsmall, MSTmedium, MSTlarge, ERsmall, ERmedium, ERlarge, and plant clades, in order of importance from top to bottom. The x-axis indicates the SHAP contribution value, i.e. by how much each predictor variable increases (SHAP > 0) or decreases (SHAP < 0) the response variable. The gradient color indicates the normalized value for the predictor variables, ranging from 0 (blue) to 1 (pink). For example, in panel (c) we can see that implosion safety is higher (more positive SHAP values) in networks with higher density of small veins (higher VDsmall indicated by points with pink color). Veins were classified in large, medium, and small classes based on their relative sizes (i.e. scaled rmin).

Finally, we can use the R-package 'iml' to calculate the strength of the pairwise interactions between the predictor variables.

```sh
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
```

Below we can see the strength of pairwise interactions between the predictor variables

<!-- FIGURE 6 -->
<br />
<div align="left">
  <a href="https://github.com/ilamatos/venation_tradeoffs">
    <img src="figures/Figure 6.png" alt="Broken_stick" width="700" height="800">
  </a>

<h3 align="left">Figure 6</h3>
Strength of interactions between predictor variables (VDsmall, VDmedium, VDlarge, MSTsmall, MSTmedium, MSTlarge, ERsmall, ERmedium, ERlarge, and clade) to determine different leaf functions (response variables) when vein sizes were classified in large, medium, and small veins using (a) scaled rmin (i.e relative veins sizes) and (b) unscaled rmin (i.e. absolute vein sizes). Interaction strength (H-stastistic) ranges from 0 (no interaction) to 1 (strongest interaction). Partial dependent co-plots for the strongest (highest H-statistic) pairwise interactions for the models with (c) scaled rmin and (b) unscaled rmin .


<!-- CONTACT -->
## Contact

Ilaine Silveira Matos - ilaine.matos@gmail.com

Project Link: [https://github.com/ilamatos/venation_tradeoffs](https://github.com/ilamatos/venation_tradeoffs)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- REFERENCES -->
## References

* []()Benjamini Y and Hochberg Y (1995) Journal of the Royal Statistical Society Series B
* []()Blonder B et al (2020) New Phytologist
* []()Blonder B et al (2018) Journal of Ecology
* []()Friedman JH and Popescu BE (2008) Annals of Applied Statistics
* []()Natekin A and Knoll A (2013) Frontiers in Neurorobotics
* []()Schrodt F (2015) Global Ecology and Biogeography
* []()Štrumbelj E and Kononenko I (2014) Knowledge and Information Systems 
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ilamatos/xylem_implosion_safety.svg?style=for-the-badge
[contributors-url]: https://github.com/ilamatos/xylem_implosion_safety/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ilamatos/xylem_implosion_safety.svg?style=for-the-badge
[forks-url]: https://github.com/ilamatos/xylem_implosion_safety/network/members
[stars-shield]: https://img.shields.io/github/stars/ilamatos/xylem_implosion_safety.svg?style=for-the-badge
[stars-url]: https://github.com/ilamatos/xylem_implosion_safety/stargazers
[issues-shield]: https://img.shields.io/github/issues/ilamatos/xylem_implosion_safety.svg?style=for-the-badge
[issues-url]: https://github.com/ilamatos/xylem_implosion_safety/issues
[license-shield]: https://img.shields.io/github/license/ilamatos/xylem_implosion_safety.svg?style=for-the-badge
[license-url]: https://github.com/ilamatos/xylem_implosion_safety/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
