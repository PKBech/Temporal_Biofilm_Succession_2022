setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/") #Perbec
setwd("/Users/nnseh/Desktop/Github/JH_biofilm2022-2/") #nasuh

######################################
## PREPROCESSING OF PHYLOSEQOBJECTS ##
######################################

# Cleaning 16S, 18S and AD data.frames

##############
## Packages ##
##############

library(tidyverse)
library("devtools")
library("Biostrings")
library("phyloseq")
library(vegan)
library(venn)
library(plotly)
library(MicEco)
library(ggbeeswarm)
library(gridExtra)
library(ggpubr)
library(grid)
library(lme4)
library(emmeans)

###############
## Load data ##
###############

#load('AD_amplicons/ps.AD.wTree.RData')
load('AmpliconAnalysis/AD/ps.AD.wTree.RData') #nasuh
ps_AD_filtered <- ps_AD_reduced

load("AmpliconAnalysis/16S/ps_16S.reduced.wTree.20062022.RData") # nasuh
#load('16S amplicon/ps_16S.reduced.wTree.20062022.RData')
ps_16S_filtered <- ps.ASV.reduced.16S

#load('18S/ps_18S.reduced.wTree.20062022.RData')
load('AmpliconAnalysis/18S/ps_18S.reduced.wTree.20062022.RData') #nasuh
ps_18S_filtered <- ps_18S

###################
## Preprocessing ##
###################


# ## Remove Phaeobacter strain from AD ASVs and Seawater-T1-3 because it is an outlier ##
# ps_AD <- subset_samples(ps_AD, Subject != "JH21" & sample_names(ps_AD) != "Seawater-T1-3")
# # Remove outliers: "Succession-T1-6" "Succession-T1-8" "Succession-T2-9" "Succession-T2-7" from 16S phyloobject
outliers <- c("Succession-T1-6","Succession-T1-8","Succession-T2-9","Succession-T2-7")

ps_AD_filtered <- subset_samples(ps_AD_filtered, !sample_names(ps_AD_filtered) %in% outliers)

ps_16S_filtered <- subset_samples(ps_16S_filtered, !sample_names(ps_16S_filtered) %in% outliers)

ps_18S_filtered <- subset_samples(ps_18S_filtered, !sample_names(ps_18S_filtered) %in% outliers)


## Subset to only bio.element study ##
ps_AD_filtered_succession <- subset_samples(ps_AD_filtered, Subject == "Succession")
ps_16S_filtered_succession <- subset_samples(ps_16S_filtered, element.type == "bioelement")
ps_18S_filtered_succession <- subset_samples(ps_18S_filtered, element.type == "bioelement")
#try filter only ADs mapped to MAGs
#ADs_to_filter
#ps_AD_succession <- subset_taxa(ps_AD, AD_ASV %in% ADs_to_filter)
#colnames(tax_table(ps_AD_succession))

# Rarefraction curves

# library("vegan")
# library("MicEco")
# 
# rcurve.AD = rcurve(ps_16S_filtered_succession.2, subsamp = seq(from = 1, to = 1000, by = 2))
# 
# ggplot(rcurve.AD, aes(Reads, Richness, color = day)) +
#   theme_bw() +
#   geom_line(size = 0.8 )+
#   xlab("\n Read") +
#   ylab("Richness\n ") +
#   theme_bw(base_size = 10)+
#   theme(axis.text.x = element_text(color = "black", face = "bold"),
#         axis.text.y = element_text(color = "black", face = "bold"),
#         axis.title = element_text(color = "black", face = "bold"))


## Scalling to total sum TSS ##
#All sample types PS
ps_AD_filter_norm = transform_sample_counts(ps_AD_filtered, function(x) 100000 * x/sum(x+0.1))
ps_16S_filter_norm = transform_sample_counts(ps_16S_filtered, function(x) 100000 * x/sum(x+0.1))
ps_18S_filter_norm = transform_sample_counts(ps_18S_filtered, function(x) 100000 * x/sum(x+0.1))

#Only bioelements PS
ps_AD_succession_filter_norm = transform_sample_counts(ps_AD_filtered_succession, function(x) 100000 * x/sum(x+0.1))
ps_16S_succession_filter_norm = transform_sample_counts(ps_16S_filtered_succession, function(x) 100000 * x/sum(x+0.1))
ps_18S_succession_filter_norm = transform_sample_counts(ps_18S_filtered_succession, function(x) 100000 * x/sum(x+0.1))


## Make phyloseq object to matrix ## 
#AD - all sample types PS
AD_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_filter_norm), "matrix")))
AD_filter_norm_dat <-  as.data.frame(round(t(AD_filter_norm_dat)))
AD_filter_norm_dat_clean <- AD_filter_norm_dat[rowSums(AD_filter_norm_dat)!=0,] # remove empty rows
str(AD_filter_norm_dat)
#AD - Only bioelements PS
AD_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_succession_filter_norm), "matrix")))
AD_succession_filter_norm_dat <-as.data.frame(round(t(AD_succession_filter_norm_dat)))
AD_succession_filter_norm_dat_clean <- AD_succession_filter_norm_dat[rowSums(AD_succession_filter_norm_dat)!=0,] #remove empty rows
str(AD_filter_norm_dat)

#16S- all sample types PS
P16S_filter_norm_dat <- as.data.frame((as(otu_table(ps_16S_filter_norm), "matrix")))
P16S_filter_norm_dat <- as.data.frame(round(t(P16S_filter_norm_dat)))
P16S_filter_norm_dat_clean <- P16S_filter_norm_dat[rowSums(P16S_filter_norm_dat)!=0,] #remove empty rows
str(P16S_filter_norm_dat)
#16S - Only bioelements PS
P16S_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_16S_succession_filter_norm), "matrix")))
P16S_succession_filter_norm_dat <- as.data.frame(round(t(P16S_succession_filter_norm_dat)))
P16S_succession_filter_norm_dat_clean <- P16S_succession_filter_norm_dat[rowSums(P16S_succession_filter_norm_dat)!=0,] #remove empty rows
str(P16S_succession_filter_norm_dat)

#18S- all sample types PS
P18S_filter_norm_dat <- as.data.frame((as(otu_table(ps_18S_filter_norm), "matrix")))
P18S_filter_norm_dat <- as.data.frame(round(t(P18S_filter_norm_dat)))
P18S_filter_norm_dat_clean <- P18S_filter_norm_dat[rowSums(P18S_filter_norm_dat)!=0,] #remove empty rows
str(P18S_filter_norm_dat)
#18S - Only bioelements PS
P18S_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_18S_succession_filter_norm), "matrix")))
P18S_succession_filter_norm_dat <- as.data.frame(round(t(P18S_succession_filter_norm_dat)))
P18S_succession_filter_norm_dat_clean <- P18S_succession_filter_norm_dat[rowSums(P18S_succession_filter_norm_dat)!=0,] # remove empty rows
str(P18S_succession_filter_norm_dat)

## Match metadata table with rownames of XX_filter_norm_dat_clean ##
#AD - all sample types PS
AD_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_filter_norm))
AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[rownames(AD_filter_norm_dat_clean), ]
#AD - Only bioelements PS
AD_succession_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_succession_filter_norm))
AD_succession_filter_norm_meta_dat <- AD_succession_filter_norm_meta_dat[rownames(AD_succession_filter_norm_dat_clean), ]

#16S- all sample types PS
P16S_filter_norm_meta_dat <- data.frame(sample_data(ps_16S_filter_norm))
P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat[rownames(P16S_filter_norm_dat_clean), ]
#16S - Only bioelements PS
P16S_succession_filter_norm_meta_dat <- data.frame(sample_data(ps_16S_succession_filter_norm))
P16S_succession_filter_norm_meta_dat <- P16S_succession_filter_norm_meta_dat[rownames(P16S_succession_filter_norm_dat_clean), ]

#18S- all sample types PS
P18S_filter_norm_meta_dat <- data.frame(sample_data(ps_18S_filter_norm))
P18S_filter_norm_meta_dat <- P18S_filter_norm_meta_dat[rownames(P18S_filter_norm_dat_clean), ]
#16S - Only bioelements PS
P18S_succession_filter_norm_meta_dat <- data.frame(sample_data(ps_18S_succession_filter_norm))
P18S_succession_filter_norm_meta_dat <- P18S_succession_filter_norm_meta_dat[rownames(P18S_succession_filter_norm_dat_clean), ]


##############
## Richness ##
##############

## Estimate Richness ##
#16S - All sample types (use PS that us TTS)

Observed_16S_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_16S_filter_norm)) ))))) 
shannon_16S_filter_norm <- diversity(round(t(otu_table(ps_16S_filter_norm))), index = "shannon")
Observed_16S_filter_norm_tab <- data.frame(cbind(Obs = Observed_16S_filter_norm$S.obs, shannon = shannon_16S_filter_norm,
                                                 Sample=rownames(Observed_16S_filter_norm) ))
colnames(Observed_16S_filter_norm_tab)[1] <- "Obs_16S"
colnames(Observed_16S_filter_norm_tab)[2] <- "Shannon_16S"

str(Observed_16S_filter_norm_tab)
#AD - All sample types (use PS that us TTS)

Observed_AD_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_AD_filter_norm))) ))))
str(Observed_AD_filter_norm)

shannon_AD_filter_norm <- diversity((round(t(otu_table(ps_AD_filter_norm)))), index = "shannon")
# str(shannon_AD_filter_norm)
Observed_AD_filter_norm_tab <- data.frame(cbind(Obs = Observed_AD_filter_norm$S.obs, shannon= shannon_AD_filter_norm, 
                                                Day=sample_data(ps_AD_filter_norm)$Day, 
                                                Timepoint=sample_data(ps_AD_filter_norm)$Timepoint, 
                                                Subject=sample_data(ps_AD_filter_norm)$element.type), 
                                          Rep_tec=sample_data(ps_AD_filter_norm)$Rep, 
                                          Sample=rownames(Observed_AD_filter_norm) )
str(Observed_AD_filter_norm_tab)
colnames(Observed_AD_filter_norm_tab)[1] <- "Obs_AD"
colnames(Observed_AD_filter_norm_tab)[2] <- "Shannon_AD"
#18S - All sample types (use PS that us TTS)
set.seed(41)
Observed_18S_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_18S_filter_norm)) )))))
shannon_18S_filter_norm <- diversity(round(t(otu_table(ps_18S_filter_norm))), index = "shannon")
Observed_18S_filter_norm_tab <- data.frame(cbind(Obs = Observed_18S_filter_norm$S.obs, shannon = shannon_18S_filter_norm,
                                                 Sample=rownames(Observed_18S_filter_norm) ))
colnames(Observed_18S_filter_norm_tab)[1] <- "Obs_18S"
colnames(Observed_18S_filter_norm_tab)[2] <- "Shannon_18S"


## Merge 16S, 18S to AD tables (to normalize the complexitiy or do correlations) ##
#AD and 16S
Observed_ADand16S_filter_norm_tab <- left_join(Observed_AD_filter_norm_tab, Observed_16S_filter_norm_tab, "Sample")
Observed_ADand16S_filter_norm_tab$Obs_AD <- as.numeric(Observed_ADand16S_filter_norm_tab$Obs_AD)
Observed_ADand16S_filter_norm_tab$Obs_16S <- as.numeric(Observed_ADand16S_filter_norm_tab$Obs_16S)
Observed_ADand16S_filter_norm_tab$Shannon_AD <- as.numeric(Observed_ADand16S_filter_norm_tab$Shannon_AD)
Observed_ADand16S_filter_norm_tab$Shannon_16S <- as.numeric(Observed_ADand16S_filter_norm_tab$Shannon_16S)
#Add metadata
Observed_ADand16S_filter_norm_tab <- left_join(Observed_ADand16S_filter_norm_tab, sample_data(ps_AD_filter_norm)[,c(1,6)], "Sample")
# add 18S
Observed_ADand16Sand18S_filter_norm_tab <- left_join(Observed_ADand16S_filter_norm_tab, Observed_18S_filter_norm_tab, "Sample")
Observed_ADand16Sand18S_filter_norm_tab$Obs_18S <- as.numeric(Observed_ADand16Sand18S_filter_norm_tab$Obs_18S)
Observed_ADand16Sand18S_filter_norm_tab$Shannon_18S <- as.numeric(Observed_ADand16Sand18S_filter_norm_tab$Shannon_18S)

str(Observed_ADand16Sand18S_filter_norm_tab)

Observed_ADand16Sand18S_filter_norm_tab$Day <- as.numeric(Observed_ADand16Sand18S_filter_norm_tab$Day)

str(Observed_ADand16Sand18S_filter_norm_tab)

## Get Day to real days and not timepoints
Observed_ADand16Sand18S_filter_norm_tab$Time = 0
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==0] <- 0
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==1] <- 3
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==2] <- 7
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==3] <- 10
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==4] <- 15
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==5] <- 23
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==6] <- 29
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==7] <- 44
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==8] <- 57
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==9] <- 71
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==10] <- 85
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==11] <- 99
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==12] <- 113

str(Observed_ADand16Sand18S_filter_norm_tab)

# Diversity


# Richness with only biofilm ####
#################################

# 16S
p_Obs_16S <- Observed_ADand16Sand18S_filter_norm_tab %>% 
  filter(Subject != "Seawater" & Subject != "Bryozor" & Subject != "JH21" & Shannon_16S > 2) %>% 
  ggplot(aes(x = Time, y = Obs_16S), fill = "black", color = "black") + 
  geom_quasirandom(dodge.width=.8, cex=1, alpha = 0.6) + 
  geom_smooth(color = "black") + 
  labs(x= "\nTime (days)", y = "\n16S richness ", subtitle = "" )+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust =1)) + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Obs_16S 

p_Obs_AD <- Observed_ADand16Sand18S_filter_norm_tab %>% 
  filter(Subject != "Seawater" & Subject != "Bryozor" & Subject != "JH21" & Shannon_AD > 2) %>% 
  ggplot(aes(x = Time, y = Obs_AD), fill = "black", color = "black") + 
  geom_quasirandom(dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth(color = "black") +
  labs(x= "\nTime (days)", y = "\nAD richness ", subtitle = "")+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust =1))+
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Obs_AD

p_Obs_18S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Seawater" & Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 2) %>% 
  ggplot(aes(x = Time, y = Obs_18S), fill = "black", color = "black") + 
  geom_quasirandom(dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth(color = "black") + ylim(0,1000) +
  labs(x= "\nTime (days)", y = "\n18S richness ", subtitle = "")+
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust =1),
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Obs_18S

Fig1D = ggarrange(p_Obs_16S, p_Obs_18S, p_Obs_AD, nrow=1)
ggsave(Fig1D, file = "Fig1D.tiff", height = 2, width =7)


# Correlation with richness ####
################################

#AD vs 16S
corr_ADvs16S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(!Subject %in% c("JH21", "Bryozor")) %>% 
  ggplot(aes(y=Obs_AD, x=Obs_16S, color = Subject)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\n16S richness ", y = "\nAD richness ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,1000))

corr_ADvs16S

#AD vs 18S
corr_ADvs18S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(!Subject %in% c("JH21", "Bryozor")) %>%   
  ggplot(aes(y=Obs_AD, x=Obs_18S, color = Subject)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\n18S richness ", y = "\nAD richness ", subtitle ="")+
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_y_continuous(limits = c(0,1200))
corr_ADvs18S

#18S vs 16S
corr_16Svs18S <- Observed_ADand16Sand18S_filter_norm_tab %>% filter(!Subject %in% c("JH21", "Bryozor")) %>%   
  ggplot(aes(y=Obs_18S, x=Obs_16S, color = Subject)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\n16S richness ", y = "\n18S richness ", subtitle ="")+
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_y_continuous(limits = c(0,1000))

corr_16Svs18S


# Richness all environments ####
################################

# 16S
p_Obs_16S.all <- Observed_ADand16Sand18S_filter_norm_tab %>% 
  filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_16S > 2) %>% 
  ggplot(aes(x = Time, y = Obs_16S, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6) + 
  geom_smooth() + 
  labs(x= "\n\nTime (days)", y = "\n16S richness ", subtitle = "" )+
  theme_bw(base_size = 8)+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Obs_16S.all 

p_Obs_AD.all <- Observed_ADand16Sand18S_filter_norm_tab %>% 
  filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_AD > 2) %>% 
  ggplot(aes(x = Time, y = Obs_AD, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = "\nAD richness", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_x_continuous(breaks = seq(0,120,by=10))
p_Obs_AD.all

p_Obs_18S.all = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 2) %>% 
  ggplot(aes(x = Time, y = Obs_18S, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = "\n18S richness", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Obs_18S.all



# Shannon all environments ####
###############################
p_Shannon_16S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_16S > 0) %>% 
  ggplot(aes(x = Time, y = Shannon_16S, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6) + 
  geom_smooth() + 
  labs(x= "\n\nTime (days)", y = " \n16S Shannon diversity\n", subtitle = "" )+
  theme_bw(base_size = 8)+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Shannon_16S 

p_Shannon_AD = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_AD > 0) %>% 
  ggplot(aes(x = Time, y = Shannon_AD, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = " \nAD Shannon diversity\n", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Shannon_AD

p_Shannon_18S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 0) %>% 
  ggplot(aes(x = Time, y = Shannon_18S, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() + ylim(0,6) +
  labs(x= "\n\nTime (days)", y = " \n18S Shannon diversity \n", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "bottom") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Shannon_18S

legend.env = as_ggplot(ggpubr::get_legend(p_Shannon_18S))
# Plot for supp figure
FigS1 = ggarrange(corr_ADvs16S,corr_ADvs18S,corr_16Svs18S,
                  p_Obs_16S.all, p_Obs_18S.all,p_Obs_AD.all,
                  p_Shannon_16S,p_Shannon_18S, p_Shannon_AD)


ggsave(FigS1, file = "FigS1.tiff", height = 16, width =18, units = "cm")

ggsave(legend.env, file = "legend.tiff", height = 1, width =7, units = "cm")



#Linear model tested by emmeans ------
##Subset to succession and seawater 
Observed_ADand16Sand18S_filter_norm_tab_filter <- Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 0 & Shannon_16S > 0 & Shannon_AD > 0)

str(Observed_ADand16Sand18S_filter_norm_tab_filter)
#Add replicate ID to dataframe
#Observed_ADand16Sand18S_filter_norm_tab_filter$ID=gsub("(.*)(-T[0-9])","\\1",Observed_ADand16Sand18S_filter_norm_tab_filter$Sample)

#Make biorep as variable for lmer random effect 
Observed_ADand16Sand18S_filter_norm_tab_filter = tidyr::separate(Observed_ADand16Sand18S_filter_norm_tab_filter,Sample, 
                                                                 c("System", "Timepoint", "BioRep"), remove = FALSE)
Observed_ADand16Sand18S_filter_norm_tab_filter$Cage = paste(Observed_ADand16Sand18S_filter_norm_tab_filter$System,
                                                            Observed_ADand16Sand18S_filter_norm_tab_filter$BioRep, sep= "-")
Observed_ADand16Sand18S_filter_norm_tab_filter$Cage = as.factor(Observed_ADand16Sand18S_filter_norm_tab_filter$Cage)

#lmer 
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD.cage=lmer(Obs_AD~Time*Subject + (1|Cage), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S.cage=lmer(Obs_16S~Time*Subject + (1|Cage), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S.cage=lmer(Obs_18S~Time*Subject + (1|Cage), Observed_ADand16Sand18S_filter_norm_tab_filter)
car::Anova(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD.cage)
car::Anova(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S.cage)
car::Anova(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S.cage)

#LMer for ADs 16S and 18S# Howver, this is only relevant for the succession not the seawater (no spatial effect from seawater samples)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD=lmer(Obs_AD~Time*Subject + (1|ID), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S=lmer(Obs_16S~Time*Subject + (1|ID), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S=lmer(Obs_18S~Time*Subject + (1|ID), Observed_ADand16Sand18S_filter_norm_tab_filter)

#LM only for observed and shannon
Observed_ADand16Sand18S_filter_norm_tab_filter_LM_ObsAD=lm(Obs_AD~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Obs16S=lm(Obs_16S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Obs18S=lm(Obs_18S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)

# Observed_ADand16Sand18S_filter_norm_tab_filter_LM_ShaAD=lm(Shannon_AD~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
# Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Sha16S=lm(Shannon_16S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
# Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Sha18S=lm(Shannon_18S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
# 

#significance for observed richness
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD_EM=emmeans(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD, ~Subject|Time)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD_EM_stat <- contrast(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD_EM, interaction = "pairwise", adjust = "Bonferroni")
#Result is that over time AD OBS richness is significant greater in succession compared to seawater 
# Time = 54.4:
#   Subject_pairwise      estimate   SE   df t.ratio p.value
# Seawater - Succession     -162 43.5 18.2  -3.727  0.0015
# 
# Degrees-of-freedom method: kenward-roger 


Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S_EM=emmeans(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S, ~Subject|Time)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S_EM_stat <- contrast(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S_EM, interaction = "pairwise", adjust = "Bonferroni")
#Result is that over time 16S OBS richness is significant greater in succession compared to seawater 
# Time = 54.4:
#   Subject_pairwise      estimate   SE   df t.ratio p.value
# Seawater - Succession     -988 99.9 18.1  -9.890  <.0001
# 
# Degrees-of-freedom method: kenward-roger 


Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S_EM=emmeans(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S, ~Subject|Time)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S_EM_stat <- contrast(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S_EM, interaction = "pairwise", adjust = "Bonferroni")
#Result is that over time 18S OBS richness is NOT significant diffrent between succession and seawater 
# Time = 54.4:
#   Subject_pairwise      estimate   SE   df t.ratio p.value
# Seawater - Succession    -14.6 62.7 17.3  -0.233  0.8186
# 
# Degrees-of-freedom method: kenward-roger






## Dunnett's test ----
#Make  Dunnest test with day 29 as "control"/ reference to all other points

#Make all combinations to individual groups
groups <- factor(paste(Observed_ADand16Sand18S_filter_norm_tab_filter$Subject,Observed_ADand16Sand18S_filter_norm_tab_filter$Time))



library(DescTools)
#ADs obs richness
DunnettTest_Day29_ObsAD <- DunnettTest(Observed_ADand16Sand18S_filter_norm_tab_filter$Obs_AD, groups, control = "Succession 29")
DunnettTest_Day29_ObsAD_dat <- as.data.frame(DunnettTest_Day29_ObsAD$`Succession 29`)
DunnettTest_Day29_ObsAD_dat$pval <- round(DunnettTest_Day29_ObsAD_dat$pval, digits = 3) 

#Results ADs
# diff    lwr.ci     upr.ci  pval
# Seawater 10-Succession 29    -365.66667 -634.7135  -96.61981 0.002
# Seawater 113-Succession 29   -367.66667 -636.7135  -98.61981 0.001
# Seawater 15-Succession 29    -318.00000 -587.0469  -48.95314 0.010
# Seawater 23-Succession 29    -327.33333 -596.3802  -58.28648 0.007
# Seawater 29-Succession 29    -339.50000 -654.9854  -24.01460 0.026
# Seawater 44-Succession 29    -308.33333 -577.3802  -39.28648 0.013
# Seawater 57-Succession 29    -309.66667 -578.7135  -40.61981 0.013
# Seawater 71-Succession 29    -388.66667 -657.7135 -119.61981 0.001
# Seawater 85-Succession 29    -355.66667 -624.7135  -86.61981 0.002
# Seawater 99-Succession 29    -374.00000 -643.0469 -104.95314 0.001
# Succession 10-Succession 29  -224.66667 -414.9115  -34.42181 0.010
# Succession 113-Succession 29 -140.00000 -330.2449   50.24486 0.324
# Succession 15-Succession 29   -42.66667 -232.9115  147.57819 1.000
# Succession 23-Succession 29   -25.88889 -216.1337  164.35597 1.000
# Succession 44-Succession 29  -237.55556 -427.8004  -47.31070 0.005
# Succession 57-Succession 29  -250.85714 -454.2374  -47.47684 0.006
# Succession 71-Succession 29  -277.77778 -468.0226  -87.53292 0.000
# Succession 85-Succession 29  -318.55556 -508.8004 -128.31070 0.000
# Succession 99-Succession 29  -324.00000 -520.0999 -127.90009 0.000

#16S
DunnettTest_Day29_Obs16S <- DunnettTest(Observed_ADand16Sand18S_filter_norm_tab_filter$Obs_16S, groups, control = "Succession 29")
DunnettTest_Day29_Obs16S_dat <- as.data.frame(DunnettTest_Day29_Obs16S$`Succession 29`)
DunnettTest_Day29_Obs16S_dat$pval <- round(DunnettTest_Day29_Obs16S_dat$pval, digits = 3) 

# diff      lwr.ci    upr.ci  pval
# Seawater 10-Succession 29    -1120.77778 -1728.76780 -512.7878 0.000
# Seawater 113-Succession 29   -1142.77778 -1750.76780 -534.7878 0.000
# Seawater 15-Succession 29    -1195.44444 -1803.43446 -587.4544 0.000
# Seawater 23-Succession 29    -1081.11111 -1689.10113 -473.1211 0.000
# Seawater 29-Succession 29    -1026.94444 -1739.87594 -314.0130 0.001
# Seawater 44-Succession 29    -1239.11111 -1847.10113 -631.1211 0.000
# Seawater 57-Succession 29     -872.44444 -1480.43446 -264.4544 0.001
# Seawater 71-Succession 29     -858.44444 -1466.43446 -250.4544 0.001
# Seawater 85-Succession 29    -1097.77778 -1705.76780 -489.7878 0.000
# Seawater 99-Succession 29    -1089.44444 -1697.43446 -481.4544 0.000
# Succession 10-Succession 29   -322.88889  -752.80276  107.0250 0.296
# Succession 113-Succession 29   450.11111    20.19724  880.0250 0.033
# Succession 15-Succession 29   -312.22222  -742.13609  117.6916 0.341
# Succession 23-Succession 29   -613.55556 -1043.46942 -183.6417 0.001
# Succession 44-Succession 29   -122.66667  -552.58053  307.2472 0.999
# Succession 57-Succession 29    232.12698  -227.47027  691.7242 0.832
# Succession 71-Succession 29    111.11111  -318.80276  541.0250 1.000
# Succession 85-Succession 29    -73.66667  -503.58053  356.2472 1.000
# Succession 99-Succession 29   -188.56944  -631.71451  254.5756 0.947

#And 18S 
DunnettTest_Day23_Sha18S <- DunnettTest(Observed_ADand16Sand18S_filter_norm_tab_filter$Shannon_18S, groups, control = "Succession 23")
DunnettTest_Day23_Sha18S_dat <- as.data.frame(DunnettTest_Day23_Sha18S$`Succession 23`)
DunnettTest_Day23_Sha18S_dat$pval <- round(DunnettTest_Day23_Sha18S_dat$pval, digits = 3) 
DunnettTest_Day23_Sha18S_dat

# diff    lwr.ci     upr.ci  pval
# Seawater 10-Succession 23    -278.111111 -617.6581   61.43587 0.194
# Seawater 113-Succession 23    -33.777778 -373.3248  305.76920 1.000
# Seawater 15-Succession 23    -302.777778 -642.3248   36.76920 0.117
# Seawater 23-Succession 23    -357.777778 -697.3248  -18.23080 0.032
# Seawater 29-Succession 23     -83.111111 -481.2652  315.04301 1.000
# Seawater 44-Succession 23    -268.444444 -607.9914   71.10253 0.233
# Seawater 57-Succession 23      37.888889 -301.6581  377.43587 1.000
# Seawater 71-Succession 23    -140.444444 -479.9914  199.10253 0.958
# Seawater 85-Succession 23       3.888889 -335.6581  343.43587 1.000
# Seawater 99-Succession 23    -117.444444 -456.9914  222.10253 0.992
# Succession 10-Succession 23   -61.333333 -301.4293  178.76264 1.000
# Succession 113-Succession 23   -8.888889 -248.9849  231.20708 1.000
# Succession 15-Succession 23    38.444444 -201.6515  278.54041 1.000
# Succession 29-Succession 23  -101.555556 -341.6515  138.54041 0.949
# Succession 44-Succession 23  -336.666667 -576.7626  -96.57070 0.001
# Succession 57-Succession 23  -123.682540 -380.3559  132.99085 0.873
# Succession 71-Succession 23  -317.444444 -557.5404  -77.34847 0.002
# Succession 85-Succession 23  -367.555556 -607.6515 -127.45959 0.000
# Succession 99-Succession 23  -145.861111 -393.3464  101.62415 0.647



#write.table(DunnettTes_Day0_dat, file = "DunnettTes_Day0_dat.csv", dec=",", sep="_")

## Tukey's test ####


#Perfom tukeys test on shannon richness for aov
aov_Sha18S_succession <- Observed_ADand16Sand18S_filter_norm_tab_filter %>%
  filter(Subject == "Succession") %>% aov(Shannon_18S ~ as.factor(Time), .)

TukeyHSD(aov_Sha18S_succession)





#### beta-diversity Only for bioelements ####

#Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions
AD_filter_norm_dat <- as.data.frame(t(as(otu_table(ps_AD_succession_filter_norm), "matrix")))
AD_filter_norm_dat <- round(AD_filter_norm_dat)
str(AD_filter_norm_dat)

P16S_filter_norm_dat <- as.data.frame((as(otu_table(ps_16S_succession_filter_norm), "matrix")))
P16S_filter_norm_dat <- as.data.frame(round(t(P16S_filter_norm_dat)))
str(P16S_filter_norm_dat)

P18S_filter_norm_dat <- as.data.frame((as(otu_table(ps_18S_succession_filter_norm), "matrix")))
P18S_filter_norm_dat <- as.data.frame(round(t(P18S_filter_norm_dat)))
str(P16S_filter_norm_dat)


#remove empty rows
AD_filter_norm_dat_clean <- AD_filter_norm_dat[rowSums(AD_filter_norm_dat)!=0,]
P16S_filter_norm_dat_clean <- P16S_filter_norm_dat[rowSums(P16S_filter_norm_dat)!=0,]
P18S_filter_norm_dat_clean <- P18S_filter_norm_dat[rowSums(P18S_filter_norm_dat)!=0,]

#Match metadata table with rownames of X_filter_norm_dat_clean

AD_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_filter_norm))
#AD_filter_norm_meta_dat <- AD_succession_filter_norm_meta_dat 
AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[rownames(AD_filter_norm_dat_clean), ]
#add real timepionts as numeric days to AD metadata frame
str(Observed_ADand16Sand18S_filter_norm_tab)
str(AD_filter_norm_meta_dat)

AD_filter_norm_meta_dat_joined <- left_join(AD_filter_norm_meta_dat, Observed_ADand16Sand18S_filter_norm_tab, "Sample")
str(AD_filter_norm_meta_dat_joined)

AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat_joined[,c(1:9,17)]
str(P16S_filter_norm_meta_dat)
## Add arbitary stages to all metadata 

stages <- as.data.frame(cbind(Sample = rownames(as.data.frame(sample_data(ps_16S_filtered_succession))), 
                              Stages = c(rep("Early", 4), rep("Late", 27),  rep("Early", 31), rep("Peak", 9),rep("Late", 27) ) ) )
stages
rownames(stages) <- stages$Sample
#Join to all dataframes
AD_filter_norm_meta_dat <- left_join(AD_filter_norm_meta_dat,stages, "Sample" )

P16S_filter_norm_meta_dat <- data.frame(sample_data(ps_16S_filter_norm))
P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat[rownames(P16S_filter_norm_dat_clean), ]
P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat %>% mutate(Sample = rownames(P16S_filter_norm_meta_dat))
P16S_filter_norm_meta_dat <- left_join(P16S_filter_norm_meta_dat,stages, "Sample" )


P18S_filter_norm_meta_dat <- data.frame(sample_data(ps_18S_filter_norm))
P18S_filter_norm_meta_dat <- P18S_filter_norm_meta_dat[rownames(P18S_filter_norm_dat_clean), ]
P18S_filter_norm_meta_dat <- P18S_filter_norm_meta_dat %>% mutate(Sample = rownames(P18S_filter_norm_meta_dat))
P18S_filter_norm_meta_dat <- left_join(P18S_filter_norm_meta_dat,stages, "Sample" )


# samples_to_move <-as.vector(which(is.na(rowSums(AD_filter_norm_dat))))

# AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[-samples_to_move,]
# 
# AD_filter_norm_dat <- round((AD_filter_norm_dat[-samples_to_move,]))

#AD_norm_dat_t <- t(AD_norm_dat)

# #Squar-root transform matrix
sqrt_AD_filter_norm_dat_clean= sqrt(AD_filter_norm_dat_clean)

head(sqrt_AD_filter_norm_dat_clean)
sqrt_P16S_filter_norm_dat_clean= sqrt(P16S_filter_norm_dat_clean)
head(sqrt_P16S_filter_norm_dat_clean)
sqrt_P18S_filter_norm_dat_clean= sqrt(P18S_filter_norm_dat_clean)
head(sqrt_P18S_filter_norm_dat_clean)


# sqrt_noPhaeo_genus_norm_dat_t = sqrt(noPhaeo_genus_norm_dat_t)

#Bray-curtis distance matrix
AD_dist_bray = as.matrix((vegdist((AD_filter_norm_dat_clean), "bray")))
str((AD_dist_bray))
P16S_dist_bray = as.matrix((vegdist((P16S_filter_norm_dat_clean), "bray")))
str(P16S_dist_bray)
P18S_dist_bray = as.matrix((vegdist((P18S_filter_norm_dat_clean), "bray")))
str(P18S_dist_bray)

#Beta-disper ----
## Calculate multivariate dispersions from all 18 combinations
AD_dist_bray_bdisp_stages <- betadisper(vegdist((AD_filter_norm_dat_clean), "bray"), 
                                        factor(paste(AD_filter_norm_meta_dat$Stages,sep=" ")))
#PCoA
plot(AD_dist_bray_bdisp_stages)

#Test betadisper significans
## Permutation test for F
permutest(AD_dist_bray_bdisp_stages, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD_AD <- TukeyHSD(AD_dist_bray_bdisp_stages))
plot(bdisp.HSD_AD)



str(AD_dist_bray_bdisp_stages)

P16S_dist_bray_bdisp_stages <- betadisper(vegdist((P16S_filter_norm_dat_clean), "bray"), 
                                          factor(paste(P16S_filter_norm_meta_dat$Stages,sep=" "))) # use this file for PCoA plot (ggplot)
plot(P16S_dist_bray_bdisp_stages)

#Test betadisper significans
## Permutation test for F
permutest(P16S_dist_bray_bdisp_stages, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD_16S <- TukeyHSD(P16S_dist_bray_bdisp_stages))
plot(bdisp.HSD_16S) 




P18S_dist_bray_bdisp_stages <- betadisper(vegdist((P18S_filter_norm_dat_clean), "bray"), 
                                          factor(paste(P18S_filter_norm_meta_dat$Stages,sep=" ")))
plot(P18S_dist_bray_bdisp_stages)

#Test betadisper significans
## Permutation test for F
permutest(P18S_dist_bray_bdisp_stages, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD_18S <- TukeyHSD(P18S_dist_bray_bdisp_stages))
plot(bdisp.HSD_18S)

# Correlation between betadisper and chla std for 18S #####



#PERMANOVA (Adonis) ----
str(P16S_dist_bray)
str(P18S_dist_bray)
str(AD_dist_bray)



adonis_P16S_dist_bray <- adonis2(P16S_dist_bray ~ day + Stages + biorep + temperature + salinity + chla, 
                               P16S_filter_norm_meta_dat)
adonis_P16S_dist_bray

adonis_P18S_dist_bray <- adonis2(P18S_dist_bray ~ day + Stages + biorep + temperature + salinity + chla, 
                                 P18S_filter_norm_meta_dat)
adonis_P18S_dist_bray

##Add temperature to AD metadata frame
colnames(AD_filter_norm_meta_dat)[3] <- "day"
AD_filter_norm_meta_dat <- left_join(AD_filter_norm_meta_dat, P18S_filter_norm_meta_dat, "Sample")

adonis_AD_dist_bray <-  adonis2(AD_dist_bray ~ day.x + Stages.x + biorep + temperature + salinity + chla, 
                                AD_filter_norm_meta_dat)
adonis_AD_dist_bray


# write.table(adonis_Genus_dat, file = "adonis_Genus_dat.csv", sep = "\t", dec = ",")
# write.table(noPhao_adonis_Genus_dat, file = "noPhao_adonis_Genus_dat.csv", sep = "\t", dec = ",")
# 



### Procruster analysis #####
data(varespec)
# vare.dist16S <- vegdist((P16S_filter_norm_dat_clean), "bray")
# vare.distAD <- vegdist((AD_filter_norm_dat_clean), "bray")
# mds.null_AD <- monoMDS(vare.distAD, y = cmdscale(vare.distAD))
# mds.alt_AD <- monoMDS(vare.distAD)
# vare.proc_AD <- procrustes(mds.alt_AD, mds.null_AD, symmetric = TRUE)
# vare.proc_AD
# summary(vare.proc_AD)
# plot(mds.alt_AD)
# plot(vare.proc_AD, kind=2)
# residuals(vare.proc_AD)
# vare.pca <- prcomp(AD_filter_norm_dat_clean)
# scores(prcomp(AD_filter_norm_dat_clean))
PCoA.ord_AD <- AD_dist_bray_bdisp_stages$vectors
PCoA.ord_16S <- P16S_dist_bray_bdisp_stages$vectors
PCoA.ord_18S <- P18S_dist_bray_bdisp_stages$vectors
#clean up rownames so it match the same sample names for all ordinations
PCoA.ord_16S <- PCoA.ord_16S[rownames(PCoA.ord_16S) %in% rownames(PCoA.ord_AD),]
PCoA.ord_AD <- PCoA.ord_AD[rownames(PCoA.ord_AD) %in% rownames(PCoA.ord_16S),]
PCoA.ord_18S <- PCoA.ord_18S[rownames(PCoA.ord_18S) %in% rownames(PCoA.ord_AD),]



cmdscale(AD_dist_bray)
str(PCoA.ord_AD)
cbind(PCoA.ord_16S, as.data.frame(AD_dist_bray_bdisp_stages$group[1:87]))

protest(X = PCoA.ord_16S, Y = PCoA.ord_AD, scores = "sites", permutations = 999)
#make data frame for procrustes residuals for 16S and AD procrustes
str(procrustes(X = PCoA.ord_AD, Y = PCoA.ord_16S,symmetric=TRUE), kind = 2)


residuals_dat_AD_16S <- as.data.frame(resid(procrustes(X = PCoA.ord_16S, Y = PCoA.ord_AD,symmetric=TRUE)))
residuals_dat_AD_16S <- cbind(residuals_dat_AD_16S, Sample = rownames(residuals_dat_AD_16S))
colnames(residuals_dat_AD_16S)[1] <- "Procrustes risidual"

residuals_dat_AD_16S <- left_join(residuals_dat_AD_16S, AD_filter_norm_meta_dat)
#Get the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.
quantiles_AD_16S <- as.data.frame(t(quantile(resid(procrustes(X = PCoA.ord_16S, Y = PCoA.ord_AD,symmetric=TRUE)))))
#Plot residuals 
residuals_dat_AD_16S %>%
  ggplot(aes(x = as.factor(day.x), y = `Procrustes risidual`, col = Stages.x, fill = Stages.x)) + 
  geom_quasirandom(aes(col = Stages.x), dodge.width=.8, cex=1, alpha = 0.6) + 
  #geom_smooth() + 
  labs(x= "\nTime (days)", subtitle = "16S and AD ordinations association" )+
  theme_bw() +
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  scale_color_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  geom_hline(yintercept=quantiles_AD_16S$`50%`) + 
  geom_hline(yintercept=quantiles_AD_16S$`25%`, linetype="dashed") +
  geom_hline(yintercept=quantiles_AD_16S$`75%`, linetype="dashed") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        #legend.position = "none"
  )
protest(X = PCoA.ord_18S, Y = PCoA.ord_16S, scores = "sites", permutations = 999)
#make data frame for procrustes residuals for 16S and 18S procrustes
str(procrustes(X = PCoA.ord_18S, Y = PCoA.ord_16S,symmetric=TRUE), kind = 2)


residuals_dat_18S_16S <- as.data.frame(resid(procrustes(X = PCoA.ord_18S, Y = PCoA.ord_16S,symmetric=TRUE)))
residuals_dat_18S_16S <- cbind(residuals_dat_18S_16S, Sample = rownames(residuals_dat_18S_16S))
colnames(residuals_dat_18S_16S)[1] <- "Procrustes risidual"

residuals_dat_18S_16S <- left_join(residuals_dat_18S_16S, AD_filter_norm_meta_dat)
#Get the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.
quantiles_18S_16S <- as.data.frame(t(quantile(resid(procrustes(X = PCoA.ord_18S, Y = PCoA.ord_16S,symmetric=TRUE)))))
#Plot residuals 
residuals_dat_18S_16S %>%
  ggplot(aes(x = as.factor(day.y), y = `Procrustes risidual`, col = Stages.x, fill = Stages.x)) + 
  geom_quasirandom(aes(col = Stages.x), dodge.width=.8, cex=1, alpha = 0.6) + 
  #geom_smooth() + 
  labs(x= "\nTime (days)", subtitle = "16S and 18S ordinations association" )+
  theme_bw() +
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  scale_color_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  geom_hline(yintercept=quantiles_18S_16S$`50%`) + 
  geom_hline(yintercept=quantiles_18S_16S$`25%`, linetype="dashed") +
  geom_hline(yintercept=quantiles_18S_16S$`75%`, linetype="dashed") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        #legend.position = "none"
        )


protest(X = PCoA.ord_18S, Y = PCoA.ord_AD, scores = "sites", permutations = 999)

##### ggplot2 - PCOA ####

# 16S (use data from betadisper)
centroids.16 = data.frame(P16S_dist_bray_bdisp_stages$centroids[,1:2])
centroids.16$Phase = rownames(centroids.16)
points.16 = data.frame(P16S_dist_bray_bdisp_stages$vectors[,1:2], Phase = P16S_dist_bray_bdisp_stages$group)

main_col <- c("#240785", "#f21395", "#e0b62b")

PCOA.16S = ggplot(centroids.16, aes(PCoA1, PCoA2, color = Phase)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.16, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = 0.3811575, y = -0.08426253, xend = PCoA1, yend = PCoA2), 
               data = points.16[points.16$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  -0.2435972, y = -0.00192419, xend = PCoA1, yend = PCoA2), 
               data = points.16[points.16$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_segment(aes(x =  0.1717421, y = 0.25524678, xend = PCoA1, yend = PCoA2), 
               data = points.16[points.16$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_point(size = 4)+
  stat_ellipse(data = points.16, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [8.2%]", y = "\nPCoA2 [2.7%]", title = "Prokaryotic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b"))
PCOA.16S

legend.pcoa = as_ggplot(ggpubr::get_legend(PCOA.16S))
ggsave(legend.pcoa, file = "legend.pcoa.pdf", width = 1, height = 1)

# AD (use data from betadisper)
centroids.AD = data.frame(AD_dist_bray_bdisp_stages$centroids[,1:2])
points.AD = data.frame(AD_dist_bray_bdisp_stages$vectors[,1:2], Phase = AD_dist_bray_bdisp_stages$group)

PCOA.AD = ggplot(centroids.AD, aes(PCoA1, PCoA2, color = rownames(centroids.AD))) +
  theme_bw(base_size = 12)+
  geom_point(data = points.AD, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = 0.31155964, y = -0.05353969, xend = PCoA1, yend = PCoA2), 
               data = points.AD[points.AD$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  -0.15593789, y = 0.01081140, xend = PCoA1, yend = PCoA2), 
               data = points.AD[points.AD$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_segment(aes(x = 0.05237552, y = 0.24482054, xend = PCoA1, yend = PCoA2), 
               data = points.AD[points.AD$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_point(size = 4)+
  stat_ellipse(data = points.AD, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [4.3%]", y = "\nPCoA2 [2.7%]", title = "AD")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b"))

# 18S (use data from betadisper)
centroids.18S = data.frame(P18S_dist_bray_bdisp_stages$centroids[,1:2])
points.18S = data.frame(P18S_dist_bray_bdisp_stages$vectors[,1:2], Phase = P18S_dist_bray_bdisp_stages$group)

PCOA.18S = ggplot(centroids.18S, aes(PCoA1, PCoA2, color = rownames(centroids.18S))) +
  theme_bw(base_size = 12)+
  geom_point(data = points.18S, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = -0.3665717, y = -0.083343689, xend = PCoA1, yend = PCoA2), 
               data = points.18S[points.18S$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  0.2549065, y = -0.005344822, xend = PCoA1, yend = PCoA2), 
               data = points.18S[points.18S$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_segment(aes(x = -0.1308241, y = 0.319936332, xend = PCoA1, yend = PCoA2), 
               data = points.18S[points.18S$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_point(size = 4)+
  stat_ellipse(data = points.18S, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [8.8%]", y = "\nPCoA2 [4.2%]", title = "Eukaryotic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b"))

PCOAs = grid.arrange(PCOA.16S, PCOA.18S, PCOA.AD, nrow = 1)
ggsave(PCOAs, file = 'PCOAs.tiff', width = 10.5, height = 3)

# NDMS AD for only bioelements-----
NMDS_ord_bray_AD = metaMDS(AD_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray_AD)

stressplot(NMDS_ord_bray_AD)

#build a data frame with NMDS coordinates and metadata
NMDS1_AD = NMDS_ord_bray_AD$points[,1]
NMDS2_AD = NMDS_ord_bray_AD$points[,2]
NMDS3_AD = NMDS_ord_bray_AD$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species
AD_filter_norm_meta_dat
P16S_filter_norm_meta_dat

NMDS_AD = data.frame(NMDS1 = NMDS1_AD, NMDS2 = NMDS2_AD, NMDS2 = NMDS3_AD, Time = AD_filter_norm_meta_dat_joined$Time, Chla = AD_filter_norm_meta_dat$Chla_suc, Stages = AD_filter_norm_meta_dat$Stages)


p_time_NMDS12_AD <- ggplot(NMDS_AD, aes(x=NMDS1, y=NMDS2, fill = Time)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(subtitle = "AD") +
  geom_point(shape = 21,size = 2,colour = "black") +
  scale_fill_gradient(low = "#9111ab", high = "#e3b039") +
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


p_bryo_NMDS12_AD <- ggplot(NMDS_AD, aes(x=NMDS1, y=NMDS2, fill = Stages)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(title = "")  +
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) + 
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 




# NDMS 16S for only bioelements-----
NMDS_ord_bray_16S = metaMDS(P16S_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray_16S)

stressplot(NMDS_ord_bray_16S)

#build a data frame with NMDS coordinates and metadata
NMDS116S = NMDS_ord_bray_16S$points[,1]
NMDS216S = NMDS_ord_bray_16S$points[,2]
NMDS316S = NMDS_ord_bray_16S$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species
AD_filter_norm_meta_dat
P16S_filter_norm_meta_dat

P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat %>% mutate(Rep_tec = rownames(P16S_filter_norm_meta_dat))
P16S_outliers <- c(names(NMDS116S[NMDS116S < -0.25]), names(NMDS116S[NMDS116S > 0.25]))
#Remove "Succession-T1-6" "Succession-T1-8" "Succession-T2-9" "Succession-T2-7" from phyloseq object earlier

NMDS_16S = data.frame(NMDS1 = NMDS116S, NMDS2 = NMDS216S, NMDS3 = NMDS316S, Time = P16S_filter_norm_meta_dat$day ,  Chla = P16S_filter_norm_meta_dat$chla , Stages =   P16S_filter_norm_meta_dat$Stages)


p_time_NMDS12_16s <- ggplot(NMDS_16S, aes(x=NMDS1, y=NMDS2, fill = Time)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(subtitle = "16S") +
  geom_point(shape = 21,size = 2,colour = "black") +
  scale_fill_gradient(low = "#9111ab", high = "#e3b039") +
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    axis.title.x=element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

p_bryo_NMDS12_16s <- ggplot(NMDS_16S, aes(x=NMDS1, y=NMDS2, fill = Stages)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(title = "")  +
  geom_point(shape = 21,size = 2,colour = "black") +
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw(base_size = 8) +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 

# 
# p_chla_NMDS12_16S <- ggplot(NMDS_16S, aes(x=NMDS1, y=NMDS2, fill = Chla)) +
#   #stat_ellipse() +
#   theme_bw(base_size = 8) +
#   labs(title = "")  +
#   geom_point(shape = 21,size = 2,colour = "black") +
#   geom_point(shape = 21, alpha = 0.25, size = 2)
# 
# library(gridExtra)
# grid.arrange(p_subject_NMDS12, p_time_NMDS12, p_subject_NMDS23, p_time_NMDS23, ncol=2, nrow = 2, 
#              layout_matrix = rbind(c(1,2), c(3,4)),
#              widths = c(2.7, 2.7), heights = c(2.7, 2.7))





# NDMS 18S for only bioelements-----
NMDS_ord_bray_18S = metaMDS(P18S_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray_18S)

stressplot(NMDS_ord_bray_18S)

#build a data frame with NMDS coordinates and metadata
NMDS118S = NMDS_ord_bray_18S$points[,1]
NMDS218S = NMDS_ord_bray_18S$points[,2]
NMDS318S = NMDS_ord_bray_18S$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species
AD_filter_norm_meta_dat
P18S_filter_norm_meta_dat

NMDS_18S = data.frame(NMDS1 = NMDS118S, NMDS2 = NMDS218S, NMDS2 = NMDS318S, Time = P18S_filter_norm_meta_dat$day , Stages = P18S_filter_norm_meta_dat$Stages)


p_time_NMDS12_18S <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, fill = Time)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(subtitle = "18S") +
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) +
  scale_fill_gradient(low = "#9111ab", high = "#e3b039") +
  theme(#axis.line = element_line(color='black'),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 

p_bryo_NMDS12_18s <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, fill = Stages)) +
  #stat_ellipse() + 
  theme_bw(base_size = 8) +
  labs(title = "")  +
  geom_point(shape = 21,size = 2,colour = "black") +
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) + 
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw(base_size = 8) +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 

# 
# p_chla_NMDS12_AD <- ggplot(NMDS_AD, aes(x=NMDS1, y=NMDS2, fill = Chla)) +
#   #stat_ellipse() +
#   theme_bw(base_size = 8) +
#   labs(title = "")  + 
#   geom_point(shape = 21,size = 2,colour = "black") +
#   geom_point(shape = 21, alpha = 0.25, size = 2)
# 
# library(gridExtra)
# grid.arrange(p_subject_NMDS12, p_time_NMDS12, p_subject_NMDS23, p_time_NMDS23, ncol=2, nrow = 2, 
#              layout_matrix = rbind(c(1,2), c(3,4)),
#              widths = c(2.7, 2.7), heights = c(2.7, 2.7))

legend_sys_time <- cowplot::get_legend(p_time_NMDS12_AD + theme(legend.position="right") )
legend_sys_bryo <- cowplot::get_legend(p_bryo_NMDS12_AD + theme(legend.position="right") )



tiff("NMDS.tiff", units="in", width=7.08661, height=3.93701, res=300)


grid.arrange(p_time_NMDS12_16s, p_time_NMDS12_18S, p_time_NMDS12_AD, legend_sys_time, 
             p_bryo_NMDS12_16s, p_bryo_NMDS12_18s, p_bryo_NMDS12_AD,  legend_sys_bryo, 
             nrow=2, ncol = 4, top = textGrob("ASV Composition",gp=gpar(fontsize=12)),
             layout_matrix = rbind(c(1,2,3,4), c(5,6,7,8)),
             widths = c(2.85, 2.7, 2.7, 1.3), heights = c(2.4, 2.5))

dev.off()


#Class: Bacterial community composition ----- 
# agglomerate at Class level
ps_16S_succession_filter_Class <- aggregate_taxa(ps_16S_succession_filter, "class")
# Transform to rel. abundance
ps_16S_succession_filter_Class_norm <- transform_sample_counts(ps_16S_succession_filter_Class, function(x) 100 * x/sum(x))
# Melt to long format
ps_16S_succession_filter_Class_norm_melt <- psmelt(ps_16S_succession_filter_Class_norm)
#Transform to percentage
ps_16S_succession_filter_Class_norm_melt <- aggregate(Abundance ~ OTU + day + Sample, ps_16S_succession_filter_Class_norm_melt, sum)
ps_16S_succession_filter_Class_norm_melt_sum <- aggregate(Abundance ~ day + Sample, ps_16S_succession_filter_Class_norm_melt, sum)
ps_16S_succession_filter_Class_norm_melt_sum <- left_join(ps_16S_succession_filter_Class_norm_melt, ps_16S_succession_filter_Class_norm_melt_sum, by=c("day", "Sample"))
ps_16S_succession_filter_Class_norm_melt_sum_pct <- ps_16S_succession_filter_Class_norm_melt_sum %>% 
  mutate(Abundance_percentage = Abundance.x / Abundance.y *100)
ps_16S_succession_filter_Class_norm_melt_sum_pct$Abundance_percentage <- round(ps_16S_succession_filter_Class_norm_melt_sum_pct$Abundance_percentage, 2)
#Find the 11 most abundant classes for the whole dataset including Phaeobacter (Only 12 colors for plotting)
Top12_Class <- ps_16S_succession_filter_Class_norm_melt_sum_pct %>% group_by(OTU) %>% filter(day < 44) %>% 
  dplyr::summarise('Abundance_percentage' = mean(Abundance_percentage)) %>% arrange(desc(Abundance_percentage)) %>% head(11)

Top12_Class_uniques <- as.factor(Top12_Class[[1]])

#Change Class names to "Others" if not among the 12 most abundant 
ps_16S_succession_filter_Class_norm_melt_sum_pct_1 <- mutate(ps_16S_succession_filter_Class_norm_melt_sum_pct, class = ifelse(!OTU %in% Top12_Class_uniques, "Others", OTU))

# #Clean up names
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- ifelse(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class=="Bacteria_Unclassified_Unclassified", "Unclassified",  phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class)
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment <- ifelse(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment=="Seawater at Day 0", "Seawater",  phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment)
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Environment_1 <- ifelse(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Environment=="Seawater at Day 0", "Planktonic suspension", phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Environment)
# #Reorder and clean up table
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Day <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Day, levels=c("0","1","4","10"))
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class , levels=c("Actinobacteria","Flavobacteriia","Planctomycetia","Sphingobacteriia","Alphaproteobacteria", "Betaproteobacteria", "dTDAproteobacteria","Epsilonproteobacteria","Gammaproteobacteria", "Others","Unclassified", "P. inhibens OTU"))
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment , levels=c("Seawater", "Control", "WT", "dTDA"))
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 <- arrange(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, Day, Treatment)

#Arrange Sample order according to Days 
#Sample_ID_order <- distinct(ps_16S_succession_filter_Class_norm_melt_sum_pct_1, Sample, .keep_all = TRUE)
#Sample_ID_order <- unique(arrange(Sample_ID_order, day)[2])

#ps_16S_succession_filter_Class_norm_melt_sum_pct_1$Sample <- factor(ps_16S_succession_filter_Class_norm_melt_sum_pct_1$Sample, levels=Sample_ID_order)
#ps_16S_succession_filter_Class_norm_melt_sum_pct_1$day <- as.numeric(ps_16S_succession_filter_Class_norm_melt_sum_pct_1$day)

library(ggh4x)

ggplot(ps_16S_succession_filter_Class_norm_melt_sum_pct_1, aes(x=Sample, y = Abundance_percentage, fill = class)) +
  geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
  scale_fill_brewer(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  facet_nested(. ~ day , scales = "free_x") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    #legend.position="none",
    text = element_text(size=10))#,









#Order: Bacterial community composition ----- 
# agglomerate at order level
ps_16S_succession_filter_order <- aggregate_taxa(ps_16S_succession_filter, "order")
# Transform to rel. abundance
ps_16S_succession_filter_order_norm <- transform_sample_counts(ps_16S_succession_filter_order, function(x) 100 * x/sum(x))
# Melt to long format
ps_16S_succession_filter_order_norm_melt <- psmelt(ps_16S_succession_filter_order_norm)
#Transform to percentage
ps_16S_succession_filter_order_norm_melt <- aggregate(Abundance ~ OTU + day + Sample, ps_16S_succession_filter_order_norm_melt, sum)
ps_16S_succession_filter_order_norm_melt_sum <- aggregate(Abundance ~ day + Sample, ps_16S_succession_filter_order_norm_melt, sum)
ps_16S_succession_filter_order_norm_melt_sum <- left_join(ps_16S_succession_filter_order_norm_melt, ps_16S_succession_filter_order_norm_melt_sum, by=c("day", "Sample"))
ps_16S_succession_filter_order_norm_melt_sum_pct <- ps_16S_succession_filter_order_norm_melt_sum %>% 
  mutate(Abundance_percentage = Abundance.x / Abundance.y *100)
ps_16S_succession_filter_order_norm_melt_sum_pct$Abundance_percentage <- round(ps_16S_succession_filter_order_norm_melt_sum_pct$Abundance_percentage, 2)
#Find the 11 most abundant orderes for the whole dataset including Phaeobacter (Only 12 colors for plotting)
Top12_order <- ps_16S_succession_filter_order_norm_melt_sum_pct %>% group_by(OTU) %>% filter(day < 44) %>% 
  dplyr::summarise('Abundance_percentage' = mean(Abundance_percentage)) %>% arrange(desc(Abundance_percentage)) %>% head(11)

Top12_order_uniques <- as.factor(Top12_order[[1]])

#Change order names to "Others" if not among the 12 most abundant 
ps_16S_succession_filter_order_norm_melt_sum_pct_1 <- mutate(ps_16S_succession_filter_order_norm_melt_sum_pct, order = ifelse(!OTU %in% Top12_order_uniques, "Others", OTU))

# #Clean up names
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$order <- ifelse(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$order=="Bacteria_Unorderified_Unorderified", "Unorderified",  phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$order)
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Treatment <- ifelse(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Treatment=="Seawater at Day 0", "Seawater",  phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Treatment)
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Environment_1 <- ifelse(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Environment=="Seawater at Day 0", "Planktonic suspension", phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Environment)
# #Reorder and clean up table
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Day <- factor(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Day, levels=c("0","1","4","10"))
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$order <- factor(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$order , levels=c("Actinobacteria","Flavobacteriia","Planctomycetia","Sphingobacteriia","Alphaproteobacteria", "Betaproteobacteria", "dTDAproteobacteria","Epsilonproteobacteria","Gammaproteobacteria", "Others","Unorderified", "P. inhibens OTU"))
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Treatment <- factor(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1$Treatment , levels=c("Seawater", "Control", "WT", "dTDA"))
# phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1 <- arrange(phylo_main_no_cont_NoOutliers_order_norm_melt_sum_pct_1, Day, Treatment)

#Arrange Sample order according to Days 
#Sample_ID_order <- distinct(ps_16S_succession_filter_order_norm_melt_sum_pct_1, Sample, .keep_all = TRUE)
#Sample_ID_order <- unique(arrange(Sample_ID_order, day)[2])

#ps_16S_succession_filter_order_norm_melt_sum_pct_1$Sample <- factor(ps_16S_succession_filter_order_norm_melt_sum_pct_1$Sample, levels=Sample_ID_order)
#ps_16S_succession_filter_order_norm_melt_sum_pct_1$day <- as.numeric(ps_16S_succession_filter_order_norm_melt_sum_pct_1$day)

library(ggh4x)

ggplot(ps_16S_succession_filter_order_norm_melt_sum_pct_1, aes(x=Sample, y = Abundance_percentage, fill = order)) +
  geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
  scale_fill_brewer(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  facet_nested(. ~ day , scales = "free_x") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    #legend.position="none",
    text = element_text(size=10))#,







#genus: Bacterial community composition ----- 
# agglomerate at genus level
#ps_16S_succession_filter_genus <- aggregate_taxa(ps_16S_succession_filter, "genus")
# Transform to rel. abundance
#ps_16S_succession_filter_genus_norm <- transform_sample_counts(ps_16S_succession_filter_genus, function(x) 100 * x/sum(x))
# Melt to long format
ps_16S_succession_filter_norm_Phaeobacter_melt <- psmelt(ps_16S_succession_filter_norm_Phaeobacter)
#Transform to percentage
ps_16S_succession_filter_norm_Phaeobacter_melt <- aggregate(Abundance ~ OTU + day + Sample, ps_16S_succession_filter_norm_Phaeobacter_melt, sum)
ps_16S_succession_filter_norm_Phaeobacter_melt_sum <- aggregate(Abundance ~ day + Sample, ps_16S_succession_filter_norm_Phaeobacter_melt, sum)
#ps_16S_succession_filter_genus_norm_melt_sum <- left_join(ps_16S_succession_filter_genus_norm_melt, ps_16S_succession_filter_genus_norm_melt_sum, by=c("day", "Sample"))
#ps_16S_succession_filter_genus_norm_melt_sum_pct <- ps_16S_succession_filter_genus_norm_melt_sum %>% 
#  mutate(Abundance_percentage = Abundance.x / Abundance.y *100)
#ps_16S_succession_filter_genus_norm_melt_sum_pct$Abundance_percentage <- round(ps_16S_succession_filter_genus_norm_melt_sum_pct$Abundance_percentage, 2)
#Find the 11 most abundant genuses for the whole dataset including Phaeobacter (Only 12 colors for plotting)
Top12_genus <- ps_16S_succession_filter_genus_norm_melt_sum_pct %>% group_by(OTU) %>% 
  dplyr::summarise('Abundance_percentage' = mean(Abundance_percentage)) %>% arrange(desc(Abundance_percentage)) %>% head(11)

Top12_genus_uniques <- as.factor(Top12_genus[[1]])

#Change genus names to "Others" if not among the 12 most abundant 
#ps_16S_succession_filter_genus_norm_melt_sum_pct_1 <- mutate(ps_16S_succession_filter_genus_norm_melt_sum_pct, genus = ifelse(!OTU %in% Top12_genus_uniques, "Others", OTU))

# #Clean up names
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$genus <- ifelse(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$genus=="Bacteria_Ungenusified_Ungenusified", "Ungenusified",  phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$genus)
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Treatment <- ifelse(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Treatment=="Seawater at Day 0", "Seawater",  phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Treatment)
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Environment_1 <- ifelse(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Environment=="Seawater at Day 0", "Planktonic suspension", phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Environment)
# #Reorder and clean up table
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Day <- factor(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Day, levels=c("0","1","4","10"))
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$genus <- factor(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$genus , levels=c("Actinobacteria","Flavobacteriia","Planctomycetia","Sphingobacteriia","Alphaproteobacteria", "Betaproteobacteria", "dTDAproteobacteria","Epsilonproteobacteria","Gammaproteobacteria", "Others","Ungenusified", "P. inhibens OTU"))
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Treatment <- factor(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1$Treatment , levels=c("Seawater", "Control", "WT", "dTDA"))
# phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1 <- arrange(phylo_main_no_cont_NoOutliers_genus_norm_melt_sum_pct_1, Day, Treatment)

#Arrange Sample order according to Days 
#Sample_ID_order <- distinct(ps_16S_succession_filter_genus_norm_melt_sum_pct_1, Sample, .keep_all = TRUE)
#Sample_ID_order <- unique(arrange(Sample_ID_order, day)[2])

#ps_16S_succession_filter_genus_norm_melt_sum_pct_1$Sample <- factor(ps_16S_succession_filter_genus_norm_melt_sum_pct_1$Sample, levels=Sample_ID_order)
#ps_16S_succession_filter_genus_norm_melt_sum_pct_1$day <- as.numeric(ps_16S_succession_filter_genus_norm_melt_sum_pct_1$day)

library(ggh4x)

ps_16S_succession_filter_norm_Phaeobacter_melt_sum %>%
ggplot(aes(x=Sample, y = Abundance/100000*100)) +
  geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
  scale_fill_brewer(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  facet_nested(. ~ day , scales = "free_x") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    #legend.position="none",
    text = element_text(size=10))#,









###Correlations between copepoda and Phaeobacter vs. Bryozoans ######
tax_table(ps_16S_succession_filter)
ps_16S_succession_filter_Phaeobacter <- subset_taxa(ps_16S_succession_filter, genus == "Phaeobacter")

ps_18S_succession_filter_tax_dat <- data.frame(tax_table(ps_18S_succession_filter))
ps_18S_succession_filter_tax_dat<-ps_18S_succession_filter_tax_dat %>% mutate(ASV = rownames(ps_18S_succession_filter_tax_dat))
tax_table(ps_18S_succession_filter) <- as.matrix(ps_18S_succession_filter_tax_dat)

ps_18S_succession_filter_Copepoda <- subset_taxa(ps_18S_succession_filter, ASV == "ASV_1" | ASV == "ASV_2")
ps_18S_succession_filter_Bryozoans <- subset_taxa(ps_18S_succession_filter, ASV == "ASV_3" | ASV == "ASV_4")

ps_16S_succession_filter_Phaeobacter_dat <- t(as.data.frame(otu_table(ps_16S_succession_filter_Phaeobacter)))

plot((ps_16S_succession_filter_Phaeobacter_dat))
ps_16S_succession_filter_Phaeobacter_dat_sum <- rowSums(ps_16S_succession_filter_Phaeobacter_dat)
ps_16S_succession_filter_Phaeobacter_dat_sum <- data.frame(cbind(rownames(ps_16S_succession_filter_Phaeobacter_dat_sum), ps_16S_succession_filter_Phaeobacter_dat_sum))
colnames(ps_16S_succession_filter_Phaeobacter_dat_sum)[1] <- "Sample"
rownames(ps_16S_succession_filter_Phaeobacter_dat_sum) <- NULL


ps_18S_succession_filter_Copepoda_dat <- t(as.data.frame(otu_table(ps_18S_succession_filter_Copepoda)))

plot((ps_18S_succession_filter_Copepoda_dat))
ps_18S_succession_filter_Copepoda_dat_sum <- rowSums(ps_18S_succession_filter_Copepoda_dat)
ps_18S_succession_filter_Copepoda_dat_sum <- data.frame(cbind(rownames(ps_18S_succession_filter_Copepoda_dat_sum), ps_18S_succession_filter_Copepoda_dat_sum))
ps_18S_succession_filter_Copepoda_dat_sum$Sample <- rownames(ps_18S_succession_filter_Copepoda_dat_sum)
rownames(ps_18S_succession_filter_Copepoda_dat_sum) <- NULL

Phaeobacter_vs_Copepoda <- merge(ps_16S_succession_filter_Phaeobacter_dat_sum,ps_18S_succession_filter_Copepoda_dat_sum, by.x="Sample", by.y="Sample", all.x = TRUE)
colnames(Phaeobacter_vs_Copepoda) <- c("Sample", "Phaeobacter", "Copepoda")


plot((Phaeobacter_vs_Copepoda[,2:3]))





###NORM: Correlations between copepoda and Phaeobacter vs. Bryozoans ######
tax_table(ps_16S_succession_filter_norm)
ps_16S_succession_filter_norm_Phaeobacter <- subset_taxa(ps_16S_succession_filter_norm, genus == "Phaeobacter")

ps_18S_succession_filter_norm_tax_dat <- data.frame(tax_table(ps_18S_succession_filter_norm))
ps_18S_succession_filter_norm_tax_dat<-ps_18S_succession_filter_tax_dat %>% mutate(ASV = rownames(ps_18S_succession_filter_norm_tax_dat))
tax_table(ps_18S_succession_filter_norm) <- as.matrix(ps_18S_succession_filter_norm_tax_dat)

ps_18S_succession_filter_norm_Copepoda <- subset_taxa(ps_18S_succession_filter_norm, ASV == "ASV_1" | ASV == "ASV_2")
ps_18S_succession_filter_norm_Bryozoans <- subset_taxa(ps_18S_succession_filter_norm, ASV == "ASV_3" | ASV == "ASV_4")

ps_16S_succession_filter_norm_Phaeobacter_dat <- t(as.data.frame(otu_table(ps_16S_succession_filter_norm_Phaeobacter)))

plot((ps_16S_succession_filter_norm_Phaeobacter_dat))
ps_16S_succession_filter_norm_Phaeobacter_dat_sum <- rowSums(ps_16S_succession_filter_norm_Phaeobacter_dat)
ps_16S_succession_filter_norm_Phaeobacter_dat_sum <- data.frame(cbind(rownames(ps_16S_succession_filter_norm_Phaeobacter_dat_sum), ps_16S_succession_filter_norm_Phaeobacter_dat_sum))
ps_16S_succession_filter_norm_Phaeobacter_dat_sum$Sample <- rownames(ps_16S_succession_filter_norm_Phaeobacter_dat_sum)
rownames(ps_16S_succession_filter_norm_Phaeobacter_dat_sum) <- NULL


ps_18S_succession_filter_norm_Copepoda_dat <- t(as.data.frame(otu_table(ps_18S_succession_filter_norm_Copepoda)))

plot((ps_18S_succession_filter_norm_Copepoda_dat))
ps_18S_succession_filter_norm_Copepoda_dat_sum <- rowSums(ps_18S_succession_filter_norm_Copepoda_dat)
ps_18S_succession_filter_norm_Copepoda_dat_sum <- data.frame(cbind(rownames(ps_18S_succession_filter_norm_Copepoda_dat_sum), ps_18S_succession_filter_norm_Copepoda_dat_sum))
ps_18S_succession_filter_norm_Copepoda_dat_sum$Sample <- rownames(ps_18S_succession_filter_norm_Copepoda_dat_sum)
rownames(ps_18S_succession_filter_norm_Copepoda_dat_sum) <- NULL

ps_18S_succession_filter_norm_Bryozoans_dat <- t(as.data.frame(otu_table(ps_18S_succession_filter_norm_Bryozoans)))

plot((ps_18S_succession_filter_norm_Bryozoans_dat))
ps_18S_succession_filter_norm_Bryozoans_dat_sum <- rowSums(ps_18S_succession_filter_norm_Bryozoans_dat)
ps_18S_succession_filter_norm_Bryozoans_dat_sum <- data.frame(cbind(rownames(ps_18S_succession_filter_norm_Bryozoans_dat_sum), ps_18S_succession_filter_norm_Bryozoans_dat_sum))
ps_18S_succession_filter_norm_Bryozoans_dat_sum$Sample <- rownames(ps_18S_succession_filter_norm_Bryozoans_dat_sum)
rownames(ps_18S_succession_filter_norm_Bryozoans_dat_sum) <- NULL

Phaeobacter_vs_Copepoda <- merge(ps_16S_succession_filter_norm_Phaeobacter_dat_sum,ps_18S_succession_filter_norm_Copepoda_dat_sum, by.x="Sample", by.y="Sample", all.x = TRUE)
Phaeobacter_vs_Copepoda_vs_Bryozoans <- merge(Phaeobacter_vs_Copepoda,ps_18S_succession_filter_norm_Bryozoans_dat_sum, by.x="Sample", by.y="Sample", all.x = TRUE)
colnames(Phaeobacter_vs_Copepoda_vs_Bryozoans) <- c("Sample", "Phaeobacter", "Copepoda", "Bryozoans")


plot((Phaeobacter_vs_Copepoda_vs_Bryozoans[,2:3]))

plot((Phaeobacter_vs_Copepoda_vs_Bryozoans[,c(2,4)]))



Phaeobacter_vs_Copepoda_vs_Bryozoans %>%   
  ggplot(aes(y=Phaeobacter, x=Bryozoans)) +
  geom_point(alpha=0.6, size = 1) +
  geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(y = "\nPhaeobacter ", x = "\nBryozoans ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")

Phaeobacter_vs_Copepoda_vs_Bryozoans %>%   
  ggplot(aes(y=Phaeobacter, x=Bryozoans)) +
  geom_point(alpha=0.6, size = 1) +
  geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01, size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(y = "\nPhaeobacter ", x = "\nBryozoans ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")


Phaeobacter_vs_Copepoda_vs_Bryozoans %>%   
  ggplot(aes(y=Phaeobacter, x=Copepoda)) +
  geom_point(alpha=0.6, size = 1) +
  geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(y = "\nPhaeobacter ", x = "\nCopepoda ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")

Phaeobacter_vs_Copepoda_vs_Bryozoans %>%   
  ggplot(aes(y=Phaeobacter, x=Copepoda)) +
  geom_point(alpha=0.6, size = 1) +
  geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01, size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(y = "\nPhaeobacter ", x = "\nCopepoda ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")


Phaeobacter_vs_Copepoda_vs_Bryozoans %>%    ###Bryozoans and copepoda is not correlated
  ggplot(aes(y=Bryozoans, x=Copepoda)) +
  geom_point(alpha=0.6, size = 1) +
  geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(y = "\nBryozoans ", x = "\nCopepoda ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")

Phaeobacter_vs_Copepoda_vs_Bryozoans %>%   
  ggplot(aes(y=Bryozoans, x=Copepoda)) +
  geom_point(alpha=0.6, size = 1) +
  geom_smooth(method = lm, se=FALSE)+ 
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01, size = 2)+
  #scale_size(range = c(0.1, 10)) +
  labs(y = "\nBryozoans ", x = "\nCopepoda ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")

