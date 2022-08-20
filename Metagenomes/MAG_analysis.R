setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/Metagenomes/function_annotations/")

library("dplyr")
library("phyloseq")
library("ADImpute")
library(rRDPData)
library(debar)
library(tidyverse)
library("devtools")
library(Rcpp)
library("dada2")
library(plyr)
library("Biostrings")
library("phyloseq")
library(ggplot2)
library( "DESeq2" )
library(vegan)
library(venn)


#Coverage <- read.table("JH21_Co-assembly-GENE-COVERAGES.txt", sep = "\t", header = TRUE)
#MAGs_mean_coverage <- read.table("MAGs_mean_coverage.txt", sep = "\t", header = TRUE)
MAGs_mean_coverage <- read.table("../../CONCOCT_Binned_DBs/SUMMARY/bins_across_samples/mean_coverage_Q2Q3.txt", sep = "\t", header = TRUE)


#Normalize with total reads depth per sample



#rownames(Coverage) <- Coverage$key
rownames(MAGs_mean_coverage) <- MAGs_mean_coverage$bins
#Coverage <- Coverage[,-1]
MAGs_mean_coverage <- MAGs_mean_coverage[,-1]
#colnames(Coverage) <- c("T2","T3","T4","T5","T6","T7","T9","T10","T11","T12")
colnames(MAGs_mean_coverage) <- c("T2","T3","T4","T5","T6","T7","T9","T10","T11","T12")

#d_genes <- (t(round(Coverage)))
d_MAGs <- (t(round(MAGs_mean_coverage)))

d_r <- d[,1:300000]
dim(d)

#rarecurve(d_MAGs, step=10, cex=0.6)
#rarecurve(d_genes, step=10, cex=0.6)


#### MAG analysis Taxonomy ######

Tax <- read.table("bins_summary.txt", sep = "\t", header = TRUE)
BGCs <- read.table("hmms_BGC_domains.txt", sep = "\t", header = TRUE)
Tax <- merge(Tax,BGCs, by.x="bins", by.y="bins", all.x = TRUE)

rownames(Tax) <- Tax$bins
Tax <- Tax[,-1]
Tax$genome_id <- rownames(Tax)

#Load total summary table of BGC counts and add to tax table
df_antismash_6.0.1_summary <- read.csv(file = 'tables/df_antismash_6.0.1_summary.csv', sep = ",", header = TRUE)
df_antismash_6.0.1_summary <- as.data.frame(df_antismash_6.0.1_summary)
rownames(df_antismash_6.0.1_summary) <- df_antismash_6.0.1_summary$genome_id
Tax <- left_join(Tax, df_antismash_6.0.1_summary, "genome_id")
rownames(Tax) <- Tax$genome_id

#CAZymes profile
all_MAGs_CAZymes_raw <- read.csv(file = 'hmmscan_out/output_parser/all_MAGs_CAZymes.csv', sep = "\t", header = FALSE)
all_MAGs_CAZymes_raw <- as.data.frame(all_MAGs_CAZymes_raw)[,c(1,11)]
str(all_MAGs_CAZymes_raw)
colnames(all_MAGs_CAZymes_raw) <- c("CAZyme_Type", "genome_id")

CAZyme_Type <- all_MAGs_CAZymes_raw$CAZyme_Type
CAZyme_Type <- sapply(strsplit(CAZyme_Type, ".hmm"), `[`, 1)
CAZyme_family <- sapply(strsplit(CAZyme_Type, split = "[0-9]"), `[`, 1)
all_MAGs_CAZymes_raw <- cbind(all_MAGs_CAZymes_raw, CAZyme_Type, CAZyme_family)[,-1]

all_MAGs_CAZymes_raw_sum <- all_MAGs_CAZymes_raw %>%
  group_by(genome_id, CAZyme_family) %>%
  dplyr::summarise(CAZyme_family_n = n())


all_MAGs_CAZymes_raw_sum_long <- reshape(data=as.data.frame(all_MAGs_CAZymes_raw_sum),idvar="genome_id",
                                         v.names = "CAZyme_family_n",
                                         timevar = "CAZyme_family",
                                         direction="wide")
colnames(all_MAGs_CAZymes_raw_sum_long) <- c("genome_id", "AA", "CBM", "CE","GH","GT","PL","dockerin","cohesin","SLH")
str(all_MAGs_CAZymes_raw_sum_long)


antiSMASH_contigs_MAGs <- readRDS('antiSMASH_contigs_MAGs.rds')
colnames(antiSMASH_contigs_MAGs[,c(1:5, 8)])
antiSMASH_contigs_MAGs <- antiSMASH_contigs_MAGs[,c(1:5, 8)]


### MAG TMP normalisation ####

size <- data.frame(Tax = rownames(Tax), Length = Tax$total_length)
colnames(size) <- c("hgnc_symbol", "transcript_length") # change header, since the TMP.functions reads only this......
df.tmp <- NormalizeTPM(MAGs_mean_coverage, tr_length = size, scale = 1e+06) # the scale is just a constant, which can be changed, i added a million, so the total sum for each sample i 1.
colSums(df.tmp) # check Sample sum



#Metadata
chla.sum <- readRDS('chla.sum.rds')
chla.sum_suc <- chla.sum %>% filter(Experiment=="Succession")
chla.sum_suc <- chla.sum_suc[,-c(2,4)] 
colnames(chla.sum_suc)[2]<- "mean.chla_suc"
chla.sum_ass <- chla.sum %>% filter(Experiment=="Assembly")
chla.sum_ass$Timepoint <- c(5,6,7,8,9,10,11,12)
chla.sum_ass <- chla.sum_ass[,-c(2,4)] 
colnames(chla.sum_ass)[2]<- "chla.sum_ass"


Meta_data <- as.data.frame(cbind(colnames(df.tmp),chla.sum_suc$mean.chla_suc[-c(1,8)]))
colnames(Meta_data) <- c("Timepoint", "ChlaC_mean_suc")
Meta_data <- as.data.frame((Meta_data))
rownames(Meta_data) <- Meta_data[,1]

chla.sum_ass <- as.data.frame(cbind(Timepoint=Meta_data$Timepoint, chla.sum_ass=c(NA, NA, NA, chla.sum_ass$chla.sum_ass[-4])))

Meta_data$chla.sum_ass <- chla.sum_ass$chla.sum_ass
colnames(Meta_data) <- c("Timepoint", "ChlaC_mean_suc", "ChlaC_mean_ass")
saveRDS(Meta_data, "Meta_data_chlaC_mean")


str(MAGs_mean_coverage)

#MAG Create Phyloseq Object######
Phylo_MAGs <- phyloseq(otu_table(MAGs_mean_coverage, taxa_are_rows = TRUE), 
                       tax_table(as.matrix(Tax)), sample_data(Meta_data))

Phylo_MAGs_norm <- phyloseq(otu_table(df.tmp, taxa_are_rows = TRUE), 
                            tax_table(as.matrix(Tax)), sample_data(Meta_data))




# Melt to long format
Phylo_MAGs_melt <- psmelt(Phylo_MAGs_norm)
Phylo_MAGs_melt$Sample <- factor(Phylo_MAGs_melt$Sample, levels=c("T2","T3", "T4", "T5", "T6", "T7", "T9", "T10","T11", "T12")) 
Phylo_MAGs_melt$AMP.binding <- as.numeric(Phylo_MAGs_melt$AMP.binding)
Phylo_MAGs_melt$Condensation<- as.numeric(Phylo_MAGs_melt$Condensation)
Phylo_MAGs_melt$PP.binding<- as.numeric(Phylo_MAGs_melt$PP.binding)
Phylo_MAGs_melt$ketoacyl.synt<- as.numeric(Phylo_MAGs_melt$ketoacyl.synt)
Phylo_MAGs_melt$total_length<- as.numeric(Phylo_MAGs_melt$total_length)
Phylo_MAGs_melt$Timepoint <- as.numeric(c(2, 3, 4, 5, 6, 7, 9, 10, 11, 12))
Phylo_MAGs_melt$N50 <- as.numeric(Phylo_MAGs_melt$N50)
Phylo_MAGs_melt$ChlaC_mean_suc <- as.numeric(Phylo_MAGs_melt$ChlaC_mean_suc)
Phylo_MAGs_melt$GC_content<- as.numeric(Phylo_MAGs_melt$GC_content)
Phylo_MAGs_melt$bgcs_count <-  as.numeric(Phylo_MAGs_melt$bgcs_count)
Phylo_MAGs_melt$bgcs_on_contig_edge <-  as.numeric(Phylo_MAGs_melt$bgcs_on_contig_edge)
Phylo_MAGs_melt$Mean_coverage <- Phylo_MAGs_melt$Abundance



#Find the 11 most abundant phyla (Only 12 colors for plotting)
Top12_Phyla <- Phylo_MAGs_melt %>% group_by(t_phylum) %>%
  dplyr::summarise('Mean_coverage' = mean(Mean_coverage)) %>% arrange(desc(Mean_coverage)) %>% head(12)
Top12_Phyla$t_phylum[7] <- "Unclassified Phylum"
Phylo_MAGs_melt$t_phylum <- ifelse(Phylo_MAGs_melt$t_phylum=="", "Unclassified Phylum", Phylo_MAGs_melt$t_phylum)
#Change Phyla names to "Others" if not among the 12 most abundant 
Phylo_MAGs_melt_1 <- mutate(Phylo_MAGs_melt, t_phylum = ifelse(t_phylum != "Deinococcota" & t_phylum != "Myxococcota"& 
                                                                 t_phylum != "Bacteroidota" & t_phylum != "Cyanobacteria"&
                                                                 t_phylum != "Proteobacteria" & t_phylum != "Actinobacteriota"&
                                                                 t_phylum != "Bdellovibrionota" & t_phylum != "Planctomycetota"&
                                                                 t_phylum != "Verrucomicrobiota" & t_phylum != "Gemmatimonadota"& 
                                                                 t_phylum != "Unclassified Phylum", "Others", t_phylum))

Phylo_MAGs_melt_1$t_phylum <- factor(Phylo_MAGs_melt_1$t_phylum, levels=c("Deinococcota","Myxococcota", 
                                                                          "Bacteroidota", "Cyanobacteria",
                                                                          "Proteobacteria", "Actinobacteriota", 
                                                                          "Bdellovibrionota" ,"Planctomycetota",
                                                                          "Verrucomicrobiota","Gemmatimonadota",
                                                                          "Unclassified Phylum","Others"))




ggplot(Phylo_MAGs_melt_1, aes(x=Sample, y = Mean_coverage, fill = Phyla)) +
  geom_bar(stat = "identity", color = "black") + labs(y = "Mean_coverage", x="Day") +
  scale_fill_brewer(palette = "Paired") +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="bottom",
    text = element_text(size=10))#,



#Barplot normalized to 100%
Phylo_MAGs_norm_melt <- psmelt(Phylo_MAGs)
Phylo_MAGs_norm_melt$Sample <- factor(Phylo_MAGs_norm_melt$Sample, levels=c("T2","T3", "T4", "T5", "T6", "T7", "T9", "T10","T11", "T12")) 
Phylo_MAGs_norm_melt$AMP.binding <- as.numeric(Phylo_MAGs_norm_melt$AMP.binding)
Phylo_MAGs_norm_melt$Condensation<- as.numeric(Phylo_MAGs_norm_melt$Condensation)
Phylo_MAGs_norm_melt$PP.binding<- as.numeric(Phylo_MAGs_norm_melt$PP.binding)
Phylo_MAGs_norm_melt$ketoacyl.synt<- as.numeric(Phylo_MAGs_norm_melt$ketoacyl.synt)
Phylo_MAGs_norm_melt$total_length<- as.numeric(Phylo_MAGs_norm_melt$total_length)
Phylo_MAGs_norm_melt$Timepoint <- as.numeric(c(2, 3, 4, 5, 6, 7, 9, 10, 11, 12))
Phylo_MAGs_norm_melt$N50 <- as.numeric(Phylo_MAGs_norm_melt$N50)
Phylo_MAGs_norm_melt$ChlaC_mean_suc <- as.numeric(Phylo_MAGs_norm_melt$ChlaC_mean_suc)
Phylo_MAGs_norm_melt$GC_content<- as.numeric(Phylo_MAGs_norm_melt$GC_content)
Phylo_MAGs_norm_melt$Mean_coverage <- Phylo_MAGs_norm_melt$Abundance
#Find the 11 most abundant phyla (Only 12 colors for plotting)
Top12_Phyla <- Phylo_MAGs_norm_melt %>% group_by(t_phylum) %>%
  dplyr::summarise('Mean_coverage' = mean(Mean_coverage)) %>% arrange(desc(Mean_coverage)) %>% head(12)
Top12_Phyla$t_phylum[7] <- "Unclassified Phylum"
Phylo_MAGs_norm_melt$t_phylum <- ifelse(Phylo_MAGs_norm_melt$t_phylum=="", "Unclassified Phylum", Phylo_MAGs_norm_melt$t_phylum)
#Change Phyla names to "Others" if not among the 12 most abundant 
Phylo_MAGs_norm_melt_1 <- mutate(Phylo_MAGs_norm_melt, t_phylum = ifelse(t_phylum != "Deinococcota" & t_phylum != "Myxococcota"& 
                                                                           t_phylum != "Bacteroidota" & t_phylum != "Cyanobacteria"&
                                                                           t_phylum != "Proteobacteria" & t_phylum != "Actinobacteriota"&
                                                                           t_phylum != "Bdellovibrionota" & t_phylum != "Planctomycetota"&
                                                                           t_phylum != "Verrucomicrobiota" & t_phylum != "Gemmatimonadota"& 
                                                                           t_phylum != "Unclassified Phylum", "Others", t_phylum))

Phylo_MAGs_norm_melt_1$t_phylum <- factor(Phylo_MAGs_norm_melt_1$t_phylum, levels=c("Deinococcota","Myxococcota", 
                                                                                    "Bacteroidota", "Cyanobacteria",
                                                                                    "Proteobacteria", "Actinobacteriota", 
                                                                                    "Bdellovibrionota" ,"Planctomycetota",
                                                                                    "Verrucomicrobiota","Gemmatimonadota",
                                                                                    "Unclassified Phylum","Others"))


Purplerain <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963","#580a82", "#f21395", "#7a202c", "#d1bb8a", "#e06f55", "#b78bd6", "#5674a8")
scales::show_col(Purplerain)

main_col <- c("#240785", "#f21395", "#e0b62b")


Phylo_MAGs_norm_melt_1$Time = 0
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T2"] <- 7
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T3"] <- 10
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T4"] <- 15
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T5"] <- 23
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T6"] <- 29
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T7"] <- 44
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T9"] <- 71
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T10"] <- 85
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T11"] <- 99
Phylo_MAGs_norm_melt_1$Time[Phylo_MAGs_norm_melt_1$Sample=="T12"] <- 113

Phylo_MAGs_norm_melt_1$Time <- factor(Phylo_MAGs_norm_melt_1$Time)

colnames(Phylo_MAGs_norm_melt_1)[14] <- "Phyla"
 

tiff("MAG_ABS_BARPLOT050722.tiff", units="in", width=4, height=4, res=300)

Phylo_MAGs_norm_melt_1 %>% filter(Time != "7") %>% 
ggplot(aes(x=Time, y = Mean_coverage, fill = Phyla)) +
  geom_bar(stat = "identity") + xlab("Day") + ylab("Q2–Q3 mean coverage") +
  scale_fill_manual(values = Purplerain) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x=element_blank(),
   )#,
dev.off()

Phylo_MAGs_norm_melt_1 %>% filter(Phyla=="Myxococcota") %>% filter(Sample == "T6", AMP.binding > 1) %>% 
  ggplot(aes(y=AMP.binding, x=Time, size = Mean_coverage, col = Phyla)) +
  geom_point(alpha=0.5) +  xlab("Day") + ylab("AMP binding domian counts") +
  scale_size(range = c(0.1, 10)) +
  scale_color_manual(values = Purplerain) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    #plot.margin = margin(0, 0, 0, -2, "cm")
    #axis.title.x=element_blank(),
  )#,
AMP.binding_plot <- Phylo_MAGs_norm_melt_1 %>% filter(OTU!="Unbinned_Contigs") %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% 
  ggplot(aes(y=AMP.binding, x=Time, size = Mean_coverage, col = Phyla)) +
  geom_point(alpha=0.5) +  xlab("Day") + ylab("AMP binding domian counts") +
  scale_size(range = c(0.1, 10)) +
  scale_color_manual(values = Purplerain) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    #plot.margin = margin(0, 0, 0, -2, "cm")
    #axis.title.x=element_blank(),
  )#,


Summary_AMPbindings <- Phylo_MAGs_norm_melt_1 %>% filter(OTU!="Unbinned_Contigs") %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% group_by(Time) %>% dplyr::count(AMP.binding)
  
Summary_AMPbindings_plot <- Summary_AMPbindings %>% group_by(Time) %>% dplyr::summarise(count = sum(AMP.binding)) %>%
  ggplot(aes(x=Time, y=count)) +
  geom_bar(stat="identity", fill = "#e0b62b", col = "black") +
  #scale_color_manual(values = Purplerain) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.x = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    #axis.text.y = element_blank()
  )#,
Phylo_MAGs_norm_melt_1 %>% filter(OTU!="Unbinned_Contigs") %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% dplyr::count(AMP.binding) %>% dplyr::summarise(median(AMP.binding))

Phylo_MAGs_norm_melt_1 %>% filter(OTU=="JH21_MAG_00132") 
# 
# unique_MAGs_23_29 <- Phylo_MAGs_norm_melt_1 %>% filter(OTU!="Unbinned_Contigs") %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% filter(Time %in% c(10, 15, 23, 29))
# unique_MAGs_44_71 <- Phylo_MAGs_norm_melt_1 %>% filter(OTU!="Unbinned_Contigs") %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% filter(Time %in% c(44, 71, 85, 99, 113))

# 
# library(venn)
# unique_MAGs_23_29_venn_unique_MAGs_44_71 <- venn(list(unique_MAGs_23_29$OTU, unique_MAGs_44_71$OTU), snames = "23_29, 44_71")
# str(unique_MAGs_23_29_venn_unique_MAGs_44_71)
# unique_MAGs_23_29_venn_unique_MAGs_44_71$`23_29`
# unique_MAGs_23_29 <- attr(unique_MAGs_23_29_venn_unique_MAGs_44_71, "intersections")[[2]]
# 
# 
# 
# Phylo_MAGs_norm_melt_1 %>% filter(OTU %in% unique_MAGs_23_29) %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% 
#   ggplot(aes(y=AMP.binding, x=Time, size = Mean_coverage, col = Phyla)) +
#   geom_point(alpha=0.5) +  xlab("Day") + ylab("AMP binding domian counts") +
#   scale_size(range = c(0.1, 10)) +
#   scale_color_manual(values = Purplerain) +
#   #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
#   #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
#   theme_bw(base_size = 8) +
#   theme(#axis.line = element_line(color='black'),
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.position = "none",
#     #plot.margin = margin(0, 0, 0, -2, "cm")
#     #axis.title.x=element_blank(),
#   )#,


Hist_plot <- Phylo_MAGs_norm_melt_1 %>% filter(OTU!="Unbinned_Contigs") %>% filter(Sample != "T2", AMP.binding > 1, Mean_coverage > 1) %>% 
ggplot( aes(x=AMP.binding)) + 
  geom_histogram(bins = 30, fill = "#e0b62b", col = "black") +
  #scale_color_manual(values = Purplerain) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw(base_size = 8) + coord_flip() + # xlim(0,130) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    #axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0, 0, 0.5, 0, "cm")
  )#,

library(cowplot)

tiff("AMPbinding.tiff", units="in", width=5, height=5, res=300)
plot_grid(Summary_AMPbindings_plot, NULL, AMP.binding_plot, Hist_plot , 
          ncol = 2, align = 'v', rel_heights = c(0.3, 1), rel_widths = c(1, 0.5))
dev.off()

legend_AMP.binding_plot <- cowplot::get_legend(AMP.binding_plot + theme(legend.position="right")  + guides(size=guide_legend(title="Q2–Q3 mean coverage")))

tiff("AMPbinding_legend.tiff", units="in", width=2, height=6, res=300)
plot_grid(legend_AMP.binding_plot)
dev.off()          
