#======================================================================
#MINI PROJECT - GROUP 1
  
#Diff Expression analysis on Plant Transcriptional Response to Viral Infection

#The Biological Question
#How does a crop plant (e.g., Tomato or Arabidopsis or any other plant of interest) alter its gene expression when infected by a plant virus (such as Tobacco Mosaic Virus or Cucumber Mosaic Virus), and which host defense pathways are triggered?
  
#The Aim : 
#Assess the gene expression change in tomato plants infected by Tomato brown rugose fruit virus (ToBRFV) using transcriptome analysis comparing virus inoculated plant and healthy plant and Identify which host defense and metabolic pathways are involved.

#Objectives : 
#1.To identify differentially expressed genes between ToBRFV-infected and mock plants, defining upregulated genes as log2FC > 1 and downregulated genes as log2FC < -1, with FDR (padj) < 0.05.
#2.To classify the up- and downregulated genes into host defense and metabolic pathways by merging the expression data with a gene-to-pathway annotation file using dplyr::inner_join().
#3.To visually the result with a volcano plot and bar chart  
#4.To analyze the result, compile the code, data and literature findings and write it in the report.
#======================================================================
###working directory 
#getwd()

###Environment Setup
#install.packages("readxl")   
library(tidyverse)
library(readxl)

###read the file 
#deg_raw <- read_excel("Table S1.xlsx")
#dim(deg_raw)
#names(deg_raw)

###Data Ingestion, Cleaning & CSV Integration
deg_clean <- deg_raw %>%
select(gene_id,log2FoldChange,gene_description,ToBRFV1,ToBRFV2,ToBRFV3,Mock1,Mock2,Mock3,pvalue,padj,gene_name)

###Log2Fold Calculation
deg_calc <- deg_clean %>%
  mutate(mean_virus = (ToBRFV1 + ToBRFV2 + ToBRFV3) / 3,
         mean_mock  = (Mock1 + Mock2 + Mock3) / 3,
         my_log2FC  = log2((mean_virus + 1) / (mean_mock + 1)))

deg_calc %>% select(gene_id, log2FoldChange, my_log2FC) %>% head(10) #comparison check for the author Log2Fold and our Log2Fold calculation

deg_calc %>%
  summarise(agree = sum(sign(log2FoldChange) == sign(my_log2FC)),
            total = n()) #comparison check for the author Log2Fold and our Log2Fold agree towards same direction 

deg_calc <- deg_calc %>%
  mutate(status = case_when(
    log2FoldChange > 1  & padj < 0.05 ~ "Up",
    log2FoldChange < -1 & padj < 0.05 ~ "Down",
    TRUE ~ "Not significant"))

table(deg_calc$status)

###Generate Volcano Plot
ggplot(deg_calc, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "ToBRFV-infected vs mock tomato leaves",
       x = "log2 fold change",
       y = "-log10(adjusted p-value)") +
  scale_color_manual(values = c("Down" = "pink", "Up" = "steelblue"))
ggsave("volcano_plot.png", width = 7, height = 5, dpi = 300)

deg_calc %>% filter(str_detect(gene_id, "novel")) %>% select(gene_id, log2FoldChange, gene_description)


deg_calc %>%
  filter(str_detect(gene_description, "phenylalanine")) %>%
  select(gene_id, gene_name, log2FoldChange)

