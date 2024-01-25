# Code to analyse the alpha and beta diversity

```R
library(tidyverse)
library(microbiome)
library(phyloseq)
library(knitr)
library(metagMisc)
library(picante)
library(data.table)
library(ggpubr)
library(vegan)
```
### Alpha div

```R
Plotting refractory ASVs
otu_tab <- t(abundances(pseq))
ASVs <- vegan::rarecurve(otu_tab, 
                         step = 50, label = FALSE, 
                         sample = min(rowSums(otu_tab), 
                                      col = "blue", cex = 0.6))
p.rar <- plot_taxa_prevalence(pseq, "Phylum")
#rarefying at 3000 keeps all the samples, so go for it
pseq.rar <- rarefy_even_depth(pseq, rngseed=1, sample.size=3000, replace=F)
#17OTUs were removed because they are no longer present in any sample after random subsampling

#Let's extract relevant index for alpha diversity:
#Diversity - Shannon
#Richness - Chao1
#Dominance/Evenness - simpson
#Rarity and low abundance

div_tab <- microbiome::alpha(pseq.rar, index = "all") %>%
  select('diversity_shannon', 'chao1', 'diversity_gini_simpson', 
         'diversity_inverse_simpson', 'observed') %>%
  rownames_to_column(var = 'sample-id')

#Now, let's calculate faith-PD
pseq.asvtab <- as.data.frame(pseq.rar@otu_table)
pseq.tree <- pseq.rar@phy_tree #check it is rooted

div_pd <- pd(t(pseq.asvtab), pseq.tree, include.root=T) %>%
  rownames_to_column(var = 'sample-id')
#t(ou_table) transposes the table for use in picante
# and the tree file comes from the first code chunck we used to read tree file
#(see making a phyloseq object section).
print(div_pd)

#Now, we modify the metadata for joining
div_meta <- rownames_to_column(metadata, var = 'sample-id')

#And joining
div_analysis <- left_join(div_meta, div_tab, by = 'sample-id') %>%
  left_join(., div_pd, by = 'sample-id')
```

#### Alpha diversity and clinical groups

```R
#first, we will check for normality of the diversity parameters
shapiro.test(div_analysis$diversity_shannon)
#W = 0.96648, p-value = 0.09288 - normal distribution
shapiro.test(div_analysis$chao1)
#W = 0.98103, p-value = 0.4624 - normal distribution
shapiro.test(div_analysis$diversity_gini_simpson)
#W = 0.88777, p-value = 4.245e-05 - not normal distribution
shapiro.test(div_analysis$PD)
#W = 0.98258, p-value = 0.5355 -normal distribution

#Now, let's check if there is any difference between the groups and div par.
#Shannon
an_sha <- aov(diversity_shannon ~ clingroup, data = div_analysis)
summary(an_sha)#p-value = 0.0398
an_sha_2 <- aov(diversity_shannon ~ clingroup_2, data = div_analysis)
summary(an_sha_2)#p-value = 0.019

#Chao1
an_chao <- aov(chao1 ~ clingroup, data = div_analysis)
summary(an_chao)#p-value = 0.0169
an_chao_2 <- aov(chao1 ~ clingroup_2, data = div_analysis)
summary(an_chao_2)#p-value = 0.00503

#Simpson
kruskal.test(diversity_gini_simpson ~ clingroup, data = div_analysis)#p-value = 0.034
kruskal.test(diversity_gini_simpson ~ clingroup_2, data = div_analysis)#p-value = 0.012

#PD
an_PD <- aov(PD ~ clingroup, data = div_analysis)
summary(an_PD)#p-value = 0.0148
an_PD_2 <- aov(PD ~ clingroup_2, data = div_analysis)
summary(an_PD_2)#p-value = 0.00821

#Adjusting p-values

pval_3group <- c(0.039, 0.017, 0.034, 0.015)
p.adjust(pval_3group, method = 'bonferroni')
#Sha-0.156; Chao-0.068; Simp-0.136; PD-0.060

pval_2group <- c(0.019, 0.005, 0.013, 0.008)
p.adjust(pval_2group, method = 'bonferroni')
#Sha-0.076; Chao-0.020; Simp-0.052; PD-0.032

#Now, lets do some plots

#Shannon
box_sha_1 <- ggplot(div_analysis, aes(x = clingroup, y = diversity_shannon, fill = clingroup)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Shannon') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

box_sha_2 <- ggplot(div_analysis, aes(x = clingroup_2, y = diversity_shannon, fill = clingroup_2)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=90, c=200) +
  theme_bw() +
  ylab('Shannon') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Chao1
box_chao_1 <- ggplot(div_analysis, aes(x = clingroup, y = chao1, fill = clingroup)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Chao1') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

box_chao_2 <- ggplot(div_analysis, aes(x = clingroup_2, y = chao1, fill = clingroup_2)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=90, c=200) +
  theme_bw() +
  ylab('Chao1') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Simpson
box_sim_1 <- ggplot(div_analysis, aes(x = clingroup, y = diversity_gini_simpson, fill = clingroup)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Simpson') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

box_sim_2 <- ggplot(div_analysis, aes(x = clingroup_2, y = diversity_gini_simpson, fill = clingroup_2)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=90, c=200) +
  theme_bw() +
  ylab('Simpson') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#PD
box_PD_1 <- ggplot(div_analysis, aes(x = clingroup, y = PD, fill = clingroup)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('PD') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

box_PD_2 <- ggplot(div_analysis, aes(x = clingroup_2, y = PD, fill = clingroup_2)) +
  geom_boxplot(color = 'black') + scale_fill_hue(l=90, c=200) +
  theme_bw() +
  ylab('PD') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Let's arrange
boxplot_div_2gruop <- ggarrange(box_sha_2, box_chao_2, box_sim_2, box_PD_2,
                                ncol = 4)
boxplot_div_3gruop <- ggarrange(box_sha_1,  box_chao_1, box_sim_1, box_PD_1, ncol = 4)
boxplot_all_div <- ggarrange(boxplot_div_2gruop, boxplot_div_3gruop, ncol = 1, labels = c('a', 'b'))#first we check for normality
div_analysis$lcn2 <- as.numeric(div_analysis$lcn2)
shapiro.test(div_analysis$lcn2)
#W = 0.86855, p-value = 9.717e-06, so we have to go for non-parametric tests

#Let's estimate the spearmann cor between var
res_shannon <-cor.test(div_analysis$lcn2, div_analysis$diversity_shannon,  method = "spearman")
res_shannon
res_chao1 <-cor.test(div_analysis$lcn2, div_analysis$chao1,  method = "spearman")
res_chao1
res_diversity_gini_simpson <-cor.test(div_analysis$lcn2, div_analysis$diversity_gini_simpson,  method = "spearman")
res_diversity_gini_simpson
res_PD <-cor.test(div_analysis$lcn2, div_analysis$PD,  method = "spearman")
res_PD

#Now, lets plot diversity estimates and lcn2

#Shannon
print(res_shannon$estimate)
print(res_shannon$p.value)

sha <- ggplot(div_analysis, aes(x = diversity_shannon, y = lcn2)) + 
  geom_point(col = "purple", size = 2) +
  ylab('LCN-2') + 
  xlab("Shannon") +
  annotate("text", x = 4.5, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.08", ", pval = 0.55")) +
  theme_light()

#Chao1
print(res_chao1$estimate)
print(res_chao1$p.value)

cha <- ggplot(div_analysis, aes(x = chao1, y = lcn2)) + 
  geom_point(col = "turquoise4", size = 2) +
  ylab('LCN-2') + 
  xlab("Chao1") +
  annotate("text", x = 140, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.09", ", pval = 0.49")) +
  theme_light()

#Gini Simpson
print(res_diversity_gini_simpson$estimate)
print(res_diversity_gini_simpson$p.value)

simp <- ggplot(div_analysis, aes(x = diversity_inverse_simpson, y = lcn2)) + 
  geom_point(col = "gold2", size = 2) +
  ylab('LCN-2') + 
  xlab("Simpson") +
  annotate("text", x = 80, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.05", ", pval = 0.68")) +
  theme_light()

#PD
print(res_PD$estimate)
print(res_PD$p.value)

pd <- ggplot(div_analysis, aes(x = diversity_shannon, y = lcn2)) + 
  geom_point(col = "orangered", size = 2) +
  ylab('LCN-2') + 
  xlab("Faith PD") +
  annotate("text", x = 4.5, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.19", ", pval = 0.13")) +
  theme_light()

ggarrange(sha, cha, simp, pd, labels = c('a', 'b', 'c', 'd'))
ggsave('Figures2/FigS1_alpha_div_lcn2.tiff', width = 9, height = 7)

#Now, we repeat the analysis but only using TD patients
#first we select ill people
div_analysis_ill <- div_analysis %>% filter(., clingroup_2 == 1)
shapiro.test(div_analysis_ill$lcn2)
#W = 0.86777, p-value = 0.0001085, so we have to go for non-parametric tests

#Let's estimate the spearmann cor between var
res_shannon <-cor.test(div_analysis_ill$lcn2, div_analysis_ill$diversity_shannon,  method = "spearman")
res_shannon
res_chao1 <-cor.test(div_analysis_ill$lcn2, div_analysis_ill$chao1,  method = "spearman")
res_chao1
res_diversity_gini_simpson <-cor.test(div_analysis_ill$lcn2, div_analysis_ill$diversity_gini_simpson,  method = "spearman")
res_diversity_gini_simpson
res_PD <-cor.test(div_analysis_ill$lcn2, div_analysis_ill$PD,  method = "spearman")
res_PD

#Now, lets plot diversity estimates and lcn2
#Shannon
print(res_shannon$estimate)
print(res_shannon$p.value)

sha <- ggplot(div_analysis_ill, aes(x = diversity_shannon, y = lcn2)) + 
  geom_point(col = "purple", size = 2) +
  ylab('LCN-2') + 
  xlab("Shannon") +
  annotate("text", x = 4.5, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.13", ", pval = 0.40")) +
  theme_light()

#Chao1
print(res_chao1$estimate)
print(res_chao1$p.value)

cha <- ggplot(div_analysis_ill, aes(x = chao1, y = lcn2)) + 
  geom_point(col = "turquoise4", size = 2) +
  ylab('LCN-2') + 
  xlab("Chao1") +
  annotate("text", x = 140, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.15", ", pval = 0.32")) +
  theme_light()

#Gini Simpson
print(res_diversity_gini_simpson$estimate)
print(res_diversity_gini_simpson$p.value)

simp <- ggplot(div_analysis_ill, aes(x = diversity_inverse_simpson, y = lcn2)) + 
  geom_point(col = "gold2", size = 2) +
  ylab('LCN-2') + 
  xlab("Simpson") +
  annotate("text", x = 80, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.10", ", pval = 0.51")) +
  theme_light()

#PD
print(res_PD$estimate)
print(res_PD$p.value)

pd <- ggplot(div_analysis_ill, aes(x = diversity_shannon, y = lcn2)) + 
  geom_point(col = "orangered", size = 2) +
  ylab('LCN-2') + 
  xlab("Faith PD") +
  annotate("text", x = 4.5, y = 17500, col = 'grey32',
           label = paste0(" rho = -0.20", ", pval = 0.18")) +
  theme_light()

ggarrange(sha, cha, simp, pd, labels = c('a', 'b', 'c', 'd'))
ggsave('Figures/FigS2_alpha_div_lcn2_ill.tiff', width = 9, height = 7)
```

### Beta diversity

#### Relation with LCN2 and clingroup

```R
#Unweighted Unifrac
ordu.unwt.uni <- ordinate(pseq.rar, "PCoA", "unifrac", weighted=F)
barplot(ordu.unwt.uni$values$Eigenvalues[1:10])

#Plotting to know the variance explained
unwt.unifrac <- plot_ordination(pseq.rar, ordu.unwt.uni, color = 'lcnq3hi') +
  ggtitle("Unweighted UniFrac") + geom_point(size = 2) + 
  theme_classic() + scale_color_brewer()


eigen_vec_unw <- rownames_to_column(data.frame(ordu.unwt.uni$vectors), var = 'sample-id')
beta_div_unw <- left_join(div_meta, eigen_vec_unw, by = 'sample-id')
beta_div_unw$lcn2 <- as.numeric(beta_div_unw$lcn2)

poca_unw_lcn2 <- ggplot(beta_div_unw, aes(x = Axis.1, y = Axis.2, col = lcn2)) + 
  geom_point(size = 3) +
  theme_light() +
  xlab('PC1 (9.4%)') +
  ylab('PC2 (6%)') + 
  scale_color_viridis()  +
  labs(color='LCN2')

poca_unw_group <- ggplot(beta_div_unw, aes(x = Axis.1, y = Axis.2, col = clingroup)) + 
  geom_point(size = 3) +
  theme_light() +
  xlab('PC1 (9.4%)') +
  ylab('PC2 (6%)') +
  labs(color='Group') +
  scale_color_manual(values=c("aquamarine4", "aquamarine1", "orange")) +
  theme(legend.key.width = unit(1.4, 'cm'))

#Weighted Unifrac
ordu.wt.uni <- ordinate(pseq.rar, "PCoA", "unifrac", weighted=T)
barplot(ordu.wt.uni$values$Eigenvalues[1:10])

#Plotting to know the variance explained
wt.unifrac <- plot_ordination(pseq, ordu.wt.uni, color = 'lcnq3hi') +
  ggtitle("Weighted UniFrac") + geom_point(size = 2) + 
  theme_classic() + scale_color_brewer()


eigen_vec_wt <- rownames_to_column(data.frame(ordu.wt.uni$vectors), var = 'sample-id')
beta_div_wt <- left_join(div_meta, eigen_vec_wt, by = 'sample-id')
beta_div_wt$lcn2 <- as.numeric(beta_div_wt$lcn2)

poca_wt_lcn2 <- ggplot(beta_div_wt, aes(x = Axis.1, y = Axis.2, col = lcn2)) + 
  geom_point(size = 3) +
  theme_light() +
  xlab('PC1 (7.4%)') +
  ylab('PC2 (6.6%)') + 
  scale_color_viridis() +
  labs(color='LCN2')

poca_wt_group <- ggplot(beta_div_wt, aes(x = Axis.1, y = Axis.2, col = clingroup)) + 
  geom_point(size = 3) +
  theme_light() +
  xlab('PC1 (7.4%)') +
  ylab('PC2 (6.6%)') +
  labs(color='Group') +
  scale_color_manual(values=c("aquamarine4", "aquamarine1", "orange")) +
  theme(legend.key.width = unit(1.4, 'cm'))

group <- ggarrange(poca_unw_group, poca_wt_group, labels = c('c', 'd'), common.legend = T, legend = 'right')
lcn2 <- ggarrange(poca_unw_lcn2, poca_wt_lcn2, labels = c('e', 'f'), common.legend = T, legend = 'right')
beta_div <- ggarrange(group, lcn2, nrow = 2)
ggarrange(boxplot_all_div, beta_div, nrow = 2)
ggsave('Figures/Fig1_alpha_beta_div.tiff', width = 9, height = 12)
```
