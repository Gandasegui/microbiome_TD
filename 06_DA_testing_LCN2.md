# Code to perfomr DA testing by soearman correlation, SLDEx and ML
#### At phylum, family and genus level

```R
library(tidyverse)
library(microbiome)
library(phyloseq)
library(metagMisc)
library(data.table)
library(ggpubr)
library(vegan)
library('ALDEx2')
library(mia)
library("SHAPforxgboost")
library("here")
library(xgboost)
library(patchwork)
```

### At PHYLUM level

```R
#Let's rename the object, just in case
ps.phy <- pseq
#pseq@sam_data[["cured_all"]] <- str_replace(pseq@sam_data[["cured_all"]], '0', 'not_cured')

#Let's check the characteristics of this new phyloseq object
microbiome::summarize_phyloseq(ps.phy)

#Composition analysis
#We have to covert the lcn2 variable to numeric
ps.phy@sam_data[["lcn2"]] <- as.numeric(ps.phy@sam_data[["lcn2"]])

# Make sure we use functions from correct package
transform <- microbiome::transform

#Now, let's see what happens at phylim level
# Merge rare taxa to speed up examples
pseq.comp <- transform(ps.phy, "compositional")
pseq.comp.phy <- aggregate_rare(pseq.comp, level = "Phylum", detection = 1/100, prevalence = 10/100)

#LCN-2
# Let's plot the composition sorted by LCN-2
comp_phy_bar <- plot_composition(pseq.comp.phy,
                             sample.sort = "lcn2",
                             transform = "compositional") +
  scale_fill_brewer("Phylum", palette = "Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic")) +
  xlab('Samples sorted by LCN-2 values')

#There is some evidence of an association between phylum and lcn-2
#Now, let's generate a new tibble for statistical analysis and representation
comp_phy_df <- t(as.data.frame(otu_table(pseq.comp.phy))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

#Now, let's join with metadata of interest
comp_phy_df <- select(div_analysis, 'sample-id','lcn2', 'clingroup', 'clingroup_2') %>%
  left_join(comp_phy_df, ., by = 'sample-id')

#Conventional sp
#Bacteroidota
res_B <-cor.test(comp_phy_df$lcn2, comp_phy_df$Bacteroidota,  method = "spearman")
res_B#rho=0.21; pval=0.106
#Firmicutes
res_F <-cor.test(comp_phy_df$lcn2, comp_phy_df$Firmicutes,  method = "spearman")
res_F#rho=-0.37; pval=0.003
#Proteobacteria
res_P <-cor.test(comp_phy_df$lcn2, comp_phy_df$Proteobacteria,  method = "spearman")
res_P#rho=0.47; pval=1.8e-4

#Adjusting p-values
pval_phy <- c(0.106, 0.003, 0.00018)
p.adjust(pval_phy, method = 'bonferroni')
#Bac-0.318; Fir-0.009; Pro-5.4e-4

#ALDE-x2
# summexp
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.phy)
#Agglomerate by genus and subset by prevalence
tse_phy <- mia::subsetByPrevalentTaxa(tse, rank = "Phylum", prevalence = 10/100)
#Transform count assay to relative abundances
tse_phy <- mia::transformCounts(tse_phy, method = "relabundance")
detach("package:mia", unload = TRUE)
#Run the test
#lcn2
cont.var <- as.numeric(tse_phy$lcn2)
x <- aldex.clr(assay(tse_phy), tse_phy$lcn2)
res_phy_aldex <- aldex.corr(x, cont.var)
#Bacteoidota
res_phy_aldex[4,(4:6)]#rho=0.10; p-val=0.585
#Firmicutes
res_phy_aldex[8,(4:6)]#rho=-0.08; p-val=0.666
#Proteobacteria
res_phy_aldex[5,(4:6)]#rho=0.38; p-val=0.041

#Scatter plots
#Scatter plot Bacterioidetes
comp_B <- ggplot(comp_phy_df, aes(x = Bacteroidota, y = lcn2)) + 
  geom_point(col = "salmon2", size = 2) +
  ylab('LCN-2') + 
  xlab("Bacteroidota Abundance") +
  theme_light()  +
  annotate("text", x = 0.5, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = 0.20, p-val = 0.318")+
  annotate("text", x = 0.5, y = 17000, col = 'black', size = 3,
           label = "ALDEx: rho = 0.10, p-val = 0.585") +
  xlim(0,0.8)
#Scatter plot Firmicutes
comp_F <- ggplot(comp_phy_df, aes(x = Firmicutes, y = lcn2)) + 
  geom_point(col = "palevioletred", size = 2) +
  ylab('LCN-2') + 
  xlab("Firmicutes Abundance") +
  theme_light()  +
  annotate("text", x = 0.5, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = -0.37, p-val = 0.009")+
  annotate("text", x = 0.5, y = 17000, col = 'black',  size = 3,
           label = "ALDEx: rho = 0.08, p-val = 0.666") +
  xlim(0,0.8)
#Scatter plot Proteobacteria
comp_P <- ggplot(comp_phy_df, aes(x = Proteobacteria, y = lcn2)) + 
  geom_point(col = "gold2", size = 2) +
  ylab('LCN-2') + 
  xlab("Proteobacteria Abundance") +
  theme_light()  +
  annotate("text", x = 0.5, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = 0.47, p-val = 5.4e-4")+
  annotate("text", x = 0.5, y = 17000, col = 'black', size = 3,
           label = "ALDEx: rho = 0.38, p-val = 0.041") +
  xlim(0,0.8)

#Let's merge
bar <- ggarrange(comp_phy_bar)
scat <- ggarrange(comp_B, comp_F, comp_P, nrow = 1)
ggarrange(bar, scat, ncol = 1, labels = c('a', 'b'))

#Let's continue with DAA based on SHAP
xg_phy <- comp_phy_df %>%
  dplyr::select('sample-id', 'lcn2', "clingroup_2", "Actinobacteriota", 
         "Bacteroidota", "Desulfobacterota", "Firmicutes", "Other", 
         "Proteobacteria", "Verrucomicrobiota")
colnames(xg_phy) <- str_replace(colnames(xg_phy), "clingroup_2", "1-No TD; 2 = TD")

phy <- c(colnames(xg_phy))[-(1:2)]


dataX <- data.matrix(xg_phy[phy])
mod <- xgboost::xgboost(data = dataX, 
                        label = xg_phy$lcn2, 
                        params = list(objective = "reg:squarederror", learning_rate = 1),
                        nrounds = 500)
# To return the SHAP values and ranked features by mean|SHAP|
shap_values <- shap.values(xgb_model = mod, X_train = dataX)
# The ranked features by mean |SHAP|
shap_values$mean_shap_score
# To prepare the long-format data:
shap_long <- shap.prep(xgb_model = mod, X_train = dataX)
# is the same as: using given shap_contrib
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)
# **SHAP summary plot**
shap_imp_phy <- shap.plot.summary(shap_long)
# SHAP first 4 var
fig_list <- lapply(names(shap_values$mean_shap_score)[1:4], 
                   shap.plot.dependence, data_long = shap_long)
shap_imp_top_phy <- gridExtra::grid.arrange(grobs = fig_list, ncol = 2)

pro <- shap.plot.dependence(data_long = shap_long, x = 'Proteobacteria', y = 'Proteobacteria', color_feature = "1-No TD; 2 = TD")
act <- shap.plot.dependence(data_long = shap_long, x = 'Actinobacteriota', y = 'Actinobacteriota', color_feature = "1-No TD; 2 = TD")
bac <- shap.plot.dependence(data_long = shap_long, x = 'Bacteroidota', y = 'Bacteroidota', color_feature = "1-No TD; 2 = TD")
des <- shap.plot.dependence(data_long = shap_long, x = 'Desulfobacterota', y = 'Desulfobacterota', color_feature = "1-No TD; 2 = TD")
shap_imp_top_phy <- ggarrange(pro, act, bac, des, common.legend = T)
#arrange
phy_comp <- ggarrange(bar, scat, labels = c('a', 'b') , ncol = 1, vjust = 1,
                      font.label = list(size = 18))
phy_shap <- ggarrange(shap_imp_phy, shap_imp_top_phy, labels = c('c', 'd'),
                      widths = c(1.2, 1), vjust = 1, font.label = list(size = 18))
ggarrange(phy_comp, phy_shap, ncol = 1, heights = c(1.3, 1))
ggsave('Figures/Fig2_DA_phylum.tiff', width = 11, height = 14)
```

### At FAMILY level

```R
#Let's rename the object, just in case
ps.fam <- pseq

#Let's check the characteristics of this new phyloseq object
microbiome::summarize_phyloseq(ps.fam)

#Composition analysis
#We have to covert the lcn2 variable to numeric
ps.fam@sam_data[["lcn2"]] <- as.numeric(ps.fam@sam_data[["lcn2"]])

# Make sure we use functions from the correct package
transform <- microbiome::transform

#Now, lets see what happens at phylim level
# Merge rare taxa to speed up examples
pseq.comp <- transform(ps.fam, "compositional")
pseq.comp.fam <- aggregate_rare(pseq.comp, level = "Family", detection = 1/100, prevalence = 10/100)
#pseq.com.fam has 28 families

#LCN-2
#Lets plot the composition sorted by LCN-2
comp_fam_bar <- plot_composition(pseq.comp.fam,
                                 sample.sort = "lcn2",
                                 transform = "compositional") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic")) +
  xlab('Samples sorted by LCN-2 values')

#There are too many families to see something clear
#Now, let's generate a new tibble for statistical analysis and representation
comp_fam_df <- t(as.data.frame(otu_table(pseq.comp.fam))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

#Now, let's join with metadata of interest
comp_fam_df <- dplyr::select(div_analysis, 'sample-id','lcn2', 'clingroup', 'clingroup_2') %>%
  left_join(comp_fam_df, ., by = 'sample-id')

#Conventional spearman
pseq_all_families <- colnames(comp_fam_df)[-1]
pseq_all_families <- pseq_all_families[-(29:31)]

comp_sp_fam = data_frame(family = character(), rho = double(), pval = double())
for (i in pseq_all_families){
  # Creating new variables makes the code clearer
  x <- comp_fam_df %>% pull(i)
  x <- as.numeric(x)
  y <- comp_fam_df %>% pull(lcn2)
  cor <- cor.test(x, y, method="spearman")
  rho <- cor$estimate
  pval <- cor$p.value
  comp_sp_fam <- comp_sp_fam %>% add_row('family' = i, 'rho' = rho, 'pval' = pval)
}
#Let's adjust by bonferroni
p_adj_bon <- comp_sp_fam %>% pull(pval)
p_adj_bon <- p.adjust(p_adj_bon, method = 'bonferroni')
comp_sp_fam$padjusted <- p_adj_bon

#Let's select those near significant
sig_fam <- comp_sp_fam %>% 
  filter(padjusted < 0.06)
#Three families are associated
sig_fam
#Enterobac - 0.39; 0.059, Oscillospi - -0.42, 0.020; Ruminococ - -0.41, 0.033

#ALDE-x2
# summexp
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.fam)
#Agglomerate by genus and subset by prevalence
tse_fam <- mia::subsetByPrevalentTaxa(tse, rank = "Family", prevalence = 10/100)
#Transform count assay to relative abundances
tse_fam <- mia::transformCounts(tse_fam, method = "relabundance")
detach("package:mia", unload = TRUE)
#Run the test
#lcn2
cont.var <- as.numeric(tse_fam$lcn2)
x <- aldex.clr(assay(tse_fam), tse_fam$lcn2)
res_fam_aldex <- aldex.corr(x, cont.var)
view(res_fam_aldex)
#Enterobac - 0.35; 0.105, Oscillospi - -0.45, 0.017; Ruminococ - -0.25, 0.367

#Let's do a table for supplementary material
aldex_sp_fam_res <- res_fam_aldex %>%
  rownames_to_column(var = 'family') %>%
  dplyr::select('family', 'spearman.erho', 'spearman.ep', 'spearman.eBH') %>%
  as_tibble() %>%
  left_join(comp_sp_fam, ., by = 'family') %>%
  arrange(., pval)

colnames(aldex_sp_fam_res) <- c('family',
                                'sp_rho', 'sp_pval', 'sp_pval_adj',
                                'ALDEx_rho', 'ALDEx_pval', 'ALDEx_pval_adj')

write_csv(aldex_sp_fam_res, 'Tables/TableS1_DA_fam.csv')

#Scatter plots
#Scatter plot Enterobacteriaceae
comp_En <- ggplot(comp_fam_df, aes(x = Enterobacteriaceae, y = lcn2)) + 
  geom_point(col = "salmon2", size = 2) +
  ylab('LCN-2') + 
  xlab("Enterobacteriaceae Abundance") +
  theme_light()  +
  annotate("text", x = 0.35, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = 0.39, p-val = 0.059")+
  annotate("text", x = 0.35, y = 17000, col = 'black', size = 3,
           label = "ALDEx: rho = 0.35, p-val = 0.105") +
  xlim(0,0.6)
#Scatter plot Oscillospiraceae
comp_Os <- ggplot(comp_fam_df, aes(x = Oscillospiraceae, y = lcn2)) + 
  geom_point(col = "palevioletred", size = 2) +
  ylab('LCN-2') + 
  xlab("Oscillospiraceae Abundance") +
  theme_light()  +
  annotate("text", x = 0.075, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = -0.42, p-val = 0.020")+
  annotate("text", x = 0.075, y = 17000, col = 'black',  size = 3,
           label = "ALDEx: rho = -045, p-val = 0.017") +
  xlim(0,0.125)
#Scatter plot Ruminococcaceae
comp_Ru <- ggplot(comp_fam_df, aes(x = Ruminococcaceae, y = lcn2)) + 
  geom_point(col = "gold2", size = 2) +
  ylab('LCN-2') + 
  xlab("Ruminococcaceae Abundance") +
  theme_light()  +
  annotate("text", x = 0.22, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = -0.41, p-val = 0.033")+
  annotate("text", x = 0.22, y = 17000, col = 'black', size = 3,
           label = "ALDEx: rho = -0.25, p-val = 0.367") +
  xlim(0,0.38)

#Let's merge
scat_fam <- ggarrange(comp_En, comp_Os, comp_Ru, nrow = 1, labels = c('a'),
                      font.label = list(size = 18))

#Let's continue with DAA based on SHAP
xg_fam <- comp_fam_df %>%
  dplyr::select('sample-id', 'lcn2', "clingroup_2", pseq_all_families)
colnames(xg_fam) <- str_replace(colnames(xg_fam), "clingroup_2", "1-No TD; 2 = TD")
colnames(xg_fam) <- str_replace(colnames(xg_fam), "_group", "")

fam <- c(colnames(xg_fam))[-(1:2)]

dataX <- data.matrix(xg_fam[fam])
mod <- xgboost::xgboost(data = dataX, 
                        label = xg_fam$lcn2, 
                        params = list(objective = "reg:squarederror", learning_rate = 1),
                        nrounds = 500)
# To return the SHAP values and ranked features by mean|SHAP|
shap_values <- shap.values(xgb_model = mod, X_train = dataX)
# The ranked features by mean |SHAP|
shap_values$mean_shap_score
# To prepare the long-format data:
shap_long <- shap.prep(xgb_model = mod, X_train = dataX)
# is the same as: using given shap_contrib
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX, top_n = 20)
# **SHAP summary plot**
shap_imp_fam <- shap.plot.summary(shap_long)
# SHAP first 4 var
fig_list <- lapply(names(shap_values$mean_shap_score)[1:4], 
                   shap.plot.dependence, data_long = shap_long)
shap_imp_top_fam <- gridExtra::grid.arrange(grobs = fig_list, ncol = 2)

shap_long_all <- shap.prep(xgb_model = mod, X_train = dataX)
shap_long_all <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)

osc <- shap.plot.dependence(data_long = shap_long_all, x = 'Oscillospiraceae', y = 'Oscillospiraceae', color_feature = "1-No TD; 2 = TD")
des <- shap.plot.dependence(data_long = shap_long_all, x = 'Desulfovibrionaceae', y = 'Desulfovibrionaceae', color_feature = "1-No TD; 2 = TD")
tan <- shap.plot.dependence(data_long = shap_long_all, x = 'Tannerellaceae', y = 'Tannerellaceae', color_feature = "1-No TD; 2 = TD")
bif <- shap.plot.dependence(data_long = shap_long_all, x = 'Bifidobacteriaceae', y = 'Bifidobacteriaceae', color_feature = "1-No TD; 2 = TD")
shap_imp_top_fam <- ggarrange(osc, des, tan, bif, common.legend = T)
#arrange
fam_shap <- ggarrange(shap_imp_fam, shap_imp_top_fam, labels = c('b', 'c'), ncol = 2,
                      widths = c(1.2, 1), heights = c(1.5, 1) ,vjust = 1, font.label = list(size = 18))
ggarrange(scat_fam, fam_shap, ncol = 1, heights = c(1, 1.75))
ggsave('Figures/Fig3_DA_family.tiff', width = 11, height = 10)
```

### At GENUS level

```R
#Let's rename the object, just in case
ps.gen <- pseq

#Let's check the characteristics of this new phyloseq object
microbiome::summarize_phyloseq(ps.gen)

#Composition analysis
#We have to covert the lcn2 variable to numeric
ps.gen@sam_data[["lcn2"]] <- as.numeric(ps.gen@sam_data[["lcn2"]])

# Make sure we use functions from the correct package
transform <- microbiome::transform

#Now, let's see what happens at phylum level
# Merge rare taxa to speed up examples
pseq.comp <- transform(ps.gen, "compositional")
pseq.comp.gen <- aggregate_rare(pseq.comp, level = "Genus", detection = 1/100, prevalence = 10/100)
#56 families in total

#Lets plot the composition sorted by LCN-2
comp_gen_bar <- plot_composition(pseq.comp.gen,
                                 sample.sort = "lcn2",
                                 transform = "compositional") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic")) +
  xlab('Samples sorted by LCN-2 values')

#There are too many genera to see something clear
#Now, let's generate a new tibble for statistical analysis and representation
comp_gen_df <- t(as.data.frame(otu_table(pseq.comp.gen))) %>%
  as.data.frame() %>% 
  rownames_to_column(., var = 'sample-id') %>%
  as_tibble()

#Now, let's join with metadata of interest
comp_gen_df <- dplyr::select(div_analysis, 'sample-id','lcn2', 'clingroup', 'clingroup_2') %>%
  left_join(comp_gen_df, ., by = 'sample-id')

#Conventional spearman
pseq_all_genus <- colnames(comp_gen_df)[-1]
pseq_all_genus <- pseq_all_genus[-(57:69)]

comp_sp_gen = data_frame(genus = character(), rho = double(), pval = double())
for (i in pseq_all_genus){
  # Creating new variables makes the code clearer
  x <- comp_gen_df %>% pull(i)
  x <- as.numeric(x)
  y <- comp_gen_df %>% pull(lcn2)
  cor <- cor.test(x, y, method="spearman")
  rho <- cor$estimate
  pval <- cor$p.value
  comp_sp_gen <- comp_sp_gen %>% add_row('genus' = i, 'rho' = rho, 'pval' = pval)
}
#Let's adjust by bonferroni
p_adj_bon <- comp_sp_gen %>% pull(pval)
p_adj_bon <- p.adjust(p_adj_bon, method = 'bonferroni')
comp_sp_gen$padjusted <- p_adj_bon
#Let's slect those near significan
sig_fam <- comp_sp_gen %>% 
  filter(padjusted < 0.06)
#Three families are associated
sig_fam
#Esch-Sig - 0.44; 0.019

#ALDE-x2
library('ALDEx2')
library(mia)
# summexp
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.fam)
#Agglomerate by genus and subset by prevalence
tse_gen <- mia::subsetByPrevalentTaxa(tse, rank = "Genus", prevalence = 10/100)
#Transform count assay to relative abundances
tse_gen <- mia::transformCounts(tse_gen, method = "relabundance")
detach("package:mia", unload = TRUE)
#Run the test
#lcn2
cont.var <- as.numeric(tse_gen$lcn2)
x <- aldex.clr(assay(tse_gen), tse_gen$lcn2)
res_gen_aldex <- aldex.corr(x, cont.var)
view(res_gen_aldex)
#Esch-Sig - 0.42; 0.054

#Let's do a table for supplementary material
aldex_sp_gen_res <- res_gen_aldex %>%
  rownames_to_column(var = 'genus') %>%
  dplyr::select('genus', 'spearman.erho', 'spearman.ep', 'spearman.eBH') %>%
  as_tibble() %>%
  left_join(comp_sp_gen, ., by = 'genus') %>%
  arrange(., pval)

colnames(aldex_sp_fam_res) <- c('genus',
                                'sp_rho', 'sp_pval', 'sp_pval_adj',
                                'ALDEx_rho', 'ALDEx_pval', 'ALDEx_pval_adj')

write_csv(aldex_sp_gen_res, 'Tables/TableS2_DA_gen.csv')

#Scatter plots
colnames(comp_gen_df) <- str_replace(colnames(comp_gen_df), 'Escherichia-Shigella', 'Escherichia_Shigella')
#Scatter plot Enterobacteriaceae
comp_Esc <- ggplot(comp_gen_df, aes(x = Escherichia_Shigella, y = lcn2)) + 
  geom_point(col = "salmon2", size = 2) +
  ylab('LCN-2') + 
  xlab("Escherichia-Shigella Abundance") +
  theme_light()  +
  annotate("text", x = 0.35, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = 0.44, p-val = 0.019")+
  annotate("text", x = 0.35, y = 17000, col = 'black', size = 3,
           label = "ALDEx: rho = 0.42, p-val = 0.054") +
  xlim(0,0.6)

#let's rename as initially
colnames(comp_gen_df) <- str_replace(colnames(comp_gen_df), 'Escherichia_Shigella', 'Escherichia-Shigella')
#Let's prepare for merging
#Let's generate empty spaces
emp_a <- plot.new()
emp_b <- plot.new()
scat_gen <- ggarrange(comp_Esc, emp_a, nrow = 1, labels = c('a'),
                      font.label = list(size = 18))

#Let's continue with DAA based on SHAP
xg_gen <- comp_gen_df %>%
  dplyr::select('sample-id', 'lcn2', "clingroup_2", pseq_all_genus)
colnames(xg_gen) <- str_replace(colnames(xg_gen), "clingroup_2", "1-No TD; 2 = TD")
colnames(xg_gen) <- str_replace(colnames(xg_gen), "_group", "")

gen <- c(colnames(xg_gen))[-(1:2)]

dataX <- data.matrix(xg_gen[gen])
mod <- xgboost::xgboost(data = dataX, 
                        label = xg_gen$lcn2, 
                        params = list(objective = "reg:squarederror", learning_rate = 1),
                        nrounds = 500)
# To return the SHAP values and ranked features by mean|SHAP|
shap_values <- shap.values(xgb_model = mod, X_train = dataX)
# The ranked features by mean |SHAP|
shap_values$mean_shap_score
# To prepare the long-format data:
shap_long <- shap.prep(xgb_model = mod, X_train = dataX)
# is the same as: using given shap_contrib
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX, top_n = 20)
# **SHAP summary plot**
shap_imp_gen <- shap.plot.summary(shap_long)
# SHAP first 4 var
fig_list <- lapply(names(shap_values$mean_shap_score)[1:4], 
                   shap.plot.dependence, data_long = shap_long)
shap_imp_top_gen <- gridExtra::grid.arrange(grobs = fig_list, ncol = 2)

shap_long_all <- shap.prep(xgb_model = mod, X_train = dataX)
shap_long_all <- shap.prep(shap_contrib = shap_values$shap_score, X_train = dataX)

esc <- shap.plot.dependence(data_long = shap_long_all, x = 'Escherichia-Shigella', y = 'Escherichia-Shigella', color_feature = "1-No TD; 2 = TD")
cri <- shap.plot.dependence(data_long = shap_long_all, x = 'Christensenellaceae_R-7', y = 'Christensenellaceae_R-7', color_feature = "1-No TD; 2 = TD")
eub <- shap.plot.dependence(data_long = shap_long_all, x = '[Eubacterium]_ventriosum', y = '[Eubacterium]_ventriosum', color_feature = "1-No TD; 2 = TD")
rug <- shap.plot.dependence(data_long = shap_long_all, x = '[Ruminococcus]_gnavus', y = '[Ruminococcus]_gnavus', color_feature = "1-No TD; 2 = TD")
shap_imp_top_gen <- ggarrange(esc, cri, eub, rug, common.legend = T)
#arrange
gen_shap <- ggarrange(shap_imp_gen, shap_imp_top_gen, labels = c('b', 'c'),
                      widths = c(1.2, 1), vjust = 1, font.label = list(size = 18))
ggarrange(scat_gen, gen_shap, ncol = 1, heights = c(1, 1.75))
ggsave('Figures/Fig4_DA_gen.tiff', width = 11, height = 10)
```

### Let's try a figure merging Fig3 and Fig4, will be Fig 3 in the manuscript

```R
shap_imp_top2_fam <- ggarrange(osc, des, common.legend = T)
fam_shap_2 <- ggarrange(shap_imp_fam, shap_imp_top2_fam, labels = c('b', 'c'), ncol = 1,
                      heights = c(1.8, 1), vjust = 1, font.label = list(size = 18))

shap_imp_top2_gen <- ggarrange(esc, cri, common.legend = T)
gen_shap_2 <- ggarrange(shap_imp_gen, shap_imp_top2_gen, labels = c('e', 'f'), ncol = 1,
                        heights = c(1.8, 1), vjust = 1, font.label = list(size = 18))

shap_fam_gen <- ggarrange(fam_shap_2, gen_shap_2, ncol = 2)

scat_gen_2 <- ggarrange(comp_Esc, labels = c('d'), font.label = list(size = 18))
sact_fam_gen <- ggarrange(scat_fam, scat_gen_2, ncol = 2, widths = c(3,1))

DA_all <- ggarrange(sact_fam_gen, shap_fam_gen, nrow = 2, heights = c(1, 3))
DA_all
ggsave('Figures/Fig3_4_DA_fam_gen.tiff', width = 14, height = 15)
```
