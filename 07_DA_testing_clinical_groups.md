# Code to perform DA on clinical groups by ANCOM-BC, ALDEx, and ML

Important note, some of the R objects used in this script were generated in the script *06_DA_testing_LCN2.md*

```R
library(tidyverse)
library(microbiome)
library(phyloseq)
library(metagMisc)
library(data.table)
library(ggpubr)
library(vegan)
library('ALDEx2')
library(ANCOMBC)
library(mia)
library("SHAPforxgboost")
library("here")
library(xgboost)
library(patchwork)
```
### At PHYLUM level using ANCOM-BC and ALDEx

```R
#Let's do boxplots to see if there are clear differences

#Three groups
comp_phy_box_1 <- comp_phy_df %>% 
  dplyr::select(., "Bacteroidota","Firmicutes","Proteobacteria",'clingroup', 'clingroup_2') %>%
  pivot_longer(., cols = c("Bacteroidota","Firmicutes","Proteobacteria"),
               names_to = "Phylum", values_to = "abundance",
               names_repair = "minimal") %>%
  ggplot(., aes(x = clingroup, y = abundance, fill = clingroup)) +
  geom_boxplot()+
  facet_wrap(~Phylum)+ scale_fill_hue(l=50, c=100) +
  theme_bw() +
  ylab('Abundance') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

#Two groups
comp_phy_box_2 <- comp_phy_df %>% 
  dplyr::select(., "Bacteroidota","Firmicutes","Proteobacteria",'clingroup', 'clingroup_2') %>%
  pivot_longer(., cols = c("Bacteroidota","Firmicutes","Proteobacteria"),
               names_to = "Phylum", values_to = "abundance",
               names_repair = "minimal") %>%
  ggplot(., aes(x = clingroup_2, y = abundance, fill = clingroup_2)) +
  geom_boxplot()+
  facet_wrap(~Phylum)+ scale_fill_hue(l=90, c=200) +
  theme_bw() +
  ylab('Abundance') + xlab('Group') +
  guides(fill=guide_legend(title="")) + theme(legend.position="none")

box_DA_phy <- ggarrange(comp_phy_box_2, comp_phy_box_1, labels = c('a', 'b'),
          nrow = 2)

#ANCOMBC
#summexp
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
#Three groups
ancom_phy_clingroup <- ancombc2(data = tse, assay_name = "counts", tax_level = "Phylum",
                                fix_formula = "clingroup", rand_formula = NULL,
                                p_adj_method = "holm", 
                                prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                group = "clingroup", struc_zero = TRUE, neg_lb = TRUE,
                                alpha = 0.05, n_cl = 2, verbose = TRUE,
                                global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, 
                                                    verbose = TRUE),
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(),
                                mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                            nrow = 2, 
                                                                            byrow = TRUE),
                                                                     matrix(c(-1, 0, 1, -1),
                                                                            nrow = 2, 
                                                                            byrow = TRUE)),
                                                     node = list(2, 2),
                                                     solver = "ECOS",
                                                     B = 100))

res_ancom_phy_clingroup <- ancom_phy_clingroup$res
#q = 1 in all cases

#Two groups
ancom_phy_clingroup2 <- ancombc2(data = tse, assay_name = "counts", tax_level = "Phylum",
                                 fix_formula = "clingroup_2", rand_formula = NULL,
                                 p_adj_method = "holm", 
                                 prv_cut = 0.10, lib_cut = 500, s0_perc = 0.05,
                                 group = "clingroup_2", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                             nrow = 2, 
                                                                             byrow = TRUE),
                                                                      matrix(c(-1, 0, 1, -1),
                                                                             nrow = 2, 
                                                                             byrow = TRUE)),
                                                      node = list(2, 2),
                                                      solver = "ECOS",
                                                      B = 100))

res_ancom_phy_clingroup2 <- ancom_phy_clingroup2$res
#No differences here neither

#ALDEx
#Three groups - nothing significant
aldex_phy_clingroup <- aldex(assay(tse_phy), tse_phy$clingroup, denom = 'all',
                             test = 'kw')#nothing significant
#Two groups - nothing significant
aldex_phy_clingroup2 <- aldex(assay(tse_phy), tse_phy$clingroup_2, denom = 'all',
                              test = 'kw')#Nothing significant
```

### At FAMILY level using ANCOM-BC and ALDEx

```R
#ANCOMBC
# summexp
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
#Three groups
ancom_fam_clingroup <- ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                                fix_formula = "clingroup", rand_formula = NULL,
                                p_adj_method = "holm", 
                                prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                group = "clingroup", struc_zero = TRUE, neg_lb = TRUE,
                                alpha = 0.05, n_cl = 2, verbose = TRUE,
                                global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, 
                                                    verbose = TRUE),
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(),
                                mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                            nrow = 2, 
                                                                            byrow = TRUE),
                                                                     matrix(c(-1, 0, 1, -1),
                                                                            nrow = 2, 
                                                                            byrow = TRUE)),
                                                     node = list(2, 2),
                                                     solver = "ECOS",
                                                     B = 100))

res_ancom_fam_clingroup <- ancom_fam_clingroup$res
#q = 1 in all cases

#Two groups
ancom_fam_clingroup2 <- ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                                 fix_formula = "clingroup_2", rand_formula = NULL,
                                 p_adj_method = "holm", 
                                 prv_cut = 0.10, lib_cut = 500, s0_perc = 0.05,
                                 group = "clingroup_2", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                             nrow = 2, 
                                                                             byrow = TRUE),
                                                                      matrix(c(-1, 0, 1, -1),
                                                                             nrow = 2, 
                                                                             byrow = TRUE)),
                                                      node = list(2, 2),
                                                      solver = "ECOS",
                                                      B = 100))

res_ancom_fam_clingroup2 <- ancom_fam_clingroup2$res
#No differences here neither

#ALDEx
#Three groups - nothing significant
aldex_fam_clingroup <- aldex(assay(tse_fam), tse_phy$clingroup, denom = 'all',
                             test = 'kw')#nothing significant
#Two groups - nothing significant
aldex_fam_clingroup2 <- aldex(assay(tse_fam), tse_phy$clingroup_2, denom = 'all',
                              test = 'kw')#Nothing significant
```

### At GENUS level using ANCOM-BC and ALDEx

```R
#ANCOMBC
# summexp
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
#Three groups
ancom_gen_clingroup <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                                fix_formula = "clingroup", rand_formula = NULL,
                                p_adj_method = "holm", 
                                prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                group = "clingroup", struc_zero = TRUE, neg_lb = TRUE,
                                alpha = 0.05, n_cl = 2, verbose = TRUE,
                                global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, 
                                                    verbose = TRUE),
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(),
                                mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                            nrow = 2, 
                                                                            byrow = TRUE),
                                                                     matrix(c(-1, 0, 1, -1),
                                                                            nrow = 2, 
                                                                            byrow = TRUE)),
                                                     node = list(2, 2),
                                                     solver = "ECOS",
                                                     B = 100))

res_ancom_gen_clingroup <- ancom_gen_clingroup$res
#q = 1 in all cases

#Two groups
ancom_gen_clingroup2 <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                                 fix_formula = "clingroup_2", rand_formula = NULL,
                                 p_adj_method = "holm", 
                                 prv_cut = 0.10, lib_cut = 500, s0_perc = 0.05,
                                 group = "clingroup_2", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                             nrow = 2, 
                                                                             byrow = TRUE),
                                                                      matrix(c(-1, 0, 1, -1),
                                                                             nrow = 2, 
                                                                             byrow = TRUE)),
                                                      node = list(2, 2),
                                                      solver = "ECOS",
                                                      B = 100))

res_ancom_gen_clingroup2 <- ancom_gen_clingroup2$res
#No differences here neither

#ALDEx
#Three groups - nothing significant
aldex_gen_clingroup <- aldex(assay(tse_gen), tse_phy$clingroup, denom = 'all',
                             test = 'kw')#nothing significant
#Two groups - nothing significant
aldex_gen_clingroup2 <- aldex(assay(tse_gen), tse_phy$clingroup_2, denom = 'all',
                              test = 'kw')#Nothing significant
```
### Nothing significant in any DA test

#### Nevertheless, let's see what happens using the SHAP values
Seeing no differences, here-in-after we will compare only two groups

### Phylum

```R
xg_phy <- comp_phy_df %>%
  dplyr::select('sample-id', 'lcn2', "clingroup_2", "Actinobacteriota", 
                "Bacteroidota", "Desulfobacterota", "Firmicutes", "Other", 
                "Proteobacteria", "Verrucomicrobiota")
xg_phy$clingroup_2 <- str_replace(xg_phy$clingroup_2, '1', '0')
xg_phy$clingroup_2 <- str_replace(xg_phy$clingroup_2, '2', '1')
colnames(xg_phy) <- str_replace(colnames(xg_phy), "clingroup_2", "0-No TD; 1 = TD")
colnames(xg_phy) <- str_replace(colnames(xg_phy), "lcn2", "LCN2")

phy <- c(colnames(xg_phy))[-(1)]
phy <- phy[-2]
  
dataX <- data.matrix(xg_phy[phy])
mod <- xgboost::xgboost(data = dataX, 
                        label = xg_phy$`0-No TD; 1 = TD`, 
                        params = list(objective = "binary:logistic", learning_rate = 1),
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

fir <- shap.plot.dependence(data_long = shap_long, x = 'Firmicutes', y = 'Firmicutes', color_feature = "LCN2")
bac <- shap.plot.dependence(data_long = shap_long, x = 'Bacteroidota', y = 'Bacteroidota', color_feature = "LCN2")
lcn <- shap.plot.dependence(data_long = shap_long, x = 'LCN2', y = 'LCN2', color_feature = "LCN2")
act <- shap.plot.dependence(data_long = shap_long, x = 'Actinobacteriota', y = 'Actinobacteriota', color_feature = "LCN2")

shap_imp_top_phy <- ggarrange(fir, bac, lcn, act, common.legend = T)
#arrange
phy_shap_TD <- ggarrange(shap_imp_phy, shap_imp_top_phy, labels = c('e', 'f'),
                         ncol = 2, vjust = 1, font.label = list(size = 18))
ggarrange(box_DA_phy, phy_shap_TD, nrow = 2, heights = c(1, 1))
ggsave('Figures/FigS3_DA_phylum_TD.tiff', width = 11, height = 11)
```
 
### Family

```R
xg_fam <- comp_fam_df %>%
  dplyr::select('sample-id', 'lcn2', "clingroup_2", pseq_all_families)
colnames(xg_fam) <- str_replace(colnames(xg_fam), "_group", "")
xg_fam$clingroup_2 <- str_replace(xg_fam$clingroup_2, '1', '0')
xg_fam$clingroup_2 <- str_replace(xg_fam$clingroup_2, '2', '1')
colnames(xg_fam) <- str_replace(colnames(xg_fam), "clingroup_2", "0-No TD; 1 = TD")
colnames(xg_fam) <- str_replace(colnames(xg_fam), "lcn2", "LCN2")

fam <- c(colnames(xg_fam))[-1]
fam <- fam[-2]

dataX <- data.matrix(xg_fam[fam])
mod <- xgboost::xgboost(data = dataX, 
                        label = xg_phy$`0-No TD; 1 = TD`, 
                        params = list(objective = "binary:logistic", learning_rate = 1),
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

rik <- shap.plot.dependence(data_long = shap_long, x = 'Rikenellaceae', y = 'Rikenellaceae', color_feature = "LCN2")
bac <- shap.plot.dependence(data_long = shap_long, x = 'Bacteroidaceae', y = 'Bacteroidaceae', color_feature = "LCN2")
ent <- shap.plot.dependence(data_long = shap_long, x = 'Enterobacteriaceae', y = 'Enterobacteriaceae', color_feature = "LCN2")
lac <- shap.plot.dependence(data_long = shap_long, x = 'Lachnospiraceae', y = 'Lachnospiraceae', color_feature = "LCN2")

shap_imp_top_fam <- ggarrange(rik, bac, ent, lac, common.legend = T)
#arrange
shap_fam_TD <- ggarrange(shap_imp_fam, shap_imp_top_fam, labels = c('a', 'b'), ncol = 2,
                      widths = c(1.2, 1), heights = c(1.5, 1) ,vjust = 1, font.label = list(size = 18))
```

### Genus

```R
xg_gen <- comp_gen_df %>%
  dplyr::select('sample-id', 'lcn2', "clingroup_2", pseq_all_genus)
colnames(xg_gen) <- str_replace(colnames(xg_gen), "_group", "")
xg_gen$clingroup_2 <- str_replace(xg_gen$clingroup_2, '1', '0')
xg_gen$clingroup_2 <- str_replace(xg_gen$clingroup_2, '2', '1')
colnames(xg_gen) <- str_replace(colnames(xg_gen), "clingroup_2", "0-No TD; 1 = TD")
colnames(xg_gen) <- str_replace(colnames(xg_gen), "lcn2", "LCN2")

gen <- c(colnames(xg_gen))[-(1)]
gen <- gen[-2]

dataX <- data.matrix(xg_gen[gen])
mod <- xgboost::xgboost(data = dataX, 
                        label = xg_phy$`0-No TD; 1 = TD`, 
                        params = list(objective = "binary:logistic", learning_rate = 1),
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

ali <- shap.plot.dependence(data_long = shap_long, x = 'Alistipes', y = 'Alistipes', color_feature = "LCN2")
esc <- shap.plot.dependence(data_long = shap_long, x = 'Escherichia-Shigella', y = 'Escherichia-Shigella', color_feature = "LCN2")
bac <- shap.plot.dependence(data_long = shap_long, x = 'Bacteroides', y = 'Bacteroides', color_feature = "LCN2")
lac <- shap.plot.dependence(data_long = shap_long, x = 'Lachnospira', y = 'Lachnospira', color_feature = "LCN2")
shap_imp_top_gen <- ggarrange(ali, esc, bac, lac, common.legend = T)
#arrange
shap_gen_TD <- ggarrange(shap_imp_gen, shap_imp_top_gen, labels = c('c', 'd'), ncol = 2,
                         widths = c(1.2, 1), heights = c(1.5, 1) ,vjust = 1, font.label = list(size = 18))
ggarrange(shap_fam_TD, shap_gen_TD, ncol = 1, heights = c(1, 1))
ggsave('Figures/FigS4_DA_fam_gen_TD.tiff', width = 12, height = 12)
```
