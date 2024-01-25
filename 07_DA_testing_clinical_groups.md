# Code to perform DA on clinical groups by ANCOM-BC and ALDEx

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
```
## Let's do a boxplot to explore the differences

```R

```

## Let's do DA testing

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
