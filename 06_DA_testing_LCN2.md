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
library(umap)
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

#Now, lets see what happens at phylim level
# Merge rare taxa to speed up examples
pseq.comp <- transform(ps.phy, "compositional")
pseq.comp.phy <- aggregate_rare(pseq.comp, level = "Phylum", detection = 1/100, prevalence = 10/100)

#LCN-2
#Lets plot the composition sorted by LCN-2
comp_phy_bar <- plot_composition(pseq.comp.phy,
                                 sample.sort = "lcn2",
                                 transform = "compositional") +
  scale_fill_brewer("Phylum", palette = "Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic")) +
  xlab('Samples sorted by LCN-2 values')

#There is some evidence of association between phylum and lcn-2
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
library('ALDEx2')
library(mia)
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
scat <- ggarrange(comp_P, comp_F, comp_B, nrow = 1)
plot_phy <- ggarrange(bar, scat, ncol = 1, labels = c('a', 'b'),
                      font.label = list(size = 15))
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

# Make sure we use functions from correct package
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

#There are too many families to see something celear
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

#Let's slect those near significan
sig_fam <- comp_sp_fam %>% 
  filter(padjusted < 0.06)
#Three families are associated
sig_fam
#Enterobac - 0.39; 0.059, Oscillospi - -0.42, 0.020; Ruminococ - -0.41, 0.033

#ALDE-x2
library('ALDEx2')
library(mia)
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
aldex_sp_fam_res %>% dplyr::select('family',
                            'sp_rho', 'sp_pval', 'sp_pval_adj',
                            'ALDEx_rho', 'ALDEx_pval_adj') %>%
  write_csv(., 'Tables/Supplementary Table S1.csv')

#Scatter plots
#Scatter plot Oscillospiraceae
comp_Os <- ggplot(comp_fam_df, aes(x = Oscillospiraceae, y = lcn2)) + 
  geom_point(col = "orchid", size = 2) +
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
  geom_point(col = "violetred", size = 2) +
  ylab('LCN-2') + 
  xlab("Ruminococcaceae Abundance") +
  theme_light()  +
  annotate("text", x = 0.22, y = 18500, col = 'black', size = 3,
           label = "Sp cor: rho = -0.41, p-val = 0.033")+
  annotate("text", x = 0.22, y = 17000, col = 'black', size = 3,
           label = "ALDEx: rho = -0.25, p-val = 0.367") +
  xlim(0,0.38)

#Let's merge
scat_fam <- ggarrange(comp_Os, comp_Ru, nrow = 1, labels = c('c'),
                      font.label = list(size = 15))
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

# Make sure we use functions from correct package
transform <- microbiome::transform

#Now, lets see what happens at phylim level
# Merge rare taxa to speed up examples
pseq.comp <- transform(ps.gen, "compositional")
pseq.comp.gen <- aggregate_rare(pseq.comp, level = "Genus", detection = 1/100, prevalence = 10/100)
#56 families in total

#LCN-2
#Lets plot the composition sorted by LCN-2
comp_gen_bar <- plot_composition(pseq.comp.gen,
                                 sample.sort = "lcn2",
                                 transform = "compositional") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(face = "italic")) +
  xlab('Samples sorted by LCN-2 values')

#There are too many genus to see something clear
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

colnames(aldex_sp_gen_res) <- c('genus',
                                'sp_rho', 'sp_pval', 'sp_pval_adj',
                                'ALDEx_rho', 'ALDEx_pval', 'ALDEx_pval_adj')
aldex_sp_gen_res %>% dplyr::select('genus',
                                   'sp_rho', 'sp_pval', 'sp_pval_adj',
                                   'ALDEx_rho', 'ALDEx_pval_adj') %>%
  write_csv(., 'Tables/Supplementary Table S2.csv')

#Scatter plots
colnames(comp_gen_df) <- str_replace(colnames(comp_gen_df), 'Escherichia-Shigella', 'Escherichia_Shigella')
#Scatter plot Enterobacteriaceae
comp_Esc <- ggplot(comp_gen_df, aes(x = Escherichia_Shigella, y = lcn2)) + 
  geom_point(col = "darkorange3", size = 2) +
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
scat_gen <- ggarrange(comp_Esc, nrow = 1, labels = c('d'),
                      font.label = list(size = 15))
```

### Let's do Fig 2

```R
scat_fam_gen <- ggarrange(scat_fam, scat_gen, ncol = 2, widths = c(2,1))

ggarrange(plot_phy, scat_fam_gen, ncol = 1, heights = c(2,1))

ggsave('Figures/Fig2_DA_LCN2.tiff', width = 11, height = 11)
```

## UMAP for association LCN2 - bacteria abundances

```R
#First, let's generate a database with all the features
out_var <- comp_phy_df %>% dplyr::select('sample-id', 'lcn2', 'clingroup', 'clingroup_2')

comp_phy_df_no_out <- subset(comp_phy_df, select = -c(lcn2, clingroup, clingroup_2))
comp_fam_df_no_out <- subset(comp_fam_df, select = -c(lcn2, clingroup, clingroup_2))
comp_gen_df_no_out <- subset(comp_gen_df, select = -c(lcn2, clingroup, clingroup_2))

all_bact_df <- comp_phy_df_no_out %>%
  dplyr::left_join(., comp_fam_df_no_out, by = 'sample-id') %>%
  dplyr::left_join(., comp_gen_df_no_out, by = 'sample-id')

all_bact_df <- all_bact_df[,-1]

#Now, let's see what bacteria are significant for both methods

#First none-adjausted p-val
noadj_phy <- c('Firmicutes', 'Proteobacteria')

noadj_fam <- aldex_sp_fam_res %>% 
  filter(sp_pval < 0.05 | ALDEx_pval < 0.05)
noadj_fam <- as.vector(noadj_fam$family)
noadj_fam <- str_replace(noadj_fam, 'Clostridia_UCG-014', 'Clostridia_UCG-014.x')

noadj_gen <- aldex_sp_gen_res %>% 
  filter(sp_pval < 0.05 | ALDEx_pval < 0.05)
noadj_gen <- as.vector(noadj_gen$genus)
noadj_gen <- str_replace(noadj_gen, 'Clostridia_UCG-014', 'Clostridia_UCG-014.y')

noadj_all <- c(noadj_phy, noadj_fam, noadj_gen)

#Second adjusted p-val
adj_phy <- c('Firmicutes', 'Proteobacteria')

adj_fam <- aldex_sp_fam_res %>% 
  filter(sp_pval_adj < 0.05 | ALDEx_pval_adj < 0.05)
adj_fam <- as.vector(adj_fam$family)

adj_gen <- aldex_sp_gen_res %>% 
  filter(sp_pval_adj < 0.05 | ALDEx_pval_adj < 0.05)
adj_gen <- as.vector(adj_gen$genus)

adj_all <- c(adj_phy, adj_fam, adj_gen)

#Now, let's select the different databases and add lcn2 - 
#ALL bacteria
all_bact_df_norm <- scale(all_bact_df)
all_bact.umap <- umap(all_bact_df_norm)
all_bact.umap_df <- as_tibble(all_bact.umap$layout)

all_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(all_bact.umap_df, .)%>%
  as_tibble()

all_bact.umap_df <- all_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))

norm_all_log <- ggplot(all_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()

#Non adjusted p-val
noadj_bact_df <- subset(all_bact_df, select=noadj_all)
noadj_bact_df_norm <- scale(noadj_bact_df)
noadj_bact.umap <- umap(noadj_bact_df_norm)
noadj_bact.umap_df <- as_tibble(noadj_bact.umap$layout)

noadj_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(noadj_bact.umap_df, .)%>%
  as_tibble()

noadj_bact.umap_df <- noadj_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))

norm_noadj_log <- ggplot(noadj_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()

#Adjusted p-val
adj_bact_df <- subset(all_bact_df, select=adj_all)
adj_bact_df_norm <- scale(adj_bact_df)
adj_bact.umap <- umap(adj_bact_df_norm)
adj_bact.umap_df <- as_tibble(adj_bact.umap$layout)

adj_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(adj_bact.umap_df, .)%>%
  as_tibble()

adj_bact.umap_df <- adj_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))

norm_adj_log <- ggplot(adj_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()

############
norm_log <- ggarrange(norm_all_log, norm_noadj_log, norm_adj_log, ncol = 1,
                        labels = c('b', 'd', 'f'), common.legend = T)
#####################

#Let's try with sp and Aldex separately

#First none-adjausted p-val for sp
noadj_phy <- c('Firmicutes', 'Proteobacteria')

noadj_fam <- aldex_sp_fam_res %>% 
  filter(sp_pval < 0.05)
noadj_fam <- as.vector(noadj_fam$family)
noadj_fam <- str_replace(noadj_fam, 'Clostridia_UCG-014', 'Clostridia_UCG-014.x')

noadj_gen <- aldex_sp_gen_res %>% 
  filter(sp_pval < 0.05)
noadj_gen <- as.vector(noadj_gen$genus)
noadj_gen <- str_replace(noadj_gen, 'Clostridia_UCG-014', 'Clostridia_UCG-014.y')

noadj_all_sp <- c(noadj_phy, noadj_fam, noadj_gen)

#Second adjusted p-val - Sp
adj_phy <- c('Firmicutes', 'Proteobacteria')

adj_fam <- aldex_sp_fam_res %>% 
  filter(sp_pval_adj < 0.05)
adj_fam <- as.vector(adj_fam$family)

adj_gen <- aldex_sp_gen_res %>% 
  filter(sp_pval_adj < 0.05)
adj_gen <- as.vector(adj_gen$genus)

adj_all_sp <- c(adj_phy, adj_fam, adj_gen)

#
#No adjusted p-val - Sp
noadj_bact_df <- scale(subset(all_bact_df, select=noadj_all_sp))
noadj_bact.umap <- umap(noadj_bact_df)
noadj_bact.umap_df <- as_tibble(noadj_bact.umap$layout)

noadj_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(noadj_bact.umap_df, .)%>%
  as_tibble()

noadj_bact.umap_df <- noadj_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))

Sp_noadj_log <- ggplot(noadj_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()

# P-Val adjusted bacteria
adj_bact_df <- scale(subset(all_bact_df, select=adj_all_sp))
adj_bact.umap <- umap(adj_bact_df, n_neighbors = 15)
adj_bact.umap_df <- as_tibble(adj_bact.umap$layout)

adj_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(adj_bact.umap_df, .)%>%
  as_tibble()

adj_bact.umap_df <- adj_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))

Sp_adj_log <- ggplot(adj_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()

#Second none-adjausted p-val for ALDEx
noadj_phy <- c('Proteobacteria')

noadj_fam <- aldex_sp_fam_res %>% 
  filter(ALDEx_pval < 0.05)
noadj_fam <- as.vector(noadj_fam$family)

noadj_gen <- aldex_sp_gen_res %>% 
  filter(ALDEx_pval < 0.05)
noadj_gen <- as.vector(noadj_gen$genus)

noadj_all_ALDEx <- c(noadj_phy, noadj_fam, noadj_gen)

#Second adjusted p-val - ALDEx
adj_phy <- c('Proteobacteria')

adj_fam <- aldex_sp_fam_res %>% 
  filter(ALDEx_pval_adj < 0.05)
adj_fam <- as.vector(adj_fam$family)

adj_gen <- aldex_sp_gen_res %>% 
  filter(ALDEx_pval_adj < 0.05)
adj_gen <- as.vector(adj_gen$genus)

adj_all_ALDEx <- c(adj_phy, adj_fam, adj_gen)

#
#No adjusted p-val - ALDEx
noadj_bact_df <- scale(subset(all_bact_df, select=noadj_all_ALDEx))
noadj_bact.umap <- umap(noadj_bact_df)
noadj_bact.umap_df <- as_tibble(noadj_bact.umap$layout)

noadj_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(noadj_bact.umap_df, .)%>%
  as_tibble()

noadj_bact.umap_df <- noadj_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))


ALDEx_noadj_log <- ggplot(noadj_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()

# P-Val adjusted bacteria
adj_bact_df <- scale(subset(all_bact_df, select=adj_all_ALDEx))
adj_bact.umap <- umap(adj_bact_df, n_neighbors = 12)
adj_bact.umap_df <- as_tibble(adj_bact.umap$layout)

adj_bact.umap_df <- out_var %>%
  dplyr::select(., 'lcn2') %>%
  cbind(adj_bact.umap_df, .)%>%
  as_tibble()

adj_bact.umap_df <- adj_bact.umap_df %>%
  mutate(log_lcn2 = log(lcn2))

ALDEx_adj_log <- ggplot(adj_bact.umap_df, aes(x=V1, y=V2, color = log_lcn2)) +
  geom_point(size = 2.5) +
  labs(color = "log(LCN2)") +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_viridis()
#LEt's do two extra plots to see if we can say something else
adj_bact.umap_df_abun <- as_tibble(cbind(adj_bact.umap_df, subset(all_bact_df, select=adj_all_ALDEx)))

ALDEx_adj_log_ab <- ggplot(adj_bact.umap_df_abun, aes(x=V1, y=V2, color = log(Proteobacteria), size = Oscillospiraceae)) +
  geom_point(alpha = 0.7) +
  labs(color = "Proteobacteria", size = 'Oscillospiraceae') +
  theme(legend.title = element_text()) +
  theme_bw() +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA)

#For ploting, we will select 
top_plots <- ggarrange(norm_all_log, Sp_noadj_log, norm_adj_log, common.legend = T, 
                       labels = c('a', 'b', 'c'), legend = 'right', nrow = 1)
botton_plots <- ggarrange(ALDEx_adj_log, ALDEx_adj_log_ab, labels = c('d', 'e'), 
                          widths = c(1,1.0))

#LEt's do the figure
ggarrange(top_plots, botton_plots, nrow = 2, heights = c(1,1.2))
ggsave('Figures/Fig3_UMAP.tiff', width = 10, height = 7)
```
