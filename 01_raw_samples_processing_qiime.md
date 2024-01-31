# Code to import and process the raw files

```bash
#Actvate qiime
conda activate qiime2-2023.5
#Verify the qiiime is activated
qiime --help

#Let's go to the WD
cd ${WORKING_DIR}
mkdir qiime_art
mkdir fastq_DV/
#Move the fastq files to fastq_DV directory

##Importing Casava 1.8 paired-end demultiplexed fastq (https://docs.qiime2.org/2022.2/tutorials/importing/)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path fastq_DV \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end_TD.qza

#Let's see quality metrics
qiime demux summarize \
  --i-data demux-paired-end_TD.qza \
  --o-visualization demux-paired-end_TD.qzv
qiime tools view demux-paired-end_TD.qzv

##Using Dada2

##According to the plot, we will trim he sequences from 240 bp
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end_TD.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 287 \
  --p-trunc-len-r 249 \
  --o-representative-sequences qiime_art/dada_rep-seqs.qza \
  --o-table qiime_art/dada_table.qza \
  --o-denoising-stats qiime_art/dada_stats.qza

# To see dada stats
qiime metadata tabulate \
  --m-input-file qiime_art/dada_stats.qza \
  --o-visualization qiime_art/dada_stats.qzv
qiime tools view qiime_art/dada_stats.qzv

# Summarises features table samples
qiime feature-table summarize \
  --i-table qiime_art/dada_table.qza \
  --o-visualization qiime_art/dada_table.qzv \
  --m-sample-metadata-file metadata_TD.tsv
qiime tools view qiime_art/dada_table.qzv

# Summarises 
qiime feature-table tabulate-seqs \
  --i-data qiime_art/dada_rep-seqs.qza \
  --o-visualization qiime_art/dada_rep-seqs.qzv
qiime tools view qiime_art/dada_rep-seqs.qzv

#Let's generate the rooted tree - phylogenetics
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime_art/dada_rep-seqs.qza \
  --o-alignment qiime_art/aligned-rep-seqs.qza \
  --o-masked-alignment qiime_art/masked-aligned-rep-seqs.qza \
  --o-tree qiime_art/unrooted-tree.qza \
  --o-rooted-tree qiime_art/rooted-tree.qza

#Now, clasify the sequences based on taxonomy using the silva classifier (available in qiime website)
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads qiime_art/dada_rep-seqs.qza \
  --o-classification qiime_art/taxonomy.qza

qiime metadata tabulate \
  --m-input-file qiime_art/taxonomy.qza \
  --o-visualization qiime_art/taxonomy.qzv
qiime tools view qiime_art/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime_art/dada_table.qza \
  --i-taxonomy qiime_art/taxonomy.qza \
  --m-metadata-file metadata_TD.tsv \
  --o-visualization qiime_art/taxa-bar-plots_silva.qzv
qiime tools view qiime_art/taxa-bar-plots_silva.qzv


#Moving and renamning
cp qiime_art/rooted-tree.qza qiime_art/dada_table.qza qiime_art/taxonomy.qza .
mv dada_table.qza table.qza
```

#### Now, load into R the following files and generate the phyloseq object:
rooted-tree.qza

metadata_TD.tsv

table.qza

taxonomy.qza



