## A script to measure differential gene expression, gene annotation and pathway analysis

#DE analysis of SRP055063

library(readr)
library(tximport)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)


## Sample information ----

sample = read_csv('SraRunTable.txt') %>% dplyr::select(Run, `SRA Study`, `Sample Name`,`GEO_Accession (exp)`, Experiment, Cell_type, source_name, 
                                              stimulated_with, Platform)
write_csv(sample, "sample_info.csv")

View(sample)

condition = factor(c('control', rep('AGE_BSA', 2), rep('AGE_BSA_BMP2', 2)), 
                    levels = c('control', 'AGE_BSA', 'AGE_BSA_BMP2'))

sample$conditions = condition

#sample$type = factor(c('no_AGE', rep(c('1_AGE_BSA', '5_AGE_BSA'), 2)), levels = c('no_AGE', '1_AGE_BSA', '5_AGE_BSA'))





##tximport input files

#salmon quants file 
sample_files = paste0(pull(sample, Run),'/quant.sf') 
names(sample_files) = pull(sample, Run)   

#gene_map file (gile containing gene id and transcript id of HS)

gene_map = read_csv("/mnt/hdd1/heartp/pro1/down/quants/srp068307/gene_map.csv", col_names = c('enstid', 'ensgid')) 

#loading salmon data using tximport ----
count_data = tximport(files = sample_files, 
                      type = 'salmon', tx2gene = gene_map, 
                      ignoreTxVersion = TRUE)

##analyzing data ----

dds = DESeqDataSetFromTximport(txi = count_data,
                               colData = sample,
                               design = ~ conditions)

dds = DESeq(dds)
dim(dds)


#variance transformation for visualization
vst_dds = varianceStabilizingTransformation(dds)
vst_dds_mat = assay(vst_dds)
plotPCA(vst_dds, intgroup = 'conditions')

resultsNames(dds)

plotMA(dds)


##retrieve and filter results
dds_de = results(dds)
dim(dds_de)
dds_df = as.data.frame(dds_de)
dds_df = rownames_to_column(dds_df, var = 'ensgene')


#filterinf----
#filtering for NA
dds_filter = filter(data.frame(dds_df), complete.cases(dds_df))
dim(dds_filter)

write_csv(dds_filter, 'DEgene_SRP055063.csv')

#volcano plot ----
dds_filter$test = dds_filter$padj <0.05 & abs(dds_filter$log2FoldChange) >1
g = ggplot(dds_filter, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = test), size = 1) +
  theme_bw()


#filtering for pvalue
dds_filter_pv = filter(dds_filter, padj < 0.05)

#filtering for LFC

dds_filter_fc = filter(dds_filter_pv, abs(log2FoldChange) >1)

dim(dds_filter_pv)

write_csv(dds_filter_fc, 'DEgene_filtered_PVFC.csv')

#gene_annotation
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

anno_dds = getBM(attributes = c('ensembl_gene_id',
                                            'chromosome_name',
                                            'start_position',
                                            'end_position',
                                            'strand',
                                            'gene_biotype',
                                            'external_gene_name',
                                            'description'),
                             filters = c('ensembl_gene_id'),
                             values = dds_filter$ensgene, mart = ensembl)

anno_deg_na = left_join(dds_filter, anno_dds, by = c('ensgene' = 'ensembl_gene_id'))

write_csv(anno_deg, 'gene_annotation_filtered.csv')

anno_deg = filter(anno_deg_na, padj < 0.05)
anno_deg = filter(anno_deg, abs(log2FoldChange) > 1)

#heatmap --- 

annot_hm = anno_deg$ensgene

data_anno_hm = vst_dds_mat[annot_hm, ]
rownames(data_anno_hm) = anno_deg$external_gene_name


blues = colorRampPalette(brewer.pal(9, "Blues"))(100)
pheatmap(data_anno_hm, fontsize_row = 4, scale = 'row', color = blues, cutree_cols = 2)


##GO enrichment----

ent_gene = getBM(attributes ='entrezgene_id', filters = 'ensembl_gene_id', 
                 values = anno_deg$ensgene, mart = ensembl)

ent_gene = ent_gene$entrezgene_id
ent_gene = as.character(ent_gene)
ent_uni = getBM(attributes ='entrezgene_id', filters = 'ensembl_gene_id', 
                values = dds_filter$ensgene, mart = ensembl)

ent_uni = ent_uni$entrezgene_id
ent_uni = as.character(ent_uni)

eGO = enrichGO(gene = ent_gene, OrgDb = org.Hs.eg.db,
               ont = "BP",
               universe = ent_uni,
               readable = TRUE)

class(eGO)
write_csv(as.data.frame(eGO), 'GO_enrichment.tsv' )


View(eGO)

barplot(eGO, showCategory = 30)
dotplot(eGO, showCategory = 30)
fc = anno_deg$log2FoldChange
names(fc) = anno_deg$external_gene_name
cnetplot(eGO, showCategory = 20, foldChange = fc)
goplot(eGO)

ekg = enrichKEGG(gene = ent_gene, universe = ent_uni)

write_csv(as.data.frame(ekg), 'KEGG_enrich.csv')

View(as.data.frame(eGO))
