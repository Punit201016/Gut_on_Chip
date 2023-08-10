#### Mittia sample new march 2023 ###

library(DOSE)
library(EnhancedVolcano)
library(tidyverse)
library(readxl)
library(stringi)
library(edgeR)
library(datapasta)
library(ggplot2)
library(biomaRt)
library(RSQLite)
library(enrichplot)
library(biomaRt)



counts_matia_chip_new <- read.delim("C:/Users/Tyagi/Desktop/matia_chip/counts_matia_chip_new.txt", comment.char="#")

RNAseq_manifest_Mattia <- read_excel("C:/Users/Tyagi/Desktop/RNAseq manifest Mattia.xlsx")

colnames(counts_all_carlota) <- sub("Aligned.*", "", colnames(counts_all_carlota))

counts_mattia <- as.data.frame(counts_matia_chip_new[, c(1, 7:23)])


Epi_in <- RNAseq_manifest_Mattia[c(1:12, 14, 15, 16, 17,18), c(1:5)]


Epi_in$group <- apply(Epi_in[, c(1, 3, 4, 5)],1,paste,collapse="_") ## add all column in data frame

##remove other text in the header in the count
colnames(counts_mattia) <- sub(".Aligned.*", "", colnames(counts_mattia)) 



Experiment_time <- factor(x = Epi_in$`Experiment number`)
Condition_1 <- factor(x = Epi_in$Response)
Treatment_type <- factor(x = Epi_in$`Type of treatment`)

epi_group <- as.factor(Epi_in$group)

y_feature_1 <- DGEList(counts = data.matrix(feature_counts_1))

design_new_epi <- model.matrix(~0+Replicate, data = feature_counts_1)


feature_counts_1 <- column_to_rownames(counts_mattia, var = "Geneid")

rownames(feature_counts_1) <- feature_counts_1[,1]
feature_counts_1 <- feature_counts_1[, -1]

feature_counts_1[is.na(feature_counts_1)] <- 0




drop <- filterByExpr(y_feature_1)
y_feature_1 <- y_feature_1[drop, , keep.lib.sizes=FALSE]


dim(y_feature_1)

### Normalizing the data
y_feature_1 <- calcNormFactors(y_feature_1)


##### MDS plot ####
nb.cols <- 17
col.treatment <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols) ### Putting the colors to the treatments ##
data.frame(epi_group, col.treatment)
mds <- plotMDS(y_feature_1, col = col.treatment, labels = epi_group)
legend("topleft", cex = 0.5, fill=col.treatment,legend=levels(epi_group))
title("MDS plot Treatment")

## using ggplot
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = epi_group)
library(ggplot2)
p1 <- ggplot(toplot, aes(Dim1, Dim2, colour = epi_group)) + 
  geom_point(size=12, stroke = 2, aes(shape=epi_group, color=col.treatment, size=epi_group))+ 
  theme(legend.position = "none", text = element_text(size = 20))+
  scale_shape_manual(values = seq(0,17))



### PCA
logCPM <- cpm(y_feature_1, log=TRUE)

# Calculate rowwise variance
rv <- apply(logCPM, 1, var)
# Sort decreasingly and take top 1000
o <- order(rv, decreasing=TRUE)
top4000 <- head(o, 4000)
# From the logCPMs subset for the top-1000
logCPM_top4000 <- logCPM[top4000,]
# Run PCA
pca <- prcomp(t(logCPM_top4000))
# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y_feature_1$samples)
# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(PC1, PC2, colour=epi_group)) + 
  geom_point(size=8, stroke = 2, aes(shape=Condition_1, color=Experiment_time))+ 
  theme(legend.position = "right", text = element_text(size = 12))+
  scale_shape_manual(values = seq(0,17))+
  xlab(labs[1]) + 
  ylab(labs[2])+
  labs(shape="Group", color="Experiment_time(Batch)")


### Batch correction

logCPMs_corrected <- limma::removeBatchEffect(logCPM, batch = Experiment_time)
# repeat PCA as before, using the same genes
logCPM_corrected_top4000 <- logCPMs_corrected[top4000,]
# Run PCA
pca <- prcomp(t(logCPM_corrected_top4000))
# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y_feature_1$samples)
# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(PC1, PC2, colour=epi_group))+ 
  geom_point(size=8, stroke = 2, aes(shape=Condition_1, color=Experiment_time))+ 
  theme(legend.position = "right", text = element_text(size = 12))+
  scale_shape_manual(values = seq(0,17))+
  xlab(labs[1]) + 
  ylab(labs[2])+
  labs(shape="Group", color="Experiment_time(Batch)")
  

counts_mattia_1 <- counts_mattia %>% dplyr::select(contains(c()))


####### DEGs #####################################

V1_voom <- voom(y_feature_1, design_new_epi, plot = T)

corfit<- duplicateCorrelation(V1_voom, design_new_epi, block = Treatment_type)

V1_voom <- voom(y_feature_1, design_new_epi, 
                block = Treatment_type, 
                correlation = corfit$consensus.correlation, plot = T)



fit_1_limma <- lmFit(V1_voom, design_new_epi)


contrast_limma <- makeContrasts(Treatment_typeSupernatant - Treatment_typeLive.bacteria, 
                                levels = colnames(coef(fit_1_limma))) #### estimating contrast using make contrast function

est_contr <- contrasts.fit(fit_1_limma, contrast_limma) 

est_contr <- eBayes(est_contr) ### smoothing the standard error

result_1 <- decideTests(est_contr)


summary(result_1)




top_table_limma <- topTable(est_contr, sort.by = "P", adjust.method = "fdr", number = Inf)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

top_table_1 <- rownames_to_column(top_table_limma)

colnames(top_table_1)[1] <- "ensembl_gene_id"

cntrl_vs_nsc <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                              "description", "uniprotswissprot", "entrezgene_id"),
                               filters = "ensembl_gene_id", values = top_table_1,
                               mart = ensembl, uniqueRows = TRUE)
cntrl_vs_nsc_merge <- merge(top_table_1, cntrl_vs_nsc, by= "ensembl_gene_id", no.dups = TRUE)
cntrl_vs_nsc_merge <- cntrl_vs_nsc_merge[!duplicated(cntrl_vs_nsc_merge$ensembl_gene_id), ] ##### remove duplicate genes


volcano_data <- as.data.frame(cntrl_vs_nsc_merge[, c(2,6)], row.names = cntrl_vs_nsc_merge$external_gene_name)


EnhancedVolcano(toptable = volcano_data, lab = rownames(volcano_data), 
                x = "logFC", y = "adj.P.Val", pCutoff = 0.05, titleLabSize = 20, 
                title = "Non-Responder vs Control", FCcutoff = 1.2, labSize = 5) + 
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) + scale_y_continuous(n.breaks = 10)


volcano_data$LogFC <- as.numeric(volcano_data$LogFC)
volcano_data$padj <- as.numeric(volcano_data$padj)


EnhancedVolcano(volcano_data,
                lab = rownames(volcano_data),
                x = 'logFC',
                y = 'adj.P.Val',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 2.0,
                labSize = 6.0,
                title = "Non-Responder vs Act_control", 
                colAlpha = 1,
                legendPosition = 'left',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, 
                max.overlaps = 100)


cntrl_vs_nsc_merge <- CREBBP_df

###Filtering the data according to the logFC 1.2 and adj.pvalue 0.05 
Filter_data <-dplyr::select(.data = cntrl_vs_nsc_merge, ensembl_gene_id, logFC, adj.P.Val, external_gene_name, entrezgene_id)
Filter_data <- dplyr::filter(.data = Filter_data, logFC > 0 | logFC < -(0), adj.P.Val < 0.05) 

Filter_data_up <- dplyr::filter(.data = Filter_data, logFC > 0, adj.P.Val < 0.05) 
Filter_data_down <- dplyr::filter(.data = Filter_data, logFC < -(0), adj.P.Val < 0.05) 



count_data <- rownames_to_column(all_cntrls_1_3)
#z-score
logCPM <- cpm(y_feature_1, log=TRUE) 
z_score <- t(scale(t(logCPM)))
z_score <- rownames_to_column(as.data.frame(z_score))

out <- list("sheet1" = cntrl_vs_nsc_merge, "sheet2"= z_score, 
            "sheet3"= Filter_data, "sheet4"= Filter_data_up,"sheet5"= Filter_data_down,
            "sheet6"= count_data)


writexl::write_xlsx(out, 
                    "/Users/Tyagi/Desktop/pca_matia/pca_after_comments/Analysed_data/R_vs_Act_cntrl_Expt_1.xlsx")


write.csv(x = Epi_in, file = "/Users/Tyagi/Desktop/matia_new_metadat.csv")


go_uni <- as.character(Filter_data$ensembl_gene_id)


###GO enrichment using cluster profiler

geneList <- Filter_data[, 2] ### logFc

names(geneList) = as.character(Filter_data[, 1]) ### naming logFc

geneList <- sort(geneList, decreasing = TRUE)


### enrichment analysis using ORA
GO_enrich_1_tmrvsFT_up <- clusterProfiler::enrichGO(Filter_data_up$ensembl_gene_id, OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL', ont = "ALL", readable = TRUE, pvalueCutoff = 0.5, qvalueCutoff = 0.2, pAdjustMethod = "fdr")

GO_enrich_up_simp_up <- clusterProfiler::simplify(GO_enrich_1_tmrvsFT_up, cutoff=0.9, by="p.adjust", select_fun= min)


ego_up <- mutate(GO_enrich_up_simp_up, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))


ggplot(ego_up, showCategory = 20, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(1, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("Biological Processes NR vs Control up-regulated genes")



edox2 <- pairwise_termsim(GO_enrich_1_tmrvsFT_up, method = "JC")
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average", showCategory = 20, nCluster =10)
aplot::plot_list(p2, tag_levels='A')


write.csv(x = GO_enrich_down_simp, file = "/Users/Tyagi/Desktop/matia/RvsCntrl/GO_up.csv")



#### gene set enrichment analysis
gse_cnt_vs_GZMK_trt <- clusterProfiler::gseGO(geneList, OrgDb = "org.Hs.eg.db", 
                                              keyType = 'ENSEMBL', ont = "ALL", 
                                              nPermSimple = 100000, eps = 0, 
                                              pvalueCutoff = 0.1, pAdjustMethod = "fdr")

Gsea_enrich_up_simp <- clusterProfiler::simplify(gse_cnt_vs_GZMK_trt, 
                                                 cutoff=0.6, by="p.adjust", 
                                                 select_fun= min)

Gsea_enrich_up_simp <- clusterProfiler::setReadable(x = gse_cnt_vs_GZMK_trt, 
                                                    OrgDb = "org.Hs.eg.db", 
                                                    keyType = 'ENSEMBL')

clusterProfiler::dotplot(Gsea_enrich_up_simp, 
                         showCategory =10, split=".sign", 
                         title = "NR vs Control") + 
                          facet_grid(.~.sign) #plot to show activated and suppressed process

## GSEA plot ##
ewp_gse_gzmk <- arrange(Gsea_enrich_up_simp, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:5)

ggplot(ewp_gse_gzmk, showCategory=10,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("NR vs Control GSEA")



##### Pathway enrichment analysis #####
geneList2 <- Filter_data[, 2] ### logFc
names(geneList2) = as.character(Filter_data[, 5]) ### naming logFc
geneList2 <- sort(geneList2, decreasing = TRUE)


###### Pathway enrichment analysis using enrichkegg ###
enrich_kegg_cnt_vs_GZMK_up <- clusterProfiler::enrichKEGG(Filter_data_up$entrezgene_id, keyType = "kegg", organism = "hsa", pvalueCutoff = 0.1, pAdjustMethod = "fdr", qvalueCutoff = 0.2)

enrich_kegg_up <- clusterProfiler::setReadable(x = enrich_kegg_cnt_vs_GZMK_up, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")


enrich_kegg_up <- mutate(enrich_kegg_up, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_kegg_up, showCategory = 10, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("KEGG pathway NR vs Control up regulated genes")




###### Pathway enrichment analysis using gsekegg ###
gse_kegg_cnt_vs_GZMK <- clusterProfiler::gseKEGG(geneList= geneList2, organism = "hsa", pvalueCutoff = 0.1, verbose= FALSE, keyType = "kegg", nPermSimple = 100000, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

gse_kegg_cnt_vs_GZMK <- clusterProfiler::setReadable(x = gse_kegg_cnt_vs_GZMK, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

## GSEA plot ##
ewp_gse_kegg_gzmk <- arrange(gse_kegg_cnt_vs_GZMK, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>% dplyr::slice(1:10)

ggplot(ewp_gse_kegg_gzmk@result,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA KEGG NR vs Control")


enrich_react_1_down <- ReactomePA::enrichPathway(Filter_data_down$entrezgene_id, 
                                                 organism = "human", readable = TRUE, 
                                                 pvalueCutoff = 0.5)

enrich_react_down <- mutate(enrich_react_1_down, 
                            richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_react_down, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("Reactome NR vs Control down regulated genes")


#### gsea

gsea_react_1 <- ReactomePA::gsePathway(geneList2, organism = "human", pvalueCutoff = 0.5, verbose= FALSE, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

Gsea_react_2 <- clusterProfiler::setReadable(x = gsea_react_1, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

clusterProfiler::dotplot(Gsea_react_2, showCategory =10, split=".sign", title = "GSEA reactome NR vs Control") + facet_grid(.~.sign) #plot to show activated and suppressed process


## GSEA plot ##
ewp_gse_react <- arrange(Gsea_react_2, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:10)

ggplot(ewp_gse_react, showCategory=20,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA reactome NR vs Control")


##### Merge csv files to one xls ##

count_data <- rownames_to_column(all_cntrls_1_3)

out <- list("sheet1" = count_data,"sheet2" = cntrl_vs_nsc_merge, 
            "sheet3"= Filter_data, "sheet4"= Filter_data_up,"sheet5"= Filter_data_down,
            "sheet6" = GO_enrich_down_simp_up@result, "sheet7" = GO_enrich_down_simp_down@result, 
            "sheet8" = Gsea_enrich_up_simp@result,"sheet9"= enrich_kegg_up@result, "sheet10"= enrich_kegg_down@result,
            "sheet11" = gse_kegg_cnt_vs_GZMK@result, "sheet12" = enrich_react_1_up@result, "sheet13" = enrich_react_1_down@result, 
            "sheet14" = Gsea_react_2@result)
writexl::write_xlsx(out, "/Users/Tyagi/Desktop/matia_chip/DEGs_1/grp_5_NR_vs_Ctr/NR_vs_cnt_1.2/outfile_NR_vsCntl_ex_1.2.xlsx")


count_data <- rownames_to_column(all_cntrls_1_3)
#z-score
logCPM <- cpm(y_feature_1, log=TRUE) 
z_score <- t(scale(t(logCPM)))
z_score <- rownames_to_column(as.data.frame(z_score))

out <- list("sheet1" = cntrl_vs_nsc_merge, "sheet2"= z_score, 
            "sheet3"= Filter_data, "sheet4"= Filter_data_up,"sheet5"= Filter_data_down,
            "sheet6"= count_data)


writexl::write_xlsx(out, 
                    "/Users/Tyagi/Desktop/pca_matia/pca_after_comments/Analysed_data/R_vs_Cntrl_Expt_2_3_cmn_pt_ex_pr.xlsx")































#### DEseq pca ######

counts_data_nsc_deseq <- as.data.frame(counts_mattia)

colnames(counts_data_nsc_deseq)[1:57] <- RNAseq_nsc_meta_sort$`Sample name`

counts_data_nsc_deseq <- as.matrix(counts_data_nsc_deseq)

Epi_in_deseq <- column_to_rownames(Epi_in, var = "Sample ID")

metadata_car <- cbind(metadata_car, epi_group)

Epi_in_deseq <- Epi_in_deseq[, c(1, 2)]
all(rownames(RNAseq_nsc_meta_sort_2) %in% colnames(counts_data_nsc_deseq))
all(rownames(RNAseq_nsc_meta_sort_2) == colnames(counts_data_nsc_deseq))

dds <- DESeqDataSetFromMatrix(countData = y_feature_1, 
                              colData = metadata_car , 
                              design= ~0+condition)

dds <- estimateSizeFactors(dds)



###### heatmap ####

dds_norm <- rlog(dds)
# Calculate the variance for each gene
variances <- apply(assay(dds_norm), 1, var)

# Determine the upper quartile variance cutoff value
upper_var <- quantile(variances, 0.75)

# Filter the data choosing only genes whose variances are in the upper quartile
df_by_var <- data.frame(assay(dds_norm)) %>%
  dplyr::filter(variances > upper_var)

# Create and store the heatmap object
heatmap <- pheatmap(
  df_by_var,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = FALSE, # Cluster the columns of the heatmap (samples),
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap",
  colorRampPalette(c(
    "green",
    "black",
    "red"
  ))(25
  ),
  scale = "row" # Scale values in the direction of genes (rows)
)




###PCA
vsdata <- vst(dds)


pca_deseq <- plotPCA(vsdata, intgroup= "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_deseq, "percentVar"))


p3 <- ggplot(pca_deseq, aes(PC1, PC2, colour= Patient__ID)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(legend.text = element_text("PatientID"))+
  coord_fixed()

p3

## batch correction
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), design=~0+chemotherapy+replicate)

pca_deseq <- plotPCA(vsdata, intgroup= "chemotherapy", returnData = TRUE)
percentVar <- round(100 * attr(pca_deseq, "percentVar"))


p4 <-ggplot(pca_deseq, aes(PC1, PC2, shape= factor(chemotherapy), colour= factor(replicate))) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



p3 | p4


rld <- rlog(dds)
plotPCA(rld, intgroup= c("Treatment_Type", "Patient_Type"))



####### after feedback from luigi on pca #####
Epi_in <- Epi_in[c(5:14), ]


feature_counts_1 <- feature_counts_1[, c(1, 6:15)]

epi_group <- Epi_in$chemotherapy

epi_group <- c("R_6",  "NR_7", "R_8", "R_9", "NR_10","NR_2","NR_2", "Act_cntrl", "Act_cntrl", "Act_cntrl")

epi_group <- as.factor(epi_group)

colnames(feature_counts_1)[2:11] <- Epi_in$`Sample ID`

Experiment_time <- factor(x = Epi_in$Experiment)
replicate <- factor(x = Epi_in$replicate)
Condition_1 <- factor(x = Epi_in$chemotherapy)


design_new_epi <- model.matrix(~0+Condition_1+Experiment_time, data = feature_counts_1)


#feature_counts_1 <- column_to_rownames(feature_counts_1, var = "Geneid")

rownames(feature_counts_1) <- feature_counts_1[,1]
feature_counts_1 <- feature_counts_1[, -1]

feature_counts_1[is.na(feature_counts_1)] <- 0

y_feature_1 <- DGEList(counts = data.matrix(feature_counts_1))


drop <- filterByExpr(y_feature_1, group = epi_group)
y_feature_1 <- y_feature_1[drop, , keep.lib.sizes=FALSE]


dim(y_feature_1)

### Normalizing the data
y_feature_1 <- calcNormFactors(y_feature_1)


##### MDS plot ####
col.treatment <- c("Red", "Pink", "Orange", "Yellow", "Purple", "Green", "Blue")[epi_group] ### Putting the colors to the treatments ##
data.frame(epi_group, col.treatment)
mds <- plotMDS(cpm(y_feature_1, log = TRUE), col = col.treatment, labels = epi_group)
legend("topleft", cex = 0.5, fill=c("Blue", "green"),legend=levels(epi_group), pch = (Epi_in$replicate))
title("MDS plot Treatment")

## using ggplot
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = epi_group)

library(ggplot2)
p1 <- ggplot(toplot, aes(Dim1, Dim2, colour = epi_group)) +
  geom_point(size=12, stroke = 2, aes(shape=epi_group, color=epi_group, size=epi_group))+ 
  theme(legend.position = "none", text = element_text(size = 20))+
  scale_shape_manual(values = seq(0,10))
p1



### Batch correction
logCPM <- cpm(y_feature_1, log=TRUE, prior.count=0.1)

y_batch <- removeBatchEffect(logCPM, design = design_new_epi)


mds <- plotMDS(y_batch, col = col.treatment, labels = epi_group)


## using ggplot
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = epi_group)
library(ggplot2)
p2 <- ggplot(toplot, aes(Dim1, Dim2, colour = epi_group)) + geom_point(size=12, stroke = 2, aes(shape=epi_group, color=epi_group, size=epi_group))+ 
  theme(text = element_text(size = 20))+scale_shape_manual(values = seq(0,10))

p1 + ggtitle("raw") | p2 +ggtitle("Batchcorrected")

par(mfrow=c(1,2))
boxplot(as.data.frame(y_feature_1), main="raw")
boxplot(as.data.frame(y_batch), main="Corrected")




####### DEGs ############

V1_voom <- voom(y_feature_1, design_new_epi, plot = T)

#corfit<- duplicateCorrelation(V1_voom, design_new_epi, block = Epi_in$)

V1_voom <- voom(y_feature_1, design_new_epi, block = Epi_in$replicate, correlation = corfit$consensus.correlation, plot = T)



fit_1_limma <- lmFit(V1_voom, design_new_epi)


contrast_limma <- makeContrasts(Condition_1NR - Condition_1Act_cntrl, levels = colnames(coef(fit_1_limma))) #### estimating contrast using make contrast function

est_contr <- contrasts.fit(fit_1_limma, contrast_limma) 

est_contr <- eBayes(est_contr) ### smoothing the standard error

result_1 <- decideTests(est_contr)


summary(result_1)




top_table_limma <- topTable(est_contr, sort.by = "P", adjust.method = "fdr", number = Inf)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl" )

top_table_1 <- rownames_to_column(top_table_limma)

colnames(top_table_1)[1] <- "ensembl_gene_id"

cntrl_vs_nsc <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "uniprotswissprot", "entrezgene_id"), filters = "ensembl_gene_id", values = top_table_1, mart = ensembl, uniqueRows = TRUE, )
cntrl_vs_nsc_merge <- merge(top_table_1, cntrl_vs_nsc, by= "ensembl_gene_id", no.dups = TRUE)
cntrl_vs_nsc_merge <- cntrl_vs_nsc_merge[!duplicated(cntrl_vs_nsc_merge$ensembl_gene_id), ] ##### remove duplicate genes


volcano_data <- as.data.frame(cntrl_vs_nsc_merge[, c(2,6)], row.names = cntrl_vs_nsc_merge$external_gene_name)


EnhancedVolcano(toptable = volcano_data, lab = rownames(volcano_data), 
                x = "logFC", y = "adj.P.Val", pCutoff = 0.05, titleLabSize = 20, 
                title = "Non-Responder vs Act_control", FCcutoff = 0.1, labSize = 5) + 
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) + scale_y_continuous(n.breaks = 10)


volcano_data$LogFC <- as.numeric(volcano_data$LogFC)
volcano_data$padj <- as.numeric(volcano_data$padj)


EnhancedVolcano(volcano_data,
                lab = rownames(volcano_data),
                x = 'logFC',
                y = 'adj.P.Val',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.6,
                pointSize = 2.0,
                labSize = 6.0,
                title = "Non-Responder vs Act_control", 
                colAlpha = 1,
                legendPosition = 'left',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, 
                max.overlaps = 100)




###Filtering the data according to the logFC 1.2 and adj.pvalue 0.05 
Filter_data <-dplyr::select(.data = cntrl_vs_nsc_merge, ensembl_gene_id, logFC, adj.P.Val, external_gene_name, entrezgene_id)
Filter_data <- dplyr::filter(.data = Filter_data, logFC > 0 | logFC < -(0), adj.P.Val < 0.05) 

Filter_data_up <- dplyr::filter(.data = Filter_data, logFC > 0, adj.P.Val < 0.05) 
Filter_data_down <- dplyr::filter(.data = Filter_data, logFC < -(0), adj.P.Val < 0.05) 

write.csv(x = metadata_car, file = "/Users/Tyagi/Desktop/metadata_car.csv")


go_uni <- as.character(Filter_data$ensembl_gene_id)


###GO enrichment using cluster profiler

geneList <- Filter_data[, 2] ### logFc

names(geneList) = as.character(Filter_data[, 1]) ### naming logFc

geneList <- sort(geneList, decreasing = TRUE)


### enrichment analysis using ORA
GO_enrich_1_tmrvsFT_up <- clusterProfiler::enrichGO(Filter_data_up$ensembl_gene_id, universe = go_uni, OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL', ont = "ALL", readable = TRUE, pvalueCutoff = 0.5, qvalueCutoff = 0.2, pAdjustMethod = "fdr")

#GO_enrich_down_simp_up <- clusterProfiler::simplify(GO_enrich_1_tmrvsFT_down, cutoff=0.6, by="p.adjust", select_fun= min)


ego_up <- mutate(GO_enrich_1_tmrvsFT_down, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(ego_up, showCategory = 30, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("Biological Processes NR vs Act_control down-regulated genes")



edox2 <- pairwise_termsim(enrich_kegg_cnt_vs_GZMK_down, method = "JC")
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average", showCategory = 30, nCluster =10)
aplot::plot_list(p2, tag_levels='A')







write.csv(x = GO_enrich_down_simp, file = "/Users/Tyagi/Desktop/matia/RvsCntrl/GO_up.csv")



#### gene set enrichment analysis
gse_cnt_vs_GZMK_trt <- clusterProfiler::gseGO(geneList, OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL', ont = "ALL", nPermSimple = 100000, eps = 0, pvalueCutoff = 1, pAdjustMethod = "fdr")

Gsea_enrich_up_simp <- clusterProfiler::simplify(gse_cnt_vs_GZMK_trt, cutoff=0.6, by="p.adjust", select_fun= min)

Gsea_enrich_up_simp <- clusterProfiler::setReadable(x = gse_cnt_vs_GZMK_trt, OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL')

clusterProfiler::dotplot(Gsea_enrich_up_simp, showCategory =20, split=".sign", title = "NR vs Act_control") + facet_grid(.~.sign) #plot to show activated and suppressed process

## GSEA plot ##
ewp_gse_gzmk <- arrange(Gsea_enrich_up_simp, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:5)

ggplot(ewp_gse_gzmk, showCategory=10,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("NR vs Act_control GSEA")



##### Pathway enrichment analysis #####
geneList2 <- Filter_data[, 2] ### logFc
names(geneList2) = as.character(Filter_data[, 5]) ### naming logFc
geneList2 <- sort(geneList2, decreasing = TRUE)


###### Pathway enrichment analysis using enrichkegg ###
enrich_kegg_cnt_vs_GZMK_up <- clusterProfiler::enrichKEGG(Filter_data_up$entrezgene_id, keyType = "kegg", organism = "hsa", pvalueCutoff = 0.5, pAdjustMethod = "fdr", qvalueCutoff = 0.5)

enrich_kegg_up <- clusterProfiler::setReadable(x = enrich_kegg_cnt_vs_GZMK_up, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")


enrich_kegg_up <- mutate(enrich_kegg_up, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_kegg_up, showCategory = 30, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("KEGG pathway NR vs Act_control up regulated genes")




###### Pathway enrichment analysis using gsekegg ###
gse_kegg_cnt_vs_GZMK <- clusterProfiler::gseKEGG(geneList= geneList2, organism = "hsa", pvalueCutoff = 1, verbose= FALSE, keyType = "kegg", nPermSimple = 100000, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

gse_kegg_cnt_vs_GZMK <- clusterProfiler::setReadable(x = gse_kegg_cnt_vs_GZMK, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

## GSEA plot ##
ewp_gse_kegg_gzmk <- arrange(gse_kegg_cnt_vs_GZMK, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  slice(1:10)

ggplot(ewp_gse_kegg_gzmk@result, showCategory=10,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA KEGG NR vs Act control")


enrich_react_1_down <- ReactomePA::enrichPathway(Filter_data_down$entrezgene_id, organism = "human", readable = TRUE, pvalueCutoff = 1)

enrich_react_down <- mutate(enrich_react_1_down, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_react_down, showCategory = 20, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("Reactome NR vs Act control down regulated genes")


#### gsea

gsea_react_1 <- ReactomePA::gsePathway(geneList2, organism = "human", pvalueCutoff = 0.05, verbose= FALSE, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

Gsea_react_2 <- clusterProfiler::setReadable(x = gsea_react_1, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

clusterProfiler::dotplot(Gsea_react_2, showCategory =10, split=".sign", title = "GSEA reactome NR vs Act_control") + facet_grid(.~.sign) #plot to show activated and suppressed process


## GSEA plot ##
ewp_gse_react <- arrange(Gsea_react_2, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:10)

ggplot(ewp_gse_react, showCategory=20,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA reactome NR vs Act_control")


##### Merge csv files to one xls ##

count_data <- rownames_to_column(feature_counts_1)

out <- list("sheet1" = feature_counts_1,"sheet2" = cntrl_vs_nsc_merge, 
            "sheet3"= Filter_data, "sheet4"= Filter_data_up,"sheet5"= Filter_data_down,
            "sheet6" = GO_enrich_down_simp_up@result, "sheet7" = GO_enrich_down_simp_down@result, 
            "sheet8" = Gsea_enrich_up_simp@result,"sheet9"= enrich_kegg_up@result, "sheet10"= enrich_kegg_down@result,
            "sheet11" = gse_kegg_cnt_vs_GZMK@result, "sheet12" = enrich_react_1_up@result, "sheet13" = enrich_react_1_down@result, 
            "sheet14" = Gsea_react_2@result)
writexl::write_xlsx(out, "/Users/Tyagi/Desktop/matia/NR_vs_Act_cntrl/NRvs_act_cntrl_03_11/outfile_NRvsAct_cntl_ex_1.xlsx")






bar_1 <- dplyr::filter(.data = GO_enrich_1_tmrvsFT_up@result, Count > 12)

bar_1 %>%
  ggplot(aes(Count, Description, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ONTOLOGY, nrow = 3, scales = "free") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.y.right= element_text(size = 50),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 20)
        
  ) +
  geom_text(
    aes(label = round(p.adjust, 4)),
    color = "black",
    size = 4,
    hjust=-0.5,
    position = position_dodge(0.01)
    
  )









#####################################################################



######### Carina data ####################

All_LINO_vs_CNTL <- read.delim("C:/Users/Tyagi/Desktop/Carina_lino/All_LINO_vs_CNTL.txt")
L_vs_cntl <- All_LINO_vs_CNTL[, c(1, 3, 6, 7)]
L_vs_cntl <- dplyr::mutate_all(L_vs_cntl, .funs=toupper)

colnames(L_vs_cntl)[1] <- "external_gene_name"

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

L_vs_CN1 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "uniprotswissprot", "entrezgene_id"), filters = "external_gene_name", values = L_vs_cntl, mart = ensembl, uniqueRows = TRUE, )
L_vs_CN1_merge <- merge(L_vs_cntl, L_vs_CN1, by.x= "external_gene_name", no.dups = TRUE)
L_vs_CN1_merge <- L_vs_CN1_merge[!duplicated(L_vs_CN1_merge$external_gene_name), ] ##### remove duplicate genes

colnames(L_vs_cntl)[2] <- "LogFC"


volcano_data <- as.data.frame(L_vs_cntl[, c(2,4)], row.names = L_vs_cntl$external_gene_name)

volcano_data$LogFC <- as.numeric(volcano_data$LogFC)
volcano_data$padj <- as.numeric(volcano_data$padj)

EnhancedVolcano(toptable = volcano_data, lab = rownames(volcano_data), 
                x = "LogFC", y = "padj", pCutoff = 0.05, titleLabSize = 20, 
                title = "Lino vs Control", FCcutoff = 0.6, labSize = 5)+ 
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) + scale_y_continuous(n.breaks = 5)


EnhancedVolcano(volcano_data,
                lab = rownames(volcano_data),
                x = 'LogFC',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.6,
                pointSize = 2.0,
                labSize = 6.0,
                title = "Lino vs Control", 
                colAlpha = 1,
                legendPosition = 'left',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, 
                max.overlaps = 50)





###Filtering the data according to the logFC 1.2 and adj.pvalue 0.05 
Filter_data <-dplyr::select(.data = L_vs_CN1_merge, ensembl_gene_id, log2FoldChange, padj, external_gene_name, entrezgene_id)
Filter_data$log2FoldChange <- as.numeric(Filter_data$log2FoldChange)
Filter_data$padj <- as.numeric(Filter_data$padj)


Filter_data <- dplyr::filter(.data = Filter_data, log2FoldChange > 0.6 | log2FoldChange < -(0.6), padj < 0.05) 



Filter_data_up <- dplyr::filter(.data = Filter_data, log2FoldChange > 0.6, padj < 0.05) 
Filter_data_down <- dplyr::filter(.data = Filter_data, log2FoldChange < -(0.6), padj < 0.05) 

## go_universe
go_sign <- dplyr::filter(Filter_data, padj < 0.05)
go_uni <- as.character(Filter_data$ensembl_gene_id)


bar_1 <- dplyr::filter(.data = GO_enrich_1_tmrvsFT_up@result, Count > 5)

###GO enrichment using cluster profiler

geneList <- Filter_data[, 2] ### logFc
names(geneList) = as.character(Filter_data[, 4]) ### naming logFc

geneList <- sort(geneList, decreasing = TRUE)


### enrichment analysis using ORA
GO_enrich_1_tmrvsFT_up <- clusterProfiler::enrichGO(Filter_data_up$ensembl_gene_id, universe = go_uni,OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL', ont = "ALL", readable = TRUE, pvalueCutoff = 0.5, qvalueCutoff = 0.5, pAdjustMethod = "BH")

GO_enrich_down_simp_up <- clusterProfiler::simplify(GO_enrich_1_tmrvsFT_up, cutoff= 1, by="p.adjust", select_fun= min)


ego_up <- mutate(GO_enrich_1_tmrvsFT_up, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(ego_up, showCategory = 40, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(5, 10)) +
  
  theme_dose(15) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("Biological Processes Lino vs Control up-regulated genes")



edox2 <- pairwise_termsim(GO_enrich_1_tmrvsFT_up, method = "JC")
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average", showCategory = 25, nCluster =10)
aplot::plot_list(p2, tag_levels='A')




#### gene set enrichment analysis
gse_cnt_vs_GZMK_trt <- clusterProfiler::gseGO(geneList, OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL', ont = "ALL", nPermSimple = 100000, eps = 0, pvalueCutoff = 0.05, pAdjustMethod = "fdr")

Gsea_enrich_up_simp <- clusterProfiler::simplify(gse_cnt_vs_GZMK_trt, cutoff=0.6, by="p.adjust", select_fun= min)

Gsea_enrich_up_simp <- clusterProfiler::setReadable(x = Gsea_enrich_up_simp, OrgDb = "org.Hs.eg.db", keyType = 'ENSEMBL')

clusterProfiler::dotplot(Gsea_enrich_up_simp, showCategory =20, split=".sign", title = "GSEA Lino vs Control") + facet_grid(.~.sign) #plot to show activated and suppressed process

## GSEA plot ##
ewp_gse_gzmk <- arrange(Gsea_enrich_up_simp, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  slice(1:5)

ggplot(ewp_gse_gzmk, showCategory=10,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("Lino vs Control GSEA")



##### Pathway enrichment analysis #####
geneList2 <- Filter_data[, 2] ### logFc
names(geneList2) = as.character(Filter_data[, 5]) ### naming logFc
geneList2 <- sort(geneList2, decreasing = TRUE)


go_uni <- as.character(Filter_data$entrezgene_id)


###### Pathway enrichment analysis using enrichkegg ###
enrich_kegg_cnt_vs_GZMK_down <- clusterProfiler::enrichKEGG(Filter_data_down$entrezgene_id, keyType = "kegg", organism = "hsa", pvalueCutoff = 0.5, pAdjustMethod = "fdr", qvalueCutoff = 0.5)

enrich_kegg_down <- clusterProfiler::setReadable(x = enrich_kegg_cnt_vs_GZMK_down, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")


enrich_kegg_down <- mutate(enrich_kegg_down, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_kegg_down, showCategory = 30, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("KEGG pathways")+
  
  ggtitle("KEGG pathway Lino vs Control Down-regulated genes")



edox2 <- pairwise_termsim(enrich_kegg_cnt_vs_GZMK_down, method = "JC")
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average", showCategory = 25, nCluster =10)
aplot::plot_list(p2, tag_levels='A')



###### Pathway enrichment analysis using gsekegg ###
gse_kegg_cnt_vs_GZMK <- clusterProfiler::gseKEGG(geneList= geneList2, organism = "hsa", pvalueCutoff = 0.5, verbose= FALSE, keyType = "kegg", nPermSimple = 100000, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

gse_kegg_cnt_vs_GZMK <- clusterProfiler::setReadable(x = gse_kegg_cnt_vs_GZMK, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

## GSEA plot ##
ewp_gse_kegg_gzmk <- arrange(gse_kegg_cnt_vs_GZMK, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  slice(1:10)

ggplot(ewp_gse_kegg_gzmk@result, showCategory=10,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA KEGG Lino vs Control")


enrich_react_1_up <- ReactomePA::enrichPathway(Filter_data_up$entrezgene_id, organism = "human", readable = TRUE, pvalueCutoff = 0.5)

enrich_react_up <- mutate(enrich_react_1_up, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_react_up, showCategory = 30, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("Reactome terms")+
  
  ggtitle("Reactome Lino vs Control up regulated genes")


#### gsea

gsea_react_1 <- ReactomePA::gsePathway(geneList2, organism = "human", pvalueCutoff = 0.5, verbose= FALSE, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

Gsea_react_2 <- clusterProfiler::setReadable(x = gsea_react_1, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

clusterProfiler::dotplot(Gsea_react_2, showCategory =10, split=".sign", title = "GSEA reactome Lino vs Control") + facet_grid(.~.sign) #plot to show activated and suppressed process


## GSEA plot ##
ewp_gse_react <- arrange(Gsea_react_2, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:5)

ggplot(ewp_gse_react, showCategory=20,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA reactome Lino vs Control")



##### Merge csv files to one xls ##

count_data <- rownames_to_column(feature_counts_1)

out <- list("sheet1" = All_LINO_vs_CNTL,"sheet2" = L_vs_CN1_merge, 
            "sheet3"= Filter_data, "sheet4"= Filter_data_up,"sheet5"= Filter_data_down,
            "sheet6" = GO_enrich_1_tmrvsFT_up@result, "sheet7" = GO_enrich_1_tmrvsFT_down@result, 
            "sheet8" = Gsea_enrich_up_simp@result,"sheet9"= enrich_kegg_up@result, "sheet10"= enrich_kegg_down@result,
            "sheet11" = gse_kegg_cnt_vs_GZMK@result, "sheet12" = enrich_react_1_up@result, "sheet13" = enrich_react_1_down@result, 
            "sheet14" = Gsea_react_2@result)
writexl::write_xlsx(out, "/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx")






#### GO terms ####


df3 <- as.data.frame(GO_enrich_1_tmrvsFT_down@result[, c(1, 3, 10, 7)])


bar_1 <- dplyr::filter(.data = GO_enrich_1_tmrvsFT_up@result, Count > 100)

bar_1 <- dplyr::filter(.data = bar_1, Count < 140)


bar_1 %>%
  ggplot(aes(Count, Description, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ONTOLOGY, nrow = 3, scales = "free") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.y.right= element_text(size = 50),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 20)
        
  ) +
  geom_text(
    aes(label = round(p.adjust, 4)),
    color = "black",
    size = 4,
    hjust=-0.5,
    position = position_dodge(0.01)
    
  )


######Carina visualization ########

### Pathway MAPK ####

outfile_Lino_vs_Cntrl_logfc_0_6 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx", sheet = "Enrich_kegg_up")

# Columns selected from the ORA pathway output with description, geneID
pathway_data_1 <- outfile_Lino_vs_Cntrl_logfc_0_6[15, c(1, 2, 8)]
separated_pathway <- separate_rows(pathway_data_1, "geneID", sep = "/", convert = TRUE)


h_1_heatmap_4 <- merge(x = separated_pathway, y = L_vs_CN1_merge, by.x = "geneID", by.y ="external_gene_name")
h_1_heatmap_4 <- h_1_heatmap_4[, c(10,4)]
h_1_heatmap_4<- column_to_rownames(h_1_heatmap_4, var = "entrezgene_id")
h_1_heatmap_4$log2FoldChange <-as.numeric(h_1_heatmap_4$log2FoldChange) 

p2.out <- pathview(gene.data = h_1_heatmap_4, pathway.id = "hsa04010", species = "hsa", limit = list(gene = max(abs(round(h_1_heatmap_4, 0)))), bins = list(gene= 10), low = "red" , mid = "yellow", high = "blue") 


##### GO terms ####
outfile_Lino_vs_Cntrl_logfc_0_6 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx", 
                                              +     sheet = "GO_up", range = "B281:I281", 
                                              +     col_names = FALSE)


pathway_data_1 <- outfile_Lino_vs_Cntrl_logfc_0_6[, c(1, 2, 3)]
separated_pathway <- separate_rows(pathway_data_1, "...8", sep = "/", convert = TRUE)


h_1_heatmap_4 <- merge(x = separated_pathway, y = L_vs_CN1_merge, by.x = "...8", by.y ="external_gene_name")
h_1_heatmap_4 <- h_1_heatmap_4[, c(10,4)]
h_1_heatmap_4<- column_to_rownames(h_1_heatmap_4, var = "entrezgene_id")



colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[2] <- "LogFC"
colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[3] <- "padj"

volcano_data <- as.data.frame(L_vs_cntl[, c(2,4)], row.names = L_vs_cntl$external_gene_name)

volcano_data$LogFC <- as.numeric(volcano_data$LogFC)
volcano_data$padj <- as.numeric(volcano_data$padj)

EnhancedVolcano(toptable = volcano_data, lab = rownames(volcano_data), 
                x = "LogFC", y = "padj", pCutoff = 0.05, titleLabSize = 20, 
                title = "Lino vs Control", FCcutoff = 0.6, labSize = 5)+ 
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) + scale_y_continuous(n.breaks = 5)

outfile_Lino_vs_Cntrl_logfc_0_6 <- column_to_rownames(outfile_Lino_vs_Cntrl_logfc_0_6, var = "gene")

outfile_Lino_vs_Cntrl_logfc_0_6 <- as.data.frame(outfile_Lino_vs_Cntrl_logfc_0_6)
#Making names uniques for row names#

rownames(outfile_Lino_vs_Cntrl_logfc_0_6) <- make.names(outfile_Lino_vs_Cntrl_logfc_0_6[, 1], unique = TRUE)



# define different cell-types that will be shaded
ER_Ca_homeostasis <- c("Lnpk", "Cisd2", "Tmem147", "Reep4", 
                       "Reep2", "Rhog", "App", "Selenow", "Fcgr2b",
                       "Cebpb", "Bbc3", "Camk2g", "Tmco6", "Psen2", "Selenon", "H6PD")

#Lipid_catabolic_process <- c("ABCA2", "ABCD1", "ABCD4", "ABHD16A", "ABHD8", "ACADS", "ACADVL", "ACOT7", "ACOT8", "ACSBG1", "ACSF3", "AGPAT1", "AKT2", "ALKBH7", "ARMC5", "B4GALNT1", "BSCL2", "CD74", "CPT1A", "CSNK1G2", "CYP2S1", "DAB2IP", "DGAT2", "DGKZ", "DHCR24", "DPM3", "EDF1", "EGR1", "EIF6", "EPHX1", "FASN", "FDXR", "FUT7", "GPAA1", "GPX1", "GPX4", "H6PD", "IMPA2", "IRS2", "ISYNA1", "ITPKA", "KIT", "LDLRAP1", "LGMN", "LIPE", "LPCAT1", "LPCAT4", "LRP10", "LRP5", "LYPLA2", "MBOAT2", "MBOAT7", "MIF", "MTLN", "NR1H2", "NUDT8", "PDGFB", "PDK4", "PER2", "PGP", "PIAS4", "PIGS", "PIGT", "PIK3R2", "PIP5K1C", "PITPNM1", "PITPNM2", "PLBD2", "PLCXD1", "PLPPR3", "PNPLA2", "POR", "PPARD", "RXRA", "SCAP", "SCARB1", "SLC27A1", "SMPD1", "SOCS1", "SORBS1", "SORL1", "SPHK2", "SREBF1", "STUB1", "TGFB1", "TSPO", "TYSND1", "VAC14", "VLDLR", "WDTC1", "WNT10B")

#Survival_Negative_regulator_of_apoptosis_and_T_cell_proliferation <- c("AGAP2", "AKT2.1", "BCL3", "BMP1", "JUNB", "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2", "AAAS")

#Mito_Respirosome_ATP_synthase <- c("ATP5D", "ATP5E", "ATP5G1", "ATP5G2", "ATP5MPL", "ATP6V0A1", "ATPAF1")

#Anabolic_pathway <- c("PKFL", "HK1")

Glycolysis <- c("G6PD", "HK1","PFKL")

Glutaminolysis <- "gsl2"

Lipids_and_PL_metabolism <- c("Srebf1", "Fasn", "DGAT2", "AGPAT1", 
                              "MBOAT2", "MBOAT7", "VLCAD", "REEP2","ABHD5", 
                              "ALDH3B2", "AQP7", "ACSL1", "ACSL5", "ACLY", 
                              "ACSL6", "ACSL4", "ACSL3" )




# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% ER_Ca_homeostasis, 17,
  ifelse(rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% Glycolysis, 11,
         ifelse (rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% Glutaminolysis, 64,
                 ifelse (rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% Lipids_and_PL_metabolism, 2, 3))))

keyvals.shape[is.na(keyvals.shape)] <- 3
names(keyvals.shape)[keyvals.shape == 17] <- "ER_Ca_homeostasis"
names(keyvals.shape)[keyvals.shape == 11] <- "Glycolysis"
names(keyvals.shape)[keyvals.shape == 64] <- "Glutaminolysis"
names(keyvals.shape)[keyvals.shape == 2] <- "Lipids_and_PL_metabolism"



## creating color key 
keyvals.colour <- ifelse(
  outfile_Lino_vs_Cntrl_logfc_0_6$LogFC < -0.2, 'dark gray',
  ifelse(outfile_Lino_vs_Cntrl_logfc_0_6$LogFC > 0.2, 'dark gray',
         'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low'


p1 <- EnhancedVolcano(outfile_Lino_vs_Cntrl_logfc_0_6,
                      lab = rownames(outfile_Lino_vs_Cntrl_logfc_0_6),
                      x = "LogFC", y = "padj",
                      selectLab = rownames(outfile_Lino_vs_Cntrl_logfc_0_6)
                      [which(names(keyvals.shape) %in% c("ER_Ca_homeostasis", 
                                                         "Lipids_and_PL_metabolism", 
                                                         "Glutaminolysis", 
                                                         "Glycolysis"))],
                      xlab = bquote(~Log[2]~ 'fold change'),
                      ylab = bquote(""-Log[10]*adjPvalue),
                      title = NULL,
                      ylim = c(0,10),
                      xlim = c(-1, 2),
                      pCutoff = 0.05,
                      FCcutoff = 0.1,
                      pointSize = 3,
                      labSize = 4,
                      titleLabSize = NULL,
                      subtitleLabSize = NULL,
                      captionLabSize = NULL,
                      subtitle = NULL,
                      shapeCustom = keyvals.shape,
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      legendLabSize = 8,
                      legendPosition = NULL,
                      legendIconSize = 3.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.3,
                      colConnectors = 'grey50',
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'partial',
                      borderWidth = 1.2,
                      borderColour = 'black')

p1



p1[["data"]]$lab[120] <- "PKFL"






#### heatmap ###
heatmap_f3 <- lino_heatmap[, c(1:6)]


heatmap_f5 <- column_to_rownames(.data = heatmap_f3,  var = "Gene_name")
heatmap_22 <- as.matrix(heatmap_f5)
rlm <- rainbow(ncol(as.matrix(heatmap_csv)))
heatmap.2(data[1:20, ], 
          keysize = 0.9, 
          trace = "none", 
          key = TRUE,
          col = bluered(255), 
          dendrogram="none",  
          srtCol = 0, 
          adjCol = c(0.5, 0), 
          margins = c(2,20),
          offsetRow = 0,
          offsetCol = 0,
          cexRow = 1.2
            
          
              ) 
legend("left", legend = c("Sphingolipid signaling pathway", "Glycerophospholipid metabolism", "Phospholipase D signaling pathway", "Sphingolipid metabolism", "Mitochondrial Fatty Acid Beta-Oxidation")
       , col = rlm, lty= 1,lwd = 8,bty = "n") 



heatmap_f5 <- column_to_rownames(.data = Heatmap_data_1,  var = "gene")
heatmap_22 <- as.matrix(heatmap_f5)
rlm <- rainbow(ncol(heatmap_22))
heatmap.2(heatmap_22,  main = "", 
          keysize = 0.6, key.title = NA, key.xlab = NA, 
          key.ylab = NA, dendrogram = "column", col = redblue(255), Rowv = FALSE,
          Colv = TRUE, key = TRUE, trace = "none",  margins = c(3,10), 
          cexCol = NULL, cexRow = NULL, 
          srtCol = 0, adjCol = c(0.5, 0), offsetRow = 0) 
legend("bottom", legend = c("Treatment-1","Treatment-3", "Treatment-4"), col = rlm, lty= 1,lwd = 8,bty = "n") 




heatmap_f5 <- column_to_rownames(.data = heatmap_f3,  var = "external_gene_name")
heatmap_22 <- as.matrix(heatmap_f5)
rlm <- rainbow(ncol(heatmap_22))
heatmap.2(heatmap_22, ColSideColors = rlm, main = input$heat_text, 
          keysize = input$heatmap_input3, key.title = NA, 
          key.xlab = NA, key.ylab = NA, dendrogram = "none", 
          col = bluered(255), Rowv = FALSE, Colv = FALSE, key = TRUE, 
          trace = "none",  margins = c(3,22), cexCol = input$heatmap_input1, 
          cexRow = input$heatmap_input2, srtCol = 0, adjCol = c(0.5, 0), offsetRow = 0,) 
legend("topright", legend = c("Treatment-1","Treatment-3", "Treatment-4"), col = rlm, lty= 1,lwd = 8,bty = "n") 



#### gmt cp #####
cp2 <- clusterProfiler::read.gmt("C:/Users/Tyagi/Desktop/gmt_file/c2.cp.v2022.1.Hs.symbols.gmt")

cp2_gmt <- enrichR::enrichr(geneList, TERM2GENE = cp2, verbose = FALSE, pvalueCutoff = 0.05, nPermSimple = 100000, eps = 0, pAdjustMethod = "fdr", by = "fgsea")

## GSEA plot ##
ewp_gse_kegg_gzmk <- arrange(cp2_gmt, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:10)

ggplot(ewp_gse_kegg_gzmk@result, showCategory=20,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA C2 Lino vs Control")


clusterProfiler::dotplot(cp8_immuno_gmt, showCategory =10, split=".sign", title = "GSEA C8 Lino vs Control") + facet_grid(.~.sign) #plot to show activated and suppressed process

out <- list("sheet1" = cp2_gmt@result, "sheet2" = cp3_TFT_gmt@result, 
            "sheet3"= cp7_immuno_gmt@result, "sheet4"= cp8_immuno_gmt@result)
writexl::write_xlsx(out, "/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_gsea_MsigDB.xlsx")






### comments ## 

outfile_Lino_vs_Cntrl_logfc_0_6 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx")

outfile_Lino_vs_Cntrl_logfc_0_6 <- outfile_Lino_vs_Cntrl_logfc_0_6[, c(1, 2 ,7)]

outfile_Lino_vs_Cntrl_logfc_0_6 <- lapply(outfile_Lino_vs_Cntrl_logfc_0_6, toupper)

outfile_Lino_vs_Cntrl_logfc_0_6 <- as.data.frame(outfile_Lino_vs_Cntrl_logfc_0_6)

colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[2] <- "LogFC"
colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[3] <- "p_value"


outfile_Lino_vs_Cntrl_logfc_0_6 <- column_to_rownames(outfile_Lino_vs_Cntrl_logfc_0_6, var = "gene")


outfile_Lino_vs_Cntrl_logfc_0_6$LogFC <- as.numeric(outfile_Lino_vs_Cntrl_logfc_0_6$LogFC)
outfile_Lino_vs_Cntrl_logfc_0_6$p_value <- as.numeric(outfile_Lino_vs_Cntrl_logfc_0_6$p_value)

outfile_Lino_vs_Cntrl_logfc_0_6$p_value <- sort(outfile_Lino_vs_Cntrl_logfc_0_6$p_value, decreasing = TRUE)

# define different cell-types that will be shaded
ER_Ca_homeostasis <- c("Lnpk", "Cisd2", "Tmem147", "Reep4", "H6PD",
                       "Reep2", "Rhog", "App", "Selenow", "Fcgr2b",
                       "Cebpb", "Bbc3", "Camk2g", "Tmco6", "Psen2", "Selenon")

ER_Ca_homeostasis <- toupper(ER_Ca_homeostasis)
Glycolysis <- c("G6PDX", "HK1","PFKL")

Glutaminolysis <- "GLS2"

Lipids_and_PL_metabolism <- c("SREBF1", "FASN", "DGAT2", "AGPAT1", "ACCS", "SLC25A1", "DGAT2",
                              "MBOAT2", "MBOAT7", "VLCAD", "REEP2","ABHD5", "BSCL2","HILPDA","PNPLA2",
                              "ALDH3B2", "AQP7", "ACSL1", "ACSL5", "ACLY", "PLN2","LIPE", "SLC27A1", "ABCA2",
                              "ACSL6", "ACSL4", "ACSL3","CPT1A")




# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.colour <- ifelse(
  rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% Glycolysis & outfile_Lino_vs_Cntrl_logfc_0_6$LogFC > 0.1 & outfile_Lino_vs_Cntrl_logfc_0_6$p_value < 0.05, 'blue',
  ifelse(rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% Glutaminolysis & outfile_Lino_vs_Cntrl_logfc_0_6$LogFC > 0.1 & outfile_Lino_vs_Cntrl_logfc_0_6$p_value < 0.05, 'black',
         ifelse(rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% Lipids_and_PL_metabolism & outfile_Lino_vs_Cntrl_logfc_0_6$LogFC > 0.1 & outfile_Lino_vs_Cntrl_logfc_0_6$p_value < 0.05, 'green',
                ifelse(rownames(outfile_Lino_vs_Cntrl_logfc_0_6) %in% ER_Ca_homeostasis & outfile_Lino_vs_Cntrl_logfc_0_6$LogFC > 0.1 & outfile_Lino_vs_Cntrl_logfc_0_6$p_value < 0.05, 'red',
                       ifelse(outfile_Lino_vs_Cntrl_logfc_0_6$p_value < 0.05, "gray74", 'gray89')))))


names(keyvals.colour)[keyvals.colour == "gray89"] <- "Non-significant"
names(keyvals.colour)[keyvals.colour == "red"] <- "ER_Ca_homeostasis"
names(keyvals.colour)[keyvals.colour == "blue"] <- "Glycolysis"
names(keyvals.colour)[keyvals.colour == "black"] <- "Glutaminolysis"
names(keyvals.colour)[keyvals.colour == "green"] <- "Lipids_and_PL_metabolism"
names(keyvals.colour)[keyvals.colour == 'gray74'] <- 'Up/Down'


## creating color key 
keyvals.colour <- ifelse(
  outfile_Lino_vs_Cntrl_logfc_0_6$LogFC < -0.2, 'dark gray',
  ifelse(outfile_Lino_vs_Cntrl_logfc_0_6$LogFC > 0.2, 'dark gray',
         'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low'


p1 <- EnhancedVolcano(outfile_Lino_vs_Cntrl_logfc_0_6,
                      lab = rownames(outfile_Lino_vs_Cntrl_logfc_0_6),
                      x = "LogFC", y = "p_value",
                      selectLab = rownames(outfile_Lino_vs_Cntrl_logfc_0_6)
                      [which(names(keyvals.colour) %in% c("ER_Ca_homeostasis", 
                                                          "Lipids_and_PL_metabolism", 
                                                          "Glutaminolysis", 
                                                          "Glycolysis"))],
                      xlab = bquote(~Log[2]~ 'fold change'),
                      title = NULL,
                      ylim = c(0,15),
                      xlim = c(-2, 2),
                      pCutoff = 0.05,
                      FCcutoff = 0.1,
                      pointSize = 1,
                      labSize = 4,
                      titleLabSize = NULL,
                      subtitleLabSize = NULL,
                      captionLabSize = NULL,
                      subtitle = NULL,
                      shapeCustom = NULL,
                      colCustom = keyvals.colour,
                      colAlpha = 10,
                      labFace = "bold",
                      legendLabSize = 11,
                      legendPosition = "right",
                      legendIconSize = 7.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      colConnectors = 'grey50',
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'partial',
                      borderWidth = 1.2,
                      borderColour = 'black',
                      max.overlaps = Inf)

p1

ggsave("volcano.eps", width = 30, height = 20, units = "cm")

### GO term ###


outfile_Lino_vs_Cntrl_logfc_0_6 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx", 
                                              +     sheet = "GO_up")



go_terms <- c("GO:0019216", "GO:0006869", "GO:0044242", "GO:0034440", "GO:0071396", "GO:0016042", "GO:0030258", "GO:0006629", "GO:0019395")


outfile_Lino_vs_Cntrl_logfc_0_6 %>% select(all_of(go_term))

go_terms <- inner_join(outfile_Lino_vs_Cntrl_logfc_0_6, go_term, by= "ID", no.dups = TRUE)



reactome <- c("R-HSA-400206", "R-HSA-1369062", "R-HSA-1483206", 
              "R-HSA-1483257", "R-HSA-390918", "R-HSA-1660662", 
              "R-HSA-75109", "R-HSA-8979227", "R-HSA-163560", "R-HSA-428157",
              "R-HSA-897886", "R-HSA-77289", "R-HSA-75105", "R-HSA-75876", "R-HSA-390247",
              "R-HSA-77286", "R-HSA-1655829", "R-HSA-191273")

reactome <- as.data.frame(reactome)

colnames(reactome) <- "ID"
go_terms_ret <- inner_join(outfile_Lino_vs_Cntrl_logfc_0_6, reactome, by= "ID", no.dups = TRUE)

final_go <- bind_rows(go_terms, go_terms_ret)



###
final_go_plot <- read_excel("final_go_plot.xlsx", sheet = "Sheet3")


final_go_plot <- as.data.frame(final_go_plot)
outfile_Lino_vs_Cntrl_logfc_0_6 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx")

outfile_Lino_vs_Cntrl_logfc_0_6 <- outfile_Lino_vs_Cntrl_logfc_0_6[, c(1, 2 ,3)]
colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[1] <- "ID"
colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[2] <- "logFC"
colnames(outfile_Lino_vs_Cntrl_logfc_0_6)[3] <- "adj.P.Val"
library(GOplot)

circ <-circle_dat(final_go_plot, outfile_Lino_vs_Cntrl_logfc_0_6)



plot_id_1 <- c("GO:0010876","GO:0071396", "GO:0016042", "GO:0043066", 
               "GO:0006629", "GO:0005743", "GO:0005759", 
               "GO:0031227", "GO:0006874", "R-HSA-400206", 
               "R-HSA-1483257", "R-HSA-77289", "R-HSA-75105",  
               "R-HSA-428157", "hsa04660", "hsa04020")
GO_data <- GOCircle(data = circ, nsub = plot_id_1, label.size = 4)

setEPS()                                             # Set postscript arguments
postscript("Go_circle.eps", height = 13, width = 17)                           # Start graphics device driver
plot(GO_data )                                 # Create plot
dev.off()                                            # Finish export



#### Z-score and heatmap ####

final_go_plot <- read_excel("final_go_plot.xlsx", sheet = "Sheet4")
final_go_plot <- as.data.frame(final_go_plot)

final_go_plot[is.na(final_go_plot)] <- 0
final_go_plot_1 <- column_to_rownames(final_go_plot, var = "gene")
logCPM <- cpm(final_go_plot_1, log = T)
z_score <- t(scale(t(logCPM))) #row normalization

z_score_table <- rownames_to_column(z_score)
z_score_table <- lapply(z_score_table, toupper)

glutamino_table <- z_score_table[z_score_table$rowname  %in% c("GLS", "GLS2", "GLUD1", "Got1", "Got2", "Gpt"),]


ABC_cycle <- c("Tap1", "Abcb6", "Abcb11", "Abcc1", "Abcg1", "Abcg4", "Abcb6", 
               "Abca1", "Abcb6", "Abcc1", "Abcc5", "Abcg4", "Abca3", "Abca8a", 
               "Abca8b", "Abcb1a", "Abcb1b", "Abcb6", "Abcb7", "Abcb9", "Abcb11",
               "Abcc1", "Abcc2", "Abcc3", "Abcc4", "Abcc5", "Abcc6", "Abcc10", "Abcd1", 
               "Abcd4", "Abcg1", "Abcg2", "Abcg3", "Tap1", "Tap2")
ABC_cycle <- toupper(ABC_cycle)


Glycolysis_path <- c("ALDOA", "BPGM", "ENO1", "ENO2", "GAPDH", "GPI", "HK1", "HK2", "HKDC1", "PFKL", "PFKM", "PGAM1", "PGAM2", "PGAM4", "PGK1", "PKLR", "PKM", "TPI1")


Sphingolipids_table <- z_score_table[z_score_table$rowname  %in% c("AKT2", "GNAI2", "MAP2K2", "MAPK12", "PIK3R2", "RAC2", "RELA", "S1PR4", "SMPD1", "SPHK2"), ]

Glycerophospholipid <- z_score_table[z_score_table$rowname  %in% c("LPCAT1", "LPCAT4", "LYPLA2", "MBOAT2", "MBOAT7", "PLD3", "AGPAT1", "DGKZ"), ]

Phospholipase <- z_score_table[z_score_table$rowname  %in% c("AKT2", "SPHK2", "AGPAT1", "ADCY3", "PIK3R2", "PIP5K1C", "RALGDS", "DGKZ", "MAP2K2", "PDGFB", "KIT"), ]


Metabolism <- z_score_table[z_score_table$rowname %in% c(MTB), ]


TCR_signaling <- z_score_table[z_score_table$rowname %in% c("ZAP70", "MAPK12", "RELA", "NFKBIA", 
                                                            "MAP3K14", "PAK6", "PIK3R2", "JUN", 
                                                            "JUNB", "JUND"), ]
Survival <-  z_score_table[z_score_table$rowname %in% c("AGAP2", "AKT2", "BCL3", "BMP1", "JUNB", 
                                                        "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2"), ]


SRV <- c("AGAP2", "AKT2", "BCL3", "BMP1", "JUNB", 
         "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2")

MTB <- c("HK1", "HK2", "PFKM", "ALDOA", 
         "ACO2", "OGDH", "TPI1", "H6PD", 
         "PGLS", "TKT", "TALDO1", "GLS2", 
         "CAD", "DHODH", "SLC27A1", "AACS", 
         "ACAA2", "FASN", "SREBF1", "CPT1A",
         "ACADVL", "SLC25A1", "MGAT1", "AGPAT1",
         "AGPAT2", "AGPAT3", "LPCAT1", "MBOAT2",
         "MBOAT7", "BSCL2", "HILPDA", "PNPLA2", 
         "PLIN2", "LIPE", "SELENOW", "Cebpb", 
         "Bbc3", "ATF6")

MTB <- toupper(MTB)

TCR_1 <- c("ZAP70", "MAPK12", "RELA", "NFKBIA", 
           "MAP3K14", "PAK6", "PIK3R2", "JUN", 
           "JUNB", "JUND")



### heatmap generation ##
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}


heatmap_f5 <- unnest(Metabolism, cols = rowname) 
heatmap_f5 <- as.data.frame(heatmap_f5)

idx <- sapply(MTB, function(x) {which(Metabolism$rowname == x)})

Metabolism <- Metabolism[idx,]
rownames(Metabolism) <- NULL


heatmap_f5 <- column_to_rownames(.data = Metabolism, var = "rowname")


heatmap_f5$lino_1 <- as.numeric(heatmap_f5$lino_1)
heatmap_f5$lino_2 <- as.numeric(heatmap_f5$lino_2)
heatmap_f5$lino_3 <- as.numeric(heatmap_f5$lino_3)
heatmap_f5$ova_1 <- as.numeric(heatmap_f5$ova_1)
heatmap_f5$ova_2 <- as.numeric(heatmap_f5$ova_2)
heatmap_f5$ova_3 <- as.numeric(heatmap_f5$ova_3)



heatmap_22 <- as.matrix((heatmap_f5))


set.seed(123)
tf<-sample(rownames(heatmap_22))

heatmap.2(heatmap_22, ColSideColors =, main = "",colsep = 3, 
          labRow = make_bold_names(heatmap_22,rownames, tf),
          keysize = 1, key.title = "z-score", key.xlab = NA, dendrogram = "none",
          key.ylab = NA, col = bluered(255), Rowv = FALSE,
          Colv = FALSE, key = TRUE, trace = "none",  margins = c(3,20), 
          cexCol = NULL, cexRow = NULL, 
          srtCol = 0, adjCol = c(0.5, 0), offsetRow = 0)

#724, 480



#### z-score try 2 ####
library(tidyverse)
library(edgeR)
library(tidyverse)
library(dplyr)

outfile_Lino_vs_Cntrl_logfc_0_6_2 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx")
outfile_Lino_vs_Cntrl_logfc_0_6_2 <- outfile_Lino_vs_Cntrl_logfc_0_6_2[, c(1, 8,9,10,11,12,13)]

final_go_plot[is.na(final_go_plot)] <- 0
outfile_Lino_vs_Cntrl_logfc_0_6_2 <- column_to_rownames(outfile_Lino_vs_Cntrl_logfc_0_6_2, var = "gene")
logCPM_2 <- outfile_Lino_vs_Cntrl_logfc_0_6_2


logCPM <- rownames_to_column(as.data.frame(logCPM))
logCPM_2 <- rownames_to_column(as.data.frame(logCPM_2))

logCPM_2 <- column_to_rownames(as.data.frame(logCPM_2), var = "rowname") 


z_score_3 <- sapply(logCPM_2, function(logCPM_2) (logCPM_2-mean(logCPM_2))/sd(logCPM_2))

row.names(z_score_3) <- rownames(logCPM_2)
z_score_3 <- data.frame(z_score_3)

z_score_table_3 <- rownames_to_column(z_score_3)
z_score_table_3 <- lapply(z_score_table_3, toupper)

z_score_table_3 <- data.frame(z_score_table_3)
Metabolism <- z_score_table_3[z_score_table_3$rowname %in% c(MTB), ]


TCR_signaling <- z_score_table_3[z_score_table_3$rowname %in% c("ZAP70", "MAPK12", "RELA", "NFKBIA", 
                                                                "MAP3K14", "PAK6", "PIK3R2", "JUN", 
                                                                "JUNB", "JUND"), ]
Survival <-  z_score_table[z_score_table$rowname %in% c("AGAP2", "AKT2", "BCL3", "BMP1", "JUNB", 
                                                        "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2"), ]


SRV <- c("AGAP2", "AKT2", "BCL3", "BMP1", "JUNB", 
         "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2")

MTB <- c("HK1", "HK2", "PFKM", "ALDOA", 
         "ACO2", "OGDH", "TPI1", "H6PD", 
         "PGLS", "TKT", "TALDO1", "GLS2", 
         "CAD", "DHODH", "SLC27A1", "AACS", 
         "ACAA2", "FASN", "SREBF1", "CPT1A",
         "ACADVL", "SLC25A1", "MGAT1", "AGPAT1",
         "AGPAT2", "AGPAT3", "LPCAT1", "MBOAT2",
         "MBOAT7", "BSCL2", "HILPDA", "PNPLA2", 
         "PLIN2", "LIPE", "SELENOW", "Cebpb", 
         "Bbc3", "ATF6")

MTB <- toupper(MTB)

TCR_1 <- c("ZAP70", "MAPK12", "RELA", "NFKBIA", 
           "MAP3K14", "PAK6", "PIK3R2", "JUN", 
           "JUNB", "JUND")


heatmap_f5 <- unnest(Metabolism, cols = rowname) 
heatmap_f5 <- as.data.frame(heatmap_f5)

idx <- sapply(MTB, function(x) {which(Metabolism$rowname == x)})

Metabolism <- Metabolism[idx,]
rownames(Metabolism) <- NULL


heatmap_f5 <- column_to_rownames(.data = Metabolism, var = "rowname")


colnames(heatmap_f5) <- c("lino_1", "lino_2", "lino_3", "Control_1", "Control_2", "Control_3")


heatmap_f5$lino_1 <- as.numeric(heatmap_f5$lino_1)
heatmap_f5$lino_2 <- as.numeric(heatmap_f5$lino_2)
heatmap_f5$lino_3 <- as.numeric(heatmap_f5$lino_3)
heatmap_f5$Control_1 <- as.numeric(heatmap_f5$Control_1)
heatmap_f5$Control_2 <- as.numeric(heatmap_f5$Control_2)
heatmap_f5$Control_3 <- as.numeric(heatmap_f5$Control_3)



heatmap_22 <- as.matrix((heatmap_f5))


set.seed(123)
tf<-sample(rownames(heatmap_22))

heatmap.2(heatmap_22, ColSideColors =, main = "",colsep = 3, scale = "row",
          labRow = make_bold_names(heatmap_22,rownames, tf),
          keysize = 1, key.title = "z-score", key.xlab = NA, dendrogram = "none",
          key.ylab = NA, col = bluered(255), Rowv = FALSE,
          Colv = FALSE, key = TRUE, trace = "none",  margins = c(3,20), 
          cexCol = NULL, cexRow = NULL, 
          srtCol = 0, adjCol = c(0.5, 0), offsetRow = 0)

## 




#### try_3 with only z score applying ###

outfile_Lino_vs_Cntrl_logfc_0_6_2 <- read_excel("C:/Users/Tyagi/Desktop/Carina_lino/Lino_vs_Cntrl/Lino_lgfc_0.6/outfile_Lino_vs_Cntrl_logfc_0.6.xlsx")
outfile_Lino_vs_Cntrl_logfc_0_6_2 <- outfile_Lino_vs_Cntrl_logfc_0_6_2[, c(1, 8,9,10,11,12,13)]


outfile_Lino_vs_Cntrl_logfc_0_6_2 <- column_to_rownames(outfile_Lino_vs_Cntrl_logfc_0_6_2, var = "gene")

z_score_table_3_1 <- rownames_to_column(outfile_Lino_vs_Cntrl_logfc_0_6_2)
z_score_table_3_1 <- lapply(z_score_table_3_1, toupper)

z_score_table_3_1 <- data.frame(z_score_table_3_1)


colnames(z_score_table_3_1) <- c("rowname", "lino_1", "lino_2", "lino_3", "Control_1", "Control_2", "Control_3")


z_score_table_3_1 <- column_to_rownames(z_score_table_3_1, var = "rowname")


z_score_table_3_1$lino_1 <- as.numeric(z_score_table_3_1$lino_1)
z_score_table_3_1$lino_2 <- as.numeric(z_score_table_3_1$lino_2)
z_score_table_3_1$lino_3 <- as.numeric(z_score_table_3_1$lino_3)
z_score_table_3_1$Control_1 <- as.numeric(z_score_table_3_1$Control_1)
z_score_table_3_1$Control_2 <- as.numeric(z_score_table_3_1$Control_2)
z_score_table_3_1$Control_3 <- as.numeric(z_score_table_3_1$Control_3)


z_score_table_3_1 <- column_to_rownames(z_score_table_3_1, var = "rowname")


z_score_table_3_1_1 <- as.matrix(z_score_table_3_1)
z_score_table_3_1_2 <- t(scale(t(z_score_table_3_1)))

z_score_table_3_1_2 <- data.frame(z_score_table_3_1_2)
z_score_table_3_1_2 <- rownames_to_column(z_score_table_3_1_2)



Metabolism <- z_score_table_3_1_2[z_score_table_3_1_2$rowname %in% c(MTB), ]


TCR_signaling <- z_score_table_3_1_2[z_score_table_3_1_2$rowname %in% c(TCR_2), ]
Survival <-  z_score_table_3_1_2[z_score_table_3_1_2$rowname %in% c("AGAP2", "AKT2", "BCL3", "BMP1", "JUNB", 
                                                                    "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2"), ]


SRV <- c("AGAP2", "AKT2", "BCL3", "BMP1", "JUNB", 
         "JUND", "BCL2L12", "IL12RB1", "JAK3", "NCK2")

MTB <- c("HK1", "HK2", "PFKM", "ALDOA", 
         "ACO2", "OGDH", "TPI1", "H6PD", 
         "PGLS", "TKT", "TALDO1", "GLS2", 
         "CAD", "DHODH", "SLC27A1", "AACS", 
         "ACAA2", "FASN", "SREBF1", "CPT1A",
         "ACADVL", "SLC25A1", "MGAT1", "AGPAT1",
         "AGPAT2", "AGPAT3", "LPCAT1", "MBOAT2",
         "MBOAT7", "BSCL2", "HILPDA", "PNPLA2", 
         "PLIN2", "LIPE", "SELENOW", "Cebpb", 
         "Bbc3", "ATF6")

MTB <- toupper(MTB)

TCR_1 <- c("ZAP70", "MAPK12", "RELA", "NFKBIA", 
           "MAP3K14", "PAK6", "PIK3R2", "JUN", 
           "JUNB", "JUND")

TCR_2 <- c("PGLS", "TKT", "AACS", "ACAA2", "FASN", "SREBF1", 
           "CPT1A", "MGAT1", "MBOAT7", "HILPDA", "JUNB", "JUND", "PAK4", "AKT2", "BCL2L12", "BMP1")


heatmap_f5 <- unnest(Metabolism, cols = rowname) 
heatmap_f5 <- as.data.frame(heatmap_f5)

idx <- sapply(SRV, function(x) {which(Survival$rowname == x)})

Survival <- Survival[idx,]
rownames(Survival) <- NULL


heatmap_f5 <- column_to_rownames(.data = Survival, var = "rowname")



heatmap_22 <- as.matrix((heatmap_f5))


set.seed(123)
tf<-sample(rownames(heatmap_22))

heatmap.2(heatmap_22, ColSideColors =, main = "",colsep = 3, scale = "none",
          labRow = make_bold_names(heatmap_22,rownames, tf),
          keysize = 1, key.title = "z-score", key.xlab = NA, dendrogram = "none",
          key.ylab = NA, col = bluered(255), Rowv = FALSE,
          Colv = FALSE, key = TRUE, trace = "none",  margins = c(3,20), 
          cexCol = NULL, cexRow = NULL, 
          srtCol = 0, adjCol = c(0.5, 0), offsetRow = 0)
### Metabolism t3 image by this example, scale none
### TCR signaling t3 image by this example, scale none (900, 800)
### Survival t3 image by this example, scale none (900, 800)


