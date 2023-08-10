
#### overlapiing genes NR_vs_control_vs_R_Expt_2_3_cmn_pt
  
NR_vs_Contrl_cmn_pt_2_3_updated <- read_excel("C:/Users/Tyagi/Desktop/pca_matia/pca_after_comments/Analysed_data/Live_vs_Cntrl_Expt_2_3_cmn_pt_ex_pr.xlsx" , 
                                              sheet = "DEGs_obtained")

NR_vs_Contrl_cmn_pt_2_3_updated <- as.data.frame(NR_vs_Contrl_cmn_pt_2_3_updated)

NR_vs_Contrl_cmn_pt_2_3_updated <- dplyr::filter(.data = NR_vs_Contrl_cmn_pt_2_3_updated, 
                                                 logFC > 0 | logFC < -(0), P.Value < 0.01) 


NR_vs_Control_vs_R_Experiment_2_3_cmn_pt <- read.csv("C:/Users/Tyagi/Desktop/pca_matia/pca_after_comments/Analysed_data/overlap/Live_vs_Control_vs_Sup_Experiment_2_and_3_common_patinet.csv")
NR_vs_Control_vs_R_Experiment_2_3_cmn_pt <- as.data.frame(NR_vs_Control_vs_R_Experiment_2_3_cmn_pt)
NR_vs_Control_vs_R_Experiment_2_3_cmn_pt[is.na(NR_vs_Control_vs_R_Experiment_2_3_cmn_pt)] <- "" 

colnames(NR_vs_Control_vs_R_Experiment_2_3_cmn_pt)[3] <- "external_gene_name"

NR_vs_Control_vs_R_Experiment_2_3_cmn_pt <- na.omit(NR_vs_Control_vs_R_Experiment_2_3_cmn_pt)

unique_live_2960 <- semi_join(NR_vs_Contrl_cmn_pt_2_3_updated, 
                      NR_vs_Control_vs_R_Experiment_2_3_cmn_pt, 
                      by="external_gene_name")


cntrl_vs_nsc_merge <- as.data.frame(Uniq_live_2960)

volcano_data <- as.data.frame(cntrl_vs_nsc_merge[, c(2,6)], 
                              row.names= cntrl_vs_nsc_merge$external_gene_name)



EnhancedVolcano(toptable = volcano_data, lab = rownames(volcano_data), 
                x = "logFC", y = "adj.P.Val", pCutoff = 0.05, titleLabSize = 20, 
                title = "Unique NR vs Control experiment 1 and 3", FCcutoff = 0, labSize = 5) + 
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) + scale_y_continuous(n.breaks = 5)



###Filtering the data according to the logFC 1.2 and adj.pvalue 0.05 
Filter_data <-dplyr::select(.data = cntrl_vs_nsc_merge, ensembl_gene_id, logFC, P.Value, adj.P.Val, external_gene_name, entrezgene_id)
Filter_data <- dplyr::filter(.data = Filter_data, logFC > 0 | logFC < -(0), P.Value< 0.01) 

Filter_data_up <- dplyr::filter(.data = Filter_data, logFC > 0,P.Value< 0.01) 
Filter_data_down <- dplyr::filter(.data = Filter_data, logFC < -(0), P.Value< 0.01) 

write.csv(x = metadata_car, file = "/Users/Tyagi/Desktop/metadata_car.csv")


go_uni <- as.character(Filter_data$ensembl_gene_id)


###GO enrichment using cluster profiler

geneList <- Filter_data[, 2] ### logFc

names(geneList) = as.character(Filter_data[, 1]) ### naming logFc

geneList <- sort(geneList, decreasing = TRUE)


### enrichment analysis using ORA
GO_enrich_1_tmrvsFT_down <- clusterProfiler::enrichGO(Filter_data_down$ensembl_gene_id,
                                                    OrgDb = "org.Hs.eg.db",
                                                    keyType = 'ENSEMBL', ont = "ALL",
                                                    readable = TRUE, pvalueCutoff = 0.5, 
                                                    qvalueCutoff = 0.2, pAdjustMethod = "fdr")

GO_enrich_down_simp_down <- clusterProfiler::simplify(GO_enrich_1_tmrvsFT_down,
                                                  cutoff=0.5, by="p.adjust", 
                                                  select_fun=min,
                                                  measure="Wang")


ego_down <- dplyr::mutate(GO_enrich_down_simp_down, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(ego_down, showCategory = 10, aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("GO terms")+
  
  ggtitle("Biological Processes Live  down-regulated genes")



edox2 <- pairwise_termsim(GO_enrich_1_tmrvsFT_down, method = "JC")
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average", showCategory = 10, nCluster=4)
aplot::plot_list(p2, tag_levels='A')



#### gene set enrichment analysis
gse_cnt_vs_GZMK_trt <- clusterProfiler::gseGO(geneList, OrgDb = "org.Hs.eg.db",
                                              keyType = 'ENSEMBL', ont = "ALL",
                                              nPermSimple = 100000, eps = 0, 
                                              pvalueCutoff = 0.5, pAdjustMethod = "fdr")

Gsea_enrich_down_simp <- clusterProfiler::simplify(gse_cnt_vs_GZMK_trt, 
                                                 cutoff=0.5, by="p.adjust", 
                                                 select_fun= min)

Gsea_enrich_down_simp <- clusterProfiler::setReadable(x = gse_cnt_vs_GZMK_trt, 
                                                    OrgDb = "org.Hs.eg.db", 
                                                    keyType = 'ENSEMBL')

clusterProfiler::dotplot(Gsea_enrich_down_simp, showCategory =10, split=".sign",
                         title = "Live") + facet_grid(.~.sign) #plot to show activated and Rpressed process

## GSEA plot ##
ewp_gse_gzmk <- dplyr::arrange(Gsea_enrich_down_simp, desc(abs(NES))) %>%
  
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
  
  ggtitle("Live  GSEA")



##### Pathway enrichment analysis #####
geneList2 <- Filter_data[, 2] ### logFc
names(geneList2) = as.character(Filter_data[, 6]) ### naming logFc
geneList2 <- sort(geneList2, decreasing = TRUE)


###### Pathway enrichment analysis using enrichkegg ###
enrich_kegg_cnt_vs_GZMK_down <- clusterProfiler::enrichKEGG(Filter_data_down$entrezgene_id, 
                                                          keyType = "kegg", organism = "hsa", 
                                                          pvalueCutoff = 0.5,
                                                          pAdjustMethod = "fdr", 
                                                          qvalueCutoff = 0.2)


enrich_kegg_down <- clusterProfiler::setReadable(x = enrich_kegg_cnt_vs_GZMK_down,
                                               OrgDb = "org.Hs.eg.db", 
                                               keyType = "ENTREZID")


enrich_kegg_down <- dplyr::mutate(enrich_kegg_down, 
                         richFactor = Count / as.numeric(sub("/\\d+", "", 
                                                             BgRatio)))

ggplot(enrich_kegg_down@result[c(1:20), ], showCategory = 10, 
       aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("KEGG Pathways")+
  
  ggtitle("KEGG pathway Live down regulated genes")




###### Pathway enrichment analysis using gsekegg ###
gse_kegg_cnt_vs_GZMK <- clusterProfiler::gseKEGG(geneList= geneList2, 
                                                 organism = "hsa", pvalueCutoff = 0.5, 
                                                 verbose= FALSE, keyType = "kegg", 
                                                 nPermSimple = 100000, eps = 0,
                                                 pAdjustMethod = "fdr", by = "fgsea")

gse_kegg_cnt_vs_GZMK <- clusterProfiler::setReadable(x = gse_kegg_cnt_vs_GZMK, 
                                                     OrgDb = "org.Hs.eg.db",
                                                     keyType = "ENTREZID")


## GSEA plot ##
ewp_gse_kegg_gzmk <- dplyr::arrange(gse_kegg_cnt_vs_GZMK, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:10)

ggplot(ewp_gse_kegg_gzmk, showCategory=20,
       aes(NES, fct_reorder(Description, NES), fill=qvalue)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(12) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle("GSEA KEGG Live ")



####reactome###

enrich_react_1_down <- ReactomePA::enrichPathway(Filter_data_down$entrezgene_id,
                                                 organism = "human", readable = TRUE, 
                                                 pvalueCutoff = 0.5)

enrich_react_down <- dplyr::mutate(enrich_react_1_down, 
                            richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(enrich_react_down@result[c(1:20), ], showCategory = 10, 
       aes(richFactor, fct_reorder(Description, richFactor))) +
  
  geom_segment(aes(xend=0, yend = Description)) +
  
  geom_point(aes(color=p.adjust, size = Count)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(20) +
  
  xlab("Rich Factor") + ylab("Reactome pathways")+
  
  ggtitle("Reactome Live down regulated genes")


#### gsea

gsea_react_1 <- ReactomePA::gsePathway(geneList2, organism = "human", 
                                       pvalueCutoff = 0.5, verbose= FALSE, 
                                       eps = 0, pAdjustMethod = "fdr", by = "fgsea")

Gsea_react_2 <- clusterProfiler::setReadable(x = gsea_react_1, 
                                             OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

clusterProfiler::dotplot(Gsea_react_2, 
                         showCategory =5, split=".sign", 
                         title = "GSEA reactome Live ") + 
  facet_grid(.~.sign) #plot to show activated and Rpressed process


## GSEA plot ##
ewp_gse_react <- dplyr::arrange(Gsea_react_2, desc(abs(NES))) %>%
  
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
  
  ggtitle("GSEA reactome Live ")


##### Merge csv files to one xls ##

count_data <- rownames_to_column(feature_counts_1)

out <- list(Degs_obtained = cntrl_vs_nsc_merge, 
            Significant_Degs= Filter_data, 
            Up= Filter_data_up,
            Down= Filter_data_down,
            Go_up = GO_enrich_1_tmrvsFT_up@result, 
            Go_down = GO_enrich_1_tmrvsFT_down@result, 
            Go_gsea = Gsea_enrich_up_simp@result,
            Kegg_up = enrich_kegg_up@result, 
            Kegg_down = enrich_kegg_down@result,
            kegg_gsea = gse_kegg_cnt_vs_GZMK@result, 
            Reactome_up = enrich_react_1_up@result, 
            Reactome_down = enrich_react_1_down@result, 
            Reactome_gsea = Gsea_react_2@result)

out <- list(Significant_Degs= Filter_data, 
            Up= Filter_data_up,
            Down= Filter_data_down,Go_down = GO_enrich_1_tmrvsFT_down@result, 
            Kegg_down = enrich_kegg_down@result,
            Reactome_down = enrich_react_1_down@result)

writexl::write_xlsx(out, "/Users/Tyagi/Desktop/pca_matia/correctionofimage/mean_163_genes/Live/163_live_enrichment.xlsx")


#To obrain unique from the common

R_only_live_18 <- dplyr::setdiff(Sup_gene_mapped_to_R$external_gene_name, 
                                    Sup_260_union_92$external_gene_name)

R_only_live_18 <- as.data.frame(R_only_live_18)
colnames(R_only_live_18) <- "external_gene_name"
#using semi join to obtain the value

