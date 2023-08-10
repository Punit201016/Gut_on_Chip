### GSEA Andrey ####

LMD_epith_normalzied_gene_expression <- read_excel("C:/Users/Tyagi/Desktop/LMD-epith-normalzied-gene expression.xlsx")

ensembl_mouse_human_gene_symbol <- read_excel("C:/Users/Tyagi/Desktop/ensembl_mouse_human_gene_symbol.xlsx")

Mouse_to_human_symbol_lmd <- merge(LMD_epith_normalzied_gene_expression, ensembl_mouse_human_gene_symbol,
                                   by.x = "Symbol", by.y = "Mouse gene name") 
Mouse_to_human_symbol_lmd <- Mouse_to_human_symbol_lmd[c(1:24040), ]

Mouse_to_human_symbol_lmd <- Mouse_to_human_symbol_lmd[, c(26, 25, 1, 10:22)]
Mouse_to_human_symbol_lmd <- Mouse_to_human_symbol_lmd[, c(1,5:16)]
Mouse_to_human_symbol_lmd <- Mouse_to_human_symbol_lmd[, c(1:5, 12,13, 6:11)]
colnames(Mouse_to_human_symbol_lmd)[2:7] <- "C"
colnames(Mouse_to_human_symbol_lmd)[8:13] <- "GF"

#Mouse_to_human_symbol_lmd <- column_to_rownames(Mouse_to_human_symbol_lmd, var = "Gene name")

Mouse_to_human_symbol_lmd_aggre <- Mouse_to_human_symbol_lmd %>% dplyr::group_by("Gene name") %>% 
                                summarize(value=mean(Mouse_to_human_symbol_lmd[, c(2:6)]))

gsea_result <- toupper()

  
#aggregate(x = Mouse_to_human_symbol_lmd[2:13], by=list(Mouse_to_human_symbol_lmd$`Gene name`),FUN = mean)


writexl::write_xlsx(Mouse_to_human_symbol_lmd, "C:/Users/Tyagi/Desktop/lmd_new_map.xlsx")

### Live_vs_Control_go_excel_file

sheet_names1 <- excel_sheets("C:/Users/Tyagi/Desktop/pca_matia/pca_after_comments/Analysed_data/Degs_analysis/live_and_sup_vs_control/degs_live_vs_control_cmn_pt_all/degs_pathway_analaysis_live_vs_control.xlsx ")           # Get sheet names


list_all <- lapply(sheet_names1, function(x)         # Read all sheets to list
  as.data.frame(read_excel("C:/Users/Tyagi/Desktop/degs_pathway_analaysis_live_vs_control.xlsx", sheet = x)))
    
sheet_names1 <- c("go_up","Kegg_up", "Reactome_up")
sheet_names1 <- c("go_down","Kegg_down", "Reactome_down")  
sheet_names1 <- c("go_up", "go_down", "Kegg_up","Kegg_down", "Reactome_up","Reactome_down")  


sheet_names1 <- c("go_gsea","kegg_gsea","Reactome_gsea")

path_factor_go_down <- lapply(list_all, function(x) x %>% select(ID, Description,
                                                        qvalue,
                                                        geneID))
path_factor_go_up <- do.call(rbind, path_factor_go_up)
path_factor_go_down <- do.call(rbind, path_factor_go_down)


path_factor_up <- path_factor_go_up %>% filter(qvalue < 0.01)
path_factor_up$UP <- "UP"

path_factor_down <- path_factor_go_down %>% filter(qvalue < 0.01)
path_factor_down$Down <- "DOWN"

# merging the sig up BP with the mapped mouse data#
down_gea_new <- read.delim("C:/Users/Tyagi/Desktop/down_gea_new.txt")

up_gsea_new <- read.delim("C:/Users/Tyagi/Desktop/up_gsea_new.txt")

## merginf new up and down gsea result ##
final_gsea_new <- rbind(up_gsea, down_gea)

## new 
final_gsea_new_1 <- rbind(up_gsea_new, down_gea_new)

path_factor_up <- merge(path_factor_up, final_gsea_new_1, by.x = "ID", by.y = "NAME")
path_factor_up <- path_factor_up[!duplicated(path_factor_up$Description), ]
path_factor_up$ES_h <- "1" 


path_factor_down <- merge(path_factor_down, final_gsea_new_1, by.x = "ID", by.y = "NAME")
path_factor_down <- path_factor_down[!duplicated(path_factor_down$Description), ]
path_factor_down$ES_h <- "-1" 


## Selecting according to the highest enrichment score
path_factor_up_1 <- path_factor_up[order(path_factor_up$ES, decreasing = T), ] 
path_factor_up_1<-as.data.frame(path_factor_up_1[c(1:15), c(2,3,8)])
colnames(path_factor_up_1)[2] <- "qvalue"

path_factor_down_1 <- path_factor_down[order(path_factor_down$ES, decreasing = T), ] 
path_factor_down_1<-as.data.frame(path_factor_down_1[c(1:15), c(2,3,8)])
colnames(path_factor_up_1)[2] <- "qvalue"

Final_human_mouse_gsea_new <- rbind(path_factor_up_1, path_factor_down_1)

gsea_and_live_merge <- merge(final_gsea_new_1, Final_human_mouse_gsea_new,
                             by.x = "NAME", by.y = "ID")

writexl::write_xlsx(list(final_gsea_new_1,path_factor_go_up, path_factor_go_down 
                         ), "C:/Users/Tyagi/Desktop/final_mouse_gsea_p_value_1.xlsx")

writexl::write_xlsx(Final_human_mouse_gsea_new, "C:/Users/Tyagi/Desktop/final_mouse_gsea_p_value.xlsx")

#for visualization##
path_sig <-arrange(df_plot_1, desc(abs(NES))) %>%
  
  group_by(sign(NES)) %>%
  
  dplyr::slice(1:30)

path_sig <- path_sig %>% drop_na
colnames(path_sig)[3] <- "FDR"
path_sig <- path_sig[c(1:8, 31:39, 42,44,48), ]

ggplot(path_sig, showCategory=10,
       aes(NES, fct_reorder(Description.y, NES), fill=FDR)) +
  
  geom_col() +
  
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       
                       guide=guide_colorbar(reverse=TRUE)) +
  
  theme_dose(14) +
  
  xlab("Normalized Enrichment Score") +
  
  ylab(NULL) +
  
  ggtitle(label = "Biological process in Live vs Mouse GSEA")



# order by q value for whole for gsea ##
path_factor <- do.call(rbind, path_factor)
path_factor <- path_factor %>% filter(qvalue < 0.01)
writexl::write_xlsx(path_factor, "C:/Users/Tyagi/Desktop/path_factor_for_gsea.xlsx")

path_factor_2 <- lapply(list_all, function(x) x%>% select(ID, Description,
                                                        qvalue,core_enrichment,
                                                        ))
path_factor_2 <- do.call(rbind, path_factor_2)
colnames(path_factor_2)[4] <- "geneID"

path_factor <- rbind(path_factor, path_factor_2)

path_factor_new <- path_factor %>% filter(qvalue < 0.01)


##### merging gsea result of C_GF and our result from live vs control ###

path_factor_new <- lapply(path_factor_new, function(x) x %>% toupper())
path_factor_new <- as.data.frame(path_factor_new)

gsea_result_1 <- gsub("_", " ", gsea_result$NAME)
vec <- c("GOBP", "REACTOME", "HP", "GOMF", "KEGG")
gsea_result_1 <- gsub(vec, c("GOBP"="", "REACTOME"="", "HP"="", "GOMF"=""), gsea_result_1)
gsea_result_1 <- as.data.frame(gsea_result_1)


library(tidyverse)
gsea_matach <- gsea_new_1 %>% mutate(match = map2_lgl(str_extract_all(V1, "\\w+"), 
                          str_extract_all(V3, "\\w+"),  ~ any(.y %in% .x)))

gsea_matach <- path_factor_new %>% mutate(match= str_match_all(Description, "\\w+"), 
str_match_all(gsea_result_1, "\\w+"))

gsea_matach$match2 <- sapply(seq(nrow(gsea_match_11)), function(i) {res <- paste0(as.character(gsea_match_11$match[[i]]), collapse = " ")})

gsea_matach <- merge(str_match_all(path_factor_new, "\\w+"), 
                                                str_match_all(gsea_result_1, "\\w+"), by.x= "Description",
                            by.y= "NAME")





C_sf_merge <- path_factor_new[str_detect(path_factor_new$Description, gsea_react_1), ]


final_gsea <- merge(gsea_c_live, path_factor, by.x="NAME", by.y = "ID")









## R and NR GSEA using LPS ###
library(qusage)
library(msigdbr)
library(magrittr)

inflamatory_lps_gmt <- qusage::read.gmt("C:/Users/Tyagi/Desktop/live_mouse_go/ZHOU_INFLAMMATORY_RESPONSE_LPS_UP.v2023.1.Hs.gmt")

NR <- clusterProfiler::GSEA(geneList = geneList, minGSSize=25,maxGSSize=500, 
                                  seed=TRUE, eps = 0, 
                                  pvalueCutoff = 0.05, pAdjustMethod = "BH",
                            TERM2GENE = inflamatory_lps_gmt %>% dplyr::select(
                              "gene_symbol"
                            )
                                  )


##################### gsea analysis#######################
df1<-readr::read_tsv("C:/Users/Tyagi/gsea_home/output/jul09/INR_2.Gsea.1688929521351/gsea_report_for_GF_1688929521351.tsv")
df2<-readr::read_tsv("C:/Users/Tyagi/gsea_home/output/jul09/INR_2.Gsea.1688929521351/gsea_report_for_C_1688929521351.tsv")

df_1 <- rbind(df1, df2)

df_2 <- merge(df_1, path_factor_go_up, by.x = "NAME",by.y= "ID", all.x = T)
df_2_2 <- 
df_3 <-merge(df_2, path_factor_go_down, by.x = "NAME",by.y= "ID", all.x = T)
df_3$qvalue.x[is.na(df_3$qvalue.x)] <- 1
df_3$qvalue.y[is.na(df_3$qvalue.y)] <- 1

df_3$Delta <- df_3$qvalue.x-df_3$qvalue.y
df_3$up_down <- ifelse(df_3$Delta >0, -1,1)

df_4 <- df_3 %>% filter((ES > 0.5 | ES < -0.4) & up_down == 1) #change filter +1 and -1 for up and down
df_4_1 <- df_3 %>% filter((ES > 0.5 | ES < -0.4) & up_down == -1)
sum(df_4$ES > 0, na.rm = T) #up
sum(df_4$ES < 0, na.rm = T) # down

## for plot##
df_5_plot <- df_4[, c(16,6,8)]
df_6_plot <- df_4[, c(13,6,8)]
colnames(df_6_plot)[1] <- "Description.y"
df_plot_1 <- rbind(df_5_plot, df_6_plot)

write_xlsx(list(df_1, df_3, df_4, df_4_1, df_plot_1), "C:/Users/Tyagi/Desktop/final_result/GSEA_C_vs_GF.xlsx")

