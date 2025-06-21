
# === load topGO package ===
library(topGO)

# === define ontology list ===
ontologies <- c("BP", "MF", "CC")

# === define analysis function ===
run_topGO_all_ontologies <- function(geneList, gene2GO, ontologies, label) {
  all_results <- list()
  
  for (ont in ontologies) {
    message("Running ", label, " - ", ont, "...")
    
    GOdata <- new("topGOdata",
                  ontology = ont,
                  allGenes = geneList,
                  geneSel = function(p) p == 1,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO)
    
    result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    
    num_nodes <- min(100, length(usedGO(GOdata)))
    
    tab <- GenTable(GOdata,
                    weightFisher = result,
                    topNodes = num_nodes)
    
    tab$Ontology <- ont
    all_results[[ont]] <- tab
  }
  
  final <- do.call(rbind, all_results)
  return(final)
}

# === run analysis ===
result_OA_cluster1 <- run_topGO_all_ontologies(geneList_OA_1, gene2GO, ontologies, "OA_cluster1")
result_OA_cluster2 <- run_topGO_all_ontologies(geneList_OA_2, gene2GO, ontologies, "OA_cluster2")
result_TM_cluster1 <- run_topGO_all_ontologies(geneList_TM_1, gene2GO, ontologies, "TM_cluster1")
result_TM_cluster2 <- run_topGO_all_ontologies(geneList_TM_2, gene2GO, ontologies, "TM_cluster2")

# === output CSV ===

write.csv(result_OA_cluster1, file = "OA_cluster1_topGO_all.csv", row.names = FALSE)
write.csv(result_OA_cluster2, file = "OA_cluster2_topGO_all.csv", row.names = FALSE)
write.csv(result_TM_cluster1, file = "TM_cluster1_topGO_all.csv", row.names = FALSE)
write.csv(result_TM_cluster2, file = "TM_cluster2_topGO_all.csv", row.names = FALSE)

