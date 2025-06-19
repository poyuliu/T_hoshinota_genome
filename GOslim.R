data <- read.delim("interproscan_go_Terpios_hoshinota")
GO <- data$V14
GO <- strsplit(x = GO,split = "|",fixed = T)

GO <- sapply(GO, function(y) unique(sapply(y, function(x) substr(x,1,10))) )


length(unique(unlist(GO)))


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")

library(GSEABase)

myIds <- unique(unlist(GO))
myIds <- myIds[!is.na(myIds)]
myIds <- myIds[myIds!="-"]

myCollection <- GOCollection(myIds)
slim <- getOBOCollection("goslim_pir.obo")
go_slim_counts_MF <- goSlim(myCollection, slim, "MF")
go_slim_counts_BP <- goSlim(myCollection, slim, "BP")
go_slim_counts_CC <- goSlim(myCollection, slim, "CC")
go_slim_counts_MF$GO_ontology <- "MF"
go_slim_counts_BP$GO_ontology <- "BP"
go_slim_counts_CC$GO_ontology <- "CC"
go_slim_counts <- rbind(go_slim_counts_MF,go_slim_counts_BP,go_slim_counts_CC)

go_slim_counts <- go_slim_counts[go_slim_counts$Count>0,]

plot.data <- data.frame(Counts=go_slim_counts$Count,row.names =paste(rownames(go_slim_counts),go_slim_counts$Term))

# Plot a bar chart of GO Slim categories
barplot(t(plot.data), main = "GO Slim Category Distribution", las = 2, ylab = "Count", col = "lightblue")

write.csv(go_slim_counts,"go_slim_counts.csv")

sum(go_slim_counts$Count)
