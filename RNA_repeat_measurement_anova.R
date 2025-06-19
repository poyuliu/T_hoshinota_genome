rm(list=ls())
library(MARco)
files <- list.files(pattern = "htseq.txt")
htseq <- lapply(files, read.delim, header=F)

data <- merge.data.frame(htseq[[1]],htseq[[2]],by = "V1")
for(i in 3:length(htseq)){
  data <- merge.data.frame(data,htseq[[i]],by = "V1")
}

data <- data.frame(data[-c(1:5),-1],row.names = data$V1[-c(1:5)])
data <- data[rowSums(data)>0,]
names(data) <- gsub("_htseq.txt","",files)
OA <- data[,1:28]
OA <- OA[rowSums(OA)>0,]
TM <- data[,29:52]
TM <- TM[rowSums(TM)>0,]
TM.gps <- as.factor(substr(names(TM),1,7))
OA.gps <- as.factor(substr(names(OA),1,6))

all_annotations <- readRDS("all_annotations.RDS")

# TM
tm.aov.summary <- list()
for(i in 1:nrow(TM)){
  X <- data.frame(Temperature=as.factor(substr(as.character(TM.gps),1,4)),
                  Day=as.numeric(substr(as.character(TM.gps),7,7)),
                  id=factor(substr(names(TM),9,9)),
                  TM=log1p(t(TM[i,])[,1]))
  
  with(X, interaction.plot(Day, Temperature, TM,
                           lty = c(1, 12), lwd = 3,
                           ylab = "log(mean of abundance)", xlab = "Day", trace.label = "Temperature"))
  
  
  tm.aov <- aov(TM ~ Temperature * Day + Error(id),
                data = X)
  tm.aov.summary[[i]] <- summary(tm.aov)
}

int.p <- sapply(sapply(sapply(tm.aov.summary,"[[",2),"[",5),'[',3)
int.a <- p.adjust(int.p,"fdr")

TM_time.mean.sig <- MARco::aggregate2df(t(TM[which(int.p<0.05),]),TM.gps,mean)
XXX <- pheatmap::pheatmap(TM_time.mean.sig[,c(4,5,6,1,2,3)],scale = "row",
                          cellwidth = 15,show_rownames = F,
                          cutree_rows = 2,cluster_cols = F)
TM_time.mean.sig <- cbind(TM_time.mean.sig,GO=all_annotations$GO[match(rownames(TM_time.mean.sig),all_annotations$qseqid)],cluster=cutree(XXX$tree_row,k=6))
TM_time.mean.sig <- TM_time.mean.sig[XXX$tree_row$order,]
write.csv(TM_time.mean.sig,"TM_time_mean_sig_GO_05.csv")

# OA
oa.aov.summary <- list()
for(i in 1:nrow(OA)){
  X <- data.frame(OceanAcidification =as.factor(substr(as.character(OA.gps),1,3)),
                  Day=as.numeric(substr(as.character(OA.gps),6,6)),
                  id=factor(substr(names(OA),8,8)),
                  OA=log1p(t(OA[i,])[,1]))
  
  with(X, interaction.plot(Day, OceanAcidification , OA,
                           lty = c(1, 12), lwd = 3,
                           ylab = "log(mean of abundance)", xlab = "Day", trace.label = "OceanAcidification "))
  
  
  oa.aov <- aov(OA ~ OceanAcidification  * Day + Error(id),
                data = X)
  oa.aov.summary[[i]] <- summary(oa.aov)
}

int.p <- sapply(sapply(sapply(oa.aov.summary,"[[",2),"[",5),'[',3)
int.a <- p.adjust(int.p,"fdr")

OA_time.mean.sig <- MARco::aggregate2df(t(OA[which(int.p<0.05),]),OA.gps,mean)
YYY <- pheatmap::pheatmap(OA_time.mean.sig[,c(4,5,6,1,2,3)],scale = "row",
                          cellwidth = 15,show_rownames = F,
                          cutree_rows = 2,cluster_cols = F)
OA_time.mean.sig <- cbind(OA_time.mean.sig,GO=all_annotations$GO[match(rownames(OA_time.mean.sig),all_annotations$qseqid)],cluster=cutree(YYY$tree_row,k=6))
OA_time.mean.sig <- OA_time.mean.sig[YYY$tree_row$order,]
write.csv(OA_time.mean.sig,"OA_time_mean_sig_GO_05.csv")
