

## WARNING. the original script was run on nder Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01). This is a computationally demanding which is likely not supposed to run on the local computers.



##upload packages

library(openxlsx)
library(data.table)
library(magrittr)
library(GGally)

##upload packages


## set working directory

setwd("/home/dzavadska/finalrun/")

## set working directory


## upload and prepare data. more detailed comments below each line.

taxa_SMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs.xlsx", sheet = 4, startRow = 1, colNames = TRUE) %>% data.frame
#taxa_SMAGs - contains info on which SMAG belongs to which taxogroup

tmp <- data.frame(levels(factor(taxa_SMAGs[,2]))) 
#subset taxogroup colunm and assign factor to it

taxa_SMAGs_melted <- merge(taxa_SMAGs,tmp,by.x="AA_Best_taxonomy_GENRE",by.y="levels.factor.taxa_SMAGs...2...",all.x = T)
taxa_SMAGs_melted <- taxa_SMAGs_melted[,c("AA_Best_taxonomy_GENRE", "AA_SMAG")]
#Now we have a version of taxa SMAGs where each row contains a SMAG ID and a taxogroup where it is assigned

vocabulary <- fread("vocabulary_taxa.csv",sep=",")
taxa_list <- c(levels(factor(taxa_SMAGs_melted[,1]))[2:length(levels(factor(taxa_SMAGs_melted[,1])))])
list_of_all <- (c(taxa_list, c(unique(vocabulary$taxa_SMAGs))))
list_of_all <- list_of_all[duplicated(c(taxa_list, c(unique(vocabulary$taxa_SMAGs))))==T]

#uncomment to test the example barcode-MAG pairs

#i <- "New_Choanozoa_01"

#uncomment to test the example barcode-MAG pairs


#starting the iteration over taxogroup

for(i in list_of_all[1:(length(list_of_all))]) {
  
  dataSMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs.xlsx", sheet = 7, startRow = 1, colNames = TRUE) %>% data.table
  dataSMAGs <- dataSMAGs[1:713,]
  #calculating the relative abundance of SMAGs
  
  numeric_cols <- names(dataSMAGs)[sapply(dataSMAGs, is.numeric)]
  dataSMAGs_relabund <- dataSMAGs[, lapply(.SD, function(x) x/sum(x)), .SDcols = numeric_cols]
  dataSMAGs <- cbind(dataSMAGs[,1],dataSMAGs_relabund)
  #calculating the relative abundance of SMAGs
  
  
  #subsetting SMAGs dataset by the taxogroup
  
  tokeep <- c(taxa_SMAGs_melted[which(taxa_SMAGs_melted$AA_Best_taxonomy_GENRE == i),2])
  
  #uncomment to test the example barcode-MAG pairs
  #tokeep <- c("TARA_ARC_108_MAG_00240","TARA_ARC_108_MAG_00247","TARA_ARC_108_MAG_00250","TARA_ARC_108_MAG_00255",
  #            "TARA_ARC_108_MAG_00257","TARA_ARC_108_MAG_00260","TARA_MED_95_MAG_00420","TARA_MED_95_MAG_00435",
  #           "TARA_PSE_93_MAG_00240","TARA_SOC_28_MAG_00054","TARA_SOC_28_MAG_00062")
  #uncomment to test the example barcode-MAG pairs
 
  dataSMAGs <- dataSMAGs[SMAG%in%tokeep]
  
  #subsetting SMAGs dataset by the taxogroup
  
  
  #reshaping the dataframe and matching station IDs with metagenome IDs
  
  dataSMAGs_melted <- melt(dataSMAGs,id.vars="SMAG")
  
  context_SMAGs <- fread("TAGs-18S-V4_NAME-PROJ-MAT-READ-CAB_nico.list",header = F)
  tmp <- fread("TARA_CONTEXT_95_SEQ_STATIONS.txt.gz")
  tmp <- tmp[,list(Sample.ID,Sample.mat)]
  tmp[,Sample.ID:=sub("TARA_","",Sample.ID)]
  tmp <- merge(context_SMAGs,tmp,by.x="V5",by.y="Sample.ID",all.x=T)
  tmp <- tmp[,list(V1,Sample.mat)]
  tmp <- unique(tmp)
  tmp <- tmp[!duplicated(V1)]
  dataSMAGs_melted <- merge(dataSMAGs_melted,tmp,by.x="variable",by.y="V1",all.x = T)
  dataSMAGs_df <- dcast(Sample.mat~SMAG,value.var="value",data=dataSMAGs_melted[!is.na(Sample.mat)],fun.aggregate = mean) %>%
    data.frame(row.names = "Sample.mat")
  
  #reshaping the dataframe and matching station IDs with metagenome IDs
  
  #writing out intermediate dataset of SMAGs abundances
  
  fwrite(dataSMAGs_df, paste0(i, "data_SMAGs.tsv"),sep="\t",row.names = T)
  
  #writing out intermediate dataset of SMAGs abundances
  
  #reading and reshaping V9 dataset
  
  dataV9 <- fread("globaldataset.otus.v20171106.withfunctions.pub.tsv")
  names(dataV9)[-1] = sub('.*\\;', "", names(dataV9)[-1]) ####This line is actually essential for the official dataset with V9 data
  dataV9_melted <- dataV9[,.SD,.SDcols=c("md5sum",grep("TV9|TA9",colnames(dataV9),value=T))] %>%
    melt(id.vars="md5sum")
  
  #reading and reshaping V9 dataset
  
  #counting relative abundances of V9 metabarcodes
  
  dataV9_melted[,value2:=value/sum(value),by=variable]
  
  #counting relative abundances of V9 metabarcodes
  
  #subsetting V9 dataset by the abundance and taxogroup
  
  taxogroup_i <- vocabulary[which(vocabulary == i)[1],2]
  #dataV9 <- dataV9[taxogroup=="Choanoflagellatea"&pid>90]
  dataV9 <- dataV9[taxogroup==c(taxogroup_i)[1]&pid>90]
  
  x <- apply(dataV9[,.SD,.SDcols=grep("TV9|TA9",colnames(dataV9),value=T)],1,function(X) sum(X>0))
  dataV9 <- dataV9[x>100]
  dataV9_melted <- dataV9_melted[md5sum%in%dataV9[,md5sum]]
  
  #subsetting V9 dataset by the abundance and taxogroup
  
  #matching V9 metabarcoding IDs with staion/sample IDs
  
  contextV9 <- fread("TARA_CONTEXT_95_SEQ_STATIONS_V9.tsv")
  dataV9_melted <- merge(dataV9_melted,contextV9[,list(Sample_seq_id,Sample.mat)],by.x="variable",by.y="Sample_seq_id")
  dataV9_melted[,variable:=NULL]
  dataV9_df_relabund <- dcast(Sample.mat~md5sum,value.var="value2",data=dataV9_melted,fun.aggregate = mean) %>%
    data.frame(row.names = "Sample.mat",check.names = F)
  
  #matching V9 metabarcoding IDs with staion/sample IDs
  
  # writing out intermediate dataset of V9 abundances
  
  fwrite(dataV9_df_relabund,paste0(i, "data_V9_relabund.tsv"),sep="\t",row.names = T)
  
  # writing out intermediate dataset of V9 abundances
  
  
  #######################Correlation estimation##############
  tmp <- dataSMAGs_df
  
  ######################################################################################################################################3
  #tmp[tmp==0] <- NA
  ######################################################################################################################################
  
  tmp <- merge(tmp,dataV9_df_relabund,by="row.names") %>%
    data.frame(row.names = "Row.names",check.names = F)
  
  tmp <- na.omit(tmp)
  
  ###original
  #data_cor <- cor(tmp[,grep("TARA",colnames(tmp))],tmp[,grep("TARA",colnames(tmp),invert = T)],method="spearman",use="pairwise.complete.obs") %>%
  #  melt %>% data.table
  #####
  ####corrected for my machine
  data_cor <- cor(tmp[,grep("TARA",colnames(tmp))],tmp[,grep("TARA",colnames(tmp),invert = T)],method="spearman",use="pairwise.complete.obs") #%>%  melt %>% data.table
  data_cor <- as.data.table(as.table(data_cor))
  #data_cor <- melt(data_cor)
  ####
  ###original
  #data_cor <- merge(data_cor,dataV9[,list(md5sum,pid,lineage,sequence,refs,ctotab)],by.x="Var2",by.y="md5sum")
  ####
  ####corrected for my machine
  data_cor$V2 = sub('.x', "", data_cor$V2)
  data_cor$V2 = sub('.y', "", data_cor$V2)####for some reason can add those stange .x and .ys
  data_cor <- merge(data_cor,dataV9[,list(md5sum,pid,lineage,sequence,refs,ctotab)],by.x="V2",by.y="md5sum")
  ####
  
  setnames(data_cor,c("md5sum","SMAG","rho","pid","lineage","sequence","refs","abundance"))
  
  #=========Counting number of points per plot==========#
  
  #
  df1 <-data.frame(matrix(nrow=length(data_cor$md5sum)))
  #create an empty dataframe where amount of co-occurence points for each barcode-MAG pair is recorded
  if (length(levels(as.factor(data_cor$SMAG))) == 1) {
    
    for (e in c(1:length(data_cor$md5sum))) {
      x <- which(tmp[,data_cor$md5sum[e]] != 0)
      y <- which(tmp[,1] != 0)
      #length(intersect(x,y))df1
      df1$SMAG[e] <- tmp[,1]
      df1$md5sum[e] <- data_cor$md5sum[e]
      df1$stationcount[e] <- length(intersect(x,y))
    }
    
    data_cor$station_number[1] <- 1
    #create a column in the dataframe that will be written to file. The number of points will be stored in that column.
    for (e in c(1:length(data_cor$md5sum))) {
      w <- which(data_cor$md5sum==df1$md5sum[e])
      z <- 1
      a <- df1$stationcount[e]
      b <- intersect(w,z)
      data_cor$station_number[b] <- a   
    }
    
  } else{
    
    for (e in c(1:length(data_cor$md5sum))) {
      x <- which(tmp[,data_cor$md5sum[e]] != 0)
      y <- which(tmp[,data_cor$SMAG[e]] != 0)
      #length(intersect(x,y))df1
      df1$SMAG[e] <- data_cor$SMAG[e]
      df1$md5sum[e] <- data_cor$md5sum[e]
      df1$stationcount[e] <- length(intersect(x,y))
    }
    
    data_cor$station_number[1] <- 1
    #create a column in the dataframe that will be written to file. The number of points will be stored in that column.
    for (e in c(1:length(data_cor$md5sum))) {
      w <- which(data_cor$md5sum==df1$md5sum[e])
      z <- which(data_cor$SMAG==df1$SMAG[e])
      a <- df1$stationcount[e]
      b <- intersect(w,z)
      data_cor$station_number[b] <- a   
    }
    
  }
  
  
  #=========Finished counting number of points per plot==========#
  
  #writing final output table
  
  data_cor[order(rho,decreasing = T)][1:10]
  fwrite(data_cor,paste0(i, "all_scores_relabund.tsv"),sep="\t")
  split(data_cor,data_cor[,SMAG]) %>%
    lapply(function(X){
      X[order(rho,decreasing = T)][1:3]
    }) %>% rbindlist() %>% fwrite(paste0(i, "best_3_scores_relabund.tsv"),sep="\t")
  
  #writing final output table

  
}


  