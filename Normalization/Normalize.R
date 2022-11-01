#reading Data_files
Data_file=read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/GSE190555_12Z_RNA_raw_gene_counts.csv",sep=",",header = T,row.names = 1)
Data_file #printing the file

#Create a count per matrix
cpm=Data_file
for(i in 1:ncol(Data_file)){
  cpm[,i]=(Data_file[,i]/sum(Data_file[,i]))*1000000
}

#Calculate a log of cpm
log_c=log2(cpm+1)
summary(log_c)

#Calculate Zscore
library(matrixStats)     #Importing Library
Zscore = (log_c - rowMeans(log_c))/rowSds(as.matrix(log_c))[row(log_c)]
Zscore     #printing z_score

#Calculate variance using log 
var_c1 = apply(log_c, 1, var)
var_c2 = sort(var_c1,decreasing = T)
top50 = var_c2[1:50]
pmat = Zscore[names(top50),]
data_mat=data.matrix(pmat)
data_mat