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
#Create a heatmap
library(ComplexHeatmap)    #importing library
pdf('data.pdf',width=10,height=10)
Heatmap(data_mat)
dev.off()


#To identify genes which are differential in tumor vs control samples

data1=matrix(NA,ncol=4,nrow = nrow(log_c))
rownames(data1)=rownames(log_c)
colnames(data1)=c('meanTumor','meanControl','pvalue','log2FC')
data1
for(i in 1:nrow(log_c)){
  vector1 = as.numeric(log_c[i, 1:7])
  
  vector2 = as.numeric(log_c[i, 8:11])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  data1[i,1]=res$estimate[[1]]
  data1[i,2]=res$estimate[[2]]
  data1[i,3]=res$p.value
  data1[i,4]=data1[i,1]-data1[i,2]
  
}

data1=as.data.frame(data1)
num=which(is.nan(mat1$pvalue))
data1[num,'pvalue']=1

library(EnhancedVolcano) #importing library
pdf('data.pdf',width=10,height=10)
EnhancedVolcano(data1,lab = rownames(data1),x = 'log2FC' ,y ='pvalue')
dev.off()
