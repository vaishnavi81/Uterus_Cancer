#To identify genes which are differential in tumor vs control samples

data1=matrix(NA,ncol=4,nrow = nrow(log_c))
rownames(data1)=rownames(log_c)
colnames(data1)=c('meanTumor','meanControl','pvalue','log2FC')
data1
for(i in 1:nrow(log_c)){
  vector1 = as.numeric(log_c[i, 1:3])
  
  vector2 = as.numeric(log_c[i, 4:6])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  data1[i,1]=res$estimate[[1]]
  data1[i,2]=res$estimate[[2]]
  data1[i,3]=res$p.value
  data1[i,4]=data1[i,1]-data1[i,2]
  
}

data1=as.data.frame(data1)
num=which(is.nan(data1$pvalue))
data1[num,'pvalue']=1

library(EnhancedVolcano) #importing library
pdf('volcano.pdf',width=10,height=10)
EnhancedVolcano(data1,lab = rownames(data1),x = 'log2FC' ,y ='pvalue',pCutoff = 3.5)
dev.off()
saveRDS(data1,file = "C:/Users/91735/Documents/BVP material/7_cancer genomics/DEG_data.rds")
M=readRDS('DEG_data.rds')
View(M)
