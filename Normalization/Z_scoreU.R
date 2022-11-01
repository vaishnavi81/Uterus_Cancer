   #importing library
library(matrixStats) 

data=source("logcpm.rds")

data2=function(F_data){        #create_function2 and load log2cpm data into it
  for (i in 1:ncol(F_data)){
    z_score = (F_data - rowMeans(F_data))/rowSds(as.matrix(F_data))[row(F_data)]
    
  }   
  
  z_score[is.na(z_score)]=0  
  z_score#calculate Z_score
  zscore = as.matrix(z_score)
  Heatmap(zscore)
  return(Heatmap(zscore[1:20]))                   #calculate the heatmap
  
}
library(ComplexHeatmap)
pdf('data1.pdf',width=10,height=10)
data2(F_data)
dev.off()

