#Create an User Defined function to create a logCPM of an input data 

library(matrixStats)  #importing_library

data=function(x){        #create_function(X)
  data_file=x         #storing it into variable
  
  for (i in 1:ncol(x)) {
    data_file[,i] = (x[,i]/sum(x[,i]))*1000000   #storing Log2cpm data 
    print(head(data_file))
    data_file[,i]= log2(data_file[,i] +1)
    logfc=log2(data_file+1)
  }
  return(data_file)
}

#store (input-file into X and read )

x=read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/GSE190555_12Z_RNA_raw_gene_counts.csv",sep=",",header=T,row.names = 1)
F_data=data(x)

F_data    #log2cpm data stored in this variable
