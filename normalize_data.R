data=read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/uterus_cancer.csv")
data
data_1= data[,-1]   #remove 1 column
data_1
mat= as.matrix(data_1)
rownames(mat)= data[,1]    #matrix
View(mat)


for(i in 1:ncol(data_1)){
  mat[, i] = (data_1[, i]/sum(data1[,i]))*1000000   #normalize the data
}
l2cpm=log2(mat + 1)
l2cpm
