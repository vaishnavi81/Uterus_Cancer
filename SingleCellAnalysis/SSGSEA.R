#single sample geneset enrichment
#Import Libraries
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(data.table)

#Ssgea function
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}
library(org.Hs.eg.db)
library(AnnotationDbi)
#Read dataset
data1 = read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/Gene_names.csv")
rownames(data1)=make.names(data1[,1],unique = TRUE)
View(data1)
data1=data1[,-1]
head(data1)
summary(data1)

data1 = as.matrix(data1)
#Import the markers dataset
gene_set = read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/markers2Sep.csv")
head(gene_set)

#Convert to list
gene_sets = as.list(as.data.frame(gene_set))
print("genes set ready")

#Use ssgea function 
res = ssgsea(data1, gene_sets, scale = TRUE, norm = FALSE)

#transpose the res
res1 = t(res)
head(res1)

#Convert to matrix
mat = as.matrix(res1)
for(i in 1:nrow(mat))
{ 
  vec = as.numeric(c(mat[i,]))
  mat[i,1:ncol(mat)] = (vec-mean(vec))/sd(vec)
  
}

#Plot the heatmap
library(ComplexHeatmap)
Heatmap(t(mat),col = colorRamp2(c(-2,0,2),c("green","Purple","orange")))

