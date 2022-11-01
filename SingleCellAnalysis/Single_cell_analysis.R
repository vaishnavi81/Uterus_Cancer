#ScRNA-Seq
#install.packages("dplyr")
#install.packages("patchwork")
#install.packages("Seurat")
library(Seurat)
library(patchwork)
library(dplyr)
#read covearge=how manty times you need to each read again again like 10X(10 times)
# Load the PBMC dataset
pbmc.data = Read10X(data.dir = "C:/Users/91735/Documents/BVP material/7_cancer genomics/SCA")
pbmc.data
# Initialize the Seurat object with the raw (non-normalized data).

pbmc = CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
pbmc
## 50 x 10 sparse Matrix of class "dgCMatrix"
pbmc.data[1:50, 1:10]

#delete mitochondiral genes because 1 they are apoctosis and 2)replicate and it means its not a single cell
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)
#1)nCount_RNA   2)nFeature_RNA   3)percent.mt
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc

pbmc = NormalizeData(pbmc)
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot2
#downstream analyses, so that highly-expressed genes do not dominate
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

## Centering and scaling data matrix
pbmc@assays$RNA@scale.data[1:50, 1:5]

#Next we perform PCA on the scaled data. (dimentionality_reduction)
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#to know about components in your data
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## Computing nearest neighbor graph
pbmc = FindNeighbors(pbmc, dims = 1:10)

## Computing SNN
pbmc = FindClusters(pbmc, resolution = 0.5)

head(pbmc@meta.data)

pbmc = RunUMAP(pbmc, dims = 1:10)
#alwasy go hierharcical clusitering 
DimPlot(pbmc, reduction = "umap")

DimPlot(pbmc, reduction = "umap", label = T)

#only need upregulated not(downregulated beacuse of they not express)
pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(pbmc.markers)

a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a

genes = a %>% pull(gene)
genes

FeaturePlot(pbmc, features = genes[1:2])

FeaturePlot(pbmc, features = genes[1:2], cols = c("white", "red"))
