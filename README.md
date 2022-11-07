# UTERUS CANCER

Gene Expression Analysis

Gene expression analysis is most simply described as the study of the way genes are transcribed to synthesize functional 
gene products functional RNA species or protein products.
The most common laboratory methods used to measure GE levels are Northern blotting, quantitative polymerase chain reaction (qPCR), DNA microarray, and RNA-Seq.
For gene expression GEO database and my case study has Uterus cancer sequence collected from NCBI GEO database

File:  GSE190555_12Z_RNA_raw_gene_counts.csv 

1)NORMALIZATION OF THE DATA.
Data normalization is the process of reorganizing data within a database so that users can utilize it for further queries and analysis.

1)Count CPM
Formula for Cpm  cpm[,i]=(Data_file[,i]/sum(Data_file[,i]))*1000000

2)Calculate Z_score
Zscore = (log_c - rowMeans(log_c))/rowSds(as.matrix(log_c))[row(log_c)]

3)Calculate variance using log

2)HEATMAP
Heat map analysis is the process of reviewing and analyzing heat map data to gather insights about user interaction and behaviour as they engage with your product. 
This data analysis can lead to improved site designs with lower bounce rates, reduced churn, fewer drop-offs, more pageviews, and better conversion rates.

3)DEG
Differential expression analysis means taking the normalised read count data and performing statistical analysis to 
discover quantitative changes in expression levels between experimental groups.

4)CIRCOS PLOT:
 Circos is a software package for visualizing data in a circular layout. This makes Circos ideal for exploring relationships between objects or positions. Circos plots have appeared in thousands of scientific publications.

For Circos plot we use DEG file 
We build the basic chromosome ideogram and file that contains log2 fold change values and p values.
Convert Enseemble_id name with bimaRT 'hgnc_symbol' gene name
Track upsiginificant genes in the chromosome

5)SSGSEA ANALYSIS:
Single-sample GSEA (ssGSEA), an extension of Gene Set Enrichment Analysis (GSEA), calculates separate enrichment scores for each pairing of a sample and gene set. Each ssGSEA enrichment score represents the degree to which the genes in a particular gene set are coordinately up- or down-regulated within a sample.
1)Create a ssgsea function
Use log2 normalized data into an object for further analysis.
Zscore use for the ssgsea output for comparative analysis.

6)SURVIVAL ANALYSIS:
File:  Breast_Cancer.csv
Survival analysis is concerned with studying the time between entry to a study and a subsequent event. 
Originally the analysis was concerned with time from treatment until death, hence the name, but survival analysis is applicable to many areas as well as mortality.








