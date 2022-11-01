library("RCircos")
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info=UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

cyto.info = UCSC.HG38.Human.CytoBandIdeogram
cyto.info$Name = NA
cyto.info$Stain = NA
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

chr_order = unique(cyto.info$Chromosome)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

ideo = RCircos.Get.Plot.Ideogram()
ideo$BandColor = 'salmon'
num = which(ideo$Chromosome == 'chrX')
ideo[num, 'BandColor'] = 'chartreuse'

num = which(ideo$Chromosome == 'chrY')
ideo[num, 'BandColor'] = 'purple'


RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()


library(biomaRt)

mat = read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/Enhanced_DEG.csv")
rownames(mat)=make.names(mat[,1],unique = TRUE)
mat=mat[,-1]

m = useMart('ensembl', dataset='hsapiens_gene_ensembl')

coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(mat)),
               mart = m)

write.csv(coords, file = 'coords.csv')
coords$chromosome_name = paste0('chr', coords$chromosome_name)
coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)

num = which(is.na(coords$chromosome_name))
coords = coords[-num, ]

up = which((mat$pval < 0.01) &
             (mat$log2FC > 1))
upmat = mat[up, ]

num = which(coords$hgnc_symbol %in% rownames(upmat))
coords1 = coords[num, ]

RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)
genes = intersect(rownames(mat), coords$hgnc_symbol)

mat1 = mat[genes, ]
df = cbind.data.frame(rownames(mat1), mat1[, c(1,2,4)])
colnames(df)[1] = 'hgnc_symbol'

data = merge(coords, df, by = 'hgnc_symbol', all.x = T)
data=data[,c(2,3,4,1,5,6,7)]

data = data[, c('chromosome_name', 'start_position',
                'end_position', 'hgnc_symbol',
                'meanTumor', 'meanControl', 'log2FC')]
RCircos.Heatmap.Plot(data, data.col = 7, track.num = 6, side = "in",
                     min.value = -2, max.value =2, genomic.columns = 3,
                     is.sorted = F)
RC.param = RCircos.Get.Plot.Parameters()
RC.param['heatmap.color'] = "GreenWhiteRed"
RCircos.Reset.Plot.Parameters(RC.param)

RCircos.Heatmap.Plot(data, data.col = 7, track.num = 10, side = "in",
                     min.value = -2, max.value = 2,
                     is.sorted = F)
