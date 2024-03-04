setwd('/Users/tshiro/Desktop/Projects/HGEOC/Github/JGOG_HGEOC/scripts')

#https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html#introduction

library(annotatr)
library(karyoploteR)
annots = c('hg19_cpgs', 'hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annots)

# EM >FT
in_bed='../data/EM_FT_UP.bed'
out_f='../data/EM_FT_UP_promotor_genes.txt'
gr=read_regions(in_bed)
annotated = annotate_regions(
  regions = gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
df_annotated = data.frame(annotated)
df_promotor=df_annotated[grep(df_annotated$annot.id,pattern = 'promoter'),]
dim(df_promotor)
Hugo_symbol=sort(unique(df_promotor$annot.symbol))
length(Hugo_symbol)
dfw=data.frame(Hugo_symbol)
write.table(dfw,out_f,quote = F,col.names = F,row.names = F)

# FT >EM
in_bed='../data/FT_EM_UP.bed'
out_f='../data/FT_EM_UP_promotor_genes.txt'
gr2=read_regions(in_bed)
annotated = annotate_regions(
  regions = gr2,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
df_annotated = data.frame(annotated)
df_promotor=df_annotated[grep(df_annotated$annot.id,pattern = 'promoter'),]
dim(df_promotor)
Hugo_symbol=sort(unique(df_promotor$annot.symbol))
length(Hugo_symbol)
dfw=data.frame(Hugo_symbol)
write.table(dfw,out_f,quote = F,col.names = F,row.names = F)

# Plot
kp <- plotKaryotype(genome="hg19",chromosomes = paste0('chr', c(c(1:22),"X")))
kpPlotRegions(kp, data=gr, r0=0.05, r1=0.15, col='red')
kpPlotRegions(kp, data=gr2, r0=0.2, r1=0.3, col='blue')

