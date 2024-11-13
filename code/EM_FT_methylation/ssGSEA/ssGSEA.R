setwd("/Users/tshiro/Desktop/Projects/HGEOC/Github/JGOG_HGEOC/code/EM_FT_methylation/ssGSEA")
source("~/Dropbox/R/Software/ssGSEA_original/common.R")
source("~/Dropbox/R/Software/ssGSEA_original/ssGSEAProjection.Library_thres_max4010.R") 
# input
gct_f = '/Users/tshiro/Desktop/Projects/JGOG3025_HGS/JGOG3025_data/RNAseq/gene_expression_add/JGOG_HGS_RSEM_tpm_gs_merged_renamed.gct'
gct_f = "/Users/tshiro/Desktop/DRY/TCGA_analysis/OV/GDC2023/RNASeq_STAR/TCGA_OV_421_tpm_gs.gct"
gct_f = '/Users/tshiro/Desktop/DRY/TCGA_analysis/UCEC/GDC2023/gexp/TCGA_UCEC_557_tpm_gs_merged.gct'

## from GEO
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets/GSE2109_OV264_GPL570/GSE2109_ov_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets/GSE6008_OV103_GPL96/GSE6008_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets/GSE19539_OV68_HuGene-1_0-st/GSE19539-GPL6244_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets/GSE44104_OV60_GPL570/GSE44104_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets/GSE65986_TOV55_GPL570/GSE65986_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/Projects/Bev_HRD/datasets/NatCom2020_ICON7/expression/gexp_matrix/ICON7_FFPE_369_tpm_gs_merged.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets_EC/GSE2109_EC225_GPL570/GSE2109_ucec_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets_EC/GSE23518_EC_SEC/GSE23518_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets_EC/GSE24537_EC_SEC/GSE24537_exp_gs_max.gct'
gct_f='/Users/tshiro/Desktop/DRY/GEO_datasets_EC/GSE56026_EC_SEC/GSE56026_exp_gs_max.gct'

gmt_fs=list.files(".",pattern = '.gmt',full.names = T)
gmt_fs
gmt_f=gmt_fs

# output
out_f = "./JGOG_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./TCGA_OV_421_tpm_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./TCGA_UCEC_557_tpm_EM_FT_promotor_genes_ssgsea_score.txt" 

out_f = "./GSE2109_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE6008_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE19539_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE44104_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE65986_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./NatCom2020_ICON7_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE2109_EC_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE23518_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE24537_EM_FT_promotor_genes_ssgsea_score.txt" 
out_f = "./GSE56026_EM_FT_promotor_genes_ssgsea_score.txt" 

# zscaled
out_f2 = paste0(sub(pattern = ".txt",replacement = "", out_f),"_z.txt")

# ssGSEA
ssGSEA <- ssGSEA.project.dataset(
  input.ds= gct_f, 
  output.ds= "TMP.gct",
  gene.sets.database= "",
  gene.sets.dbfile.list = gmt_f,
  gene.symbol.column= "Name",
  gene.set.selection  = "ALL",
  sample.norm.type    = "rank",
  weight              = 0.75,
  combine.mode        = "combine.add",
  min.overlap         = 10)


y <- read.delim("TMP.gct", header=F, comment.char="#", na.strings="")
row_line=y[,1]
# if duplicated rows are exists, retain only the first one
if (sum (duplicated(row_line)) != 0){
  y=y[duplicated(row_line)==FALSE,]
}
rowname=y[-c(1,2),1]
colname= apply(y[2,-c(1,2)],MARGIN = 2, as.character)
dm_y=y[-c(1,2),-c(1,2)]
dm_y= apply(dm_y,MARGIN = c(1,2),as.numeric)
dm_y_w=round( dm_y )
row.names(dm_y_w)=rowname
colnames(dm_y_w)=colname
write.table(dm_y_w, out_f, append = F, quote = F, sep = "\t", eol = "\n",
            na = "", dec = "." , row.names = T, col.names = NA)

# zscale
z=function(x){
  a= (x - mean(x) )/ sd(x)
  return (a)
}
dm_y_z=apply(dm_y,MARGIN = 1,z )
dm_z=(t(dm_y_z))
dm_z=round(dm_z,3)
row.names(dm_z)=rowname
colnames(dm_z)=colname
write.table(dm_z, out_f2, append = F, quote = F, sep = "\t",
            eol = "\n", na = "", dec = "." ,
            row.names = T, col.names = NA)

