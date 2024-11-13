# combat_seq
setwd('/Users/tshiro/Desktop/Projects/HGEOC/Github/JGOG_HGEOC_20240404_dev/code/Machine_learning')

library(sva)
library(data.table)

##############################################
### TCGA-OV + JGOG
##############################################
in_df1=fread('./TCGA-OV_555_ascat_spmg.CNV48.matrix.tsv'
             ,data.table = F)
row.names(in_df1)=in_df1[,1]
in_df1=in_df1[,-1]
in_df1=t(in_df1)

in_df2=fread('../../data/JGOG_282_ascat_spmg.CNV48.matrix.txt'
             ,data.table = F)
row.names(in_df2)=in_df2[,1]
in_df2=in_df2[,-1]
in_df2=t(in_df2)

out_df1="TCGA_OV_555_ascat_spmg.CNV48.matrix.combat_to_JGOG.tsv"
out_df2="JGOG3025_282_ascat_spmg.CNV48.matrix.combat_to_TCGA-OV.tsv"

count_matrix=as.matrix( cbind(in_df1,in_df2) )

batch <- c(rep(1, dim(in_df1)[2]), rep(2, dim(in_df2)[2]))
adjusted <- ComBat_seq( count_matrix, batch=batch,shrink = T)

dfw1=t(adjusted[,colnames(in_df1)])
dfw2=t(adjusted[,colnames(in_df2)])

write.table(dfw1,out_df1,sep='\t',quote = F,row.names = T,col.names = NA)
write.table(dfw2,out_df2,sep='\t',quote = F,row.names = T,col.names = NA)


##############################################
### TCGA-OV+UCEC + JGOG
##############################################
in_df3=fread('./TCGA-UCEC_287_ascat_spmg.CNV48.matrix.tsv'
             ,data.table = F)
row.names(in_df3)=in_df3[,1]
in_df3=in_df3[,-1]
in_df3=t(in_df3)

in_df13=cbind(in_df1,in_df3)
dim(in_df13)

out_df3="./TCGA-OV-UCEC_842_ascat_spmg_CNV48.matrix.combat_to_JGOG.tsv"
out_df4="./JGOG_282_ascat_spmg.CNV48.matrix.combat_to_TCGA-OV-UCEC.tsv"

count_matrix=as.matrix( cbind(in_df13,in_df2) )

batch <- c(rep(1, dim(in_df13)[2]), rep(2, dim(in_df2)[2]))
adjusted <- ComBat_seq( count_matrix, batch=batch)

dfw3=t(adjusted[,colnames(in_df13)])
dfw4=t(adjusted[,colnames(in_df2)])

write.table(dfw3,out_df3,sep='\t',quote = F,row.names = T,col.names = NA)
write.table(dfw4,out_df4,sep='\t',quote = F,row.names = T,col.names = NA)
