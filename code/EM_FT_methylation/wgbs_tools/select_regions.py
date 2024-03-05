import pandas as pd
import sys

dmr1=pd.read_csv('EM3FT3_hyper/Markers.EM.bed',sep='\t')
dmr2=pd.read_csv('EM3FT3_hyper/Markers.FT.bed',sep='\t')
#print(dmr1.shape, dmr2.shape)

dmr1=dmr1[(dmr1['ttest']<0.00005)&(dmr1['delta_quants']>0.3)]
dmr2=dmr2[(dmr2['ttest']<0.00005)&(dmr2['delta_quants']>0.3)]
print(dmr1.shape, dmr2.shape)

dmr1.iloc[:,0:3].to_csv('EM_FT_UP.bed',sep='\t',index=False,header=False)
dmr2.iloc[:,0:3].to_csv('FT_EM_UP.bed',sep='\t',index=False,header=False)

