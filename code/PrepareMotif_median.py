#!/usr/bin/env python
# this script is used to generate sample matrix for Motif for the interested region.
# Wang, Zengmiao May 20,2019
# Note that I use the gene list for HCT116, it includes all genes.
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import gc 

Chr = sys.argv[1] 
RegionFile = sys.argv[2]
Out = sys.argv[3]
FiMOFolder = sys.argv[4]
FimoPvalue = float(sys.argv[5])

hg19 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
 'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
#--GeneList is the protein list who has motifs from CIS-BP or from factorbook.
#--I download motifs for some proteins from factorbook in GM12878 and K562 cell line.
#--GeneGM12878 is the motif list that is specific to GM12878.
GeneList = pd.read_csv(FiMOFolder+'/GeneListForHCT116.txt',sep="\t",header=None)
GeneList = GeneList[0].tolist()
TF2Motif = pd.read_csv(FiMOFolder+"/TF2Motif.txt",sep="\t")
TF2Motif.index = TF2Motif['Motif_ID']
TFName = TF2Motif['TF_Name'].unique()
MotifList = list(TF2Motif['Motif_ID'].unique()) + list(set(GeneList) - set(TFName))
GeneGM12878 = ["BCL3","WRNIP1","NFIC","MTA3","KAT2A","SUPT20H"]

Motif2Protein = TF2Motif[['Motif_ID','TF_Name']]
tmp = pd.DataFrame(data=np.zeros((len(list(set(GeneList) - set(TFName))),2),dtype = np.string_ ),columns=['Motif_ID','TF_Name'])
tmp['Motif_ID'] = list(set(GeneList) - set(TFName))
tmp['TF_Name'] = list(set(GeneList) - set(TFName))
Motif2Protein = tmp.append(Motif2Protein, ignore_index=True)

#--histon motif------
tmp = pd.read_csv(FiMOFolder+'/HistoneMotifID.txt',sep='\t',header=None)
MotifList = list(MotifList+tmp[1].tolist())
tmp1 = pd.DataFrame(data=np.zeros((len(tmp[1].tolist()),2),dtype = np.string_ ),columns=['Motif_ID','TF_Name'])
tmp1['Motif_ID'] = tmp[1].tolist()
tmp1['TF_Name'] = tmp[1].tolist()
Motif2Protein = tmp1.append(Motif2Protein, ignore_index=True)

#--methylation motif----
tmp = pd.read_csv(FiMOFolder+'/methylationMotif/PublishVersion313Motif/MethylationMotifID.txt',sep=' ',header=None)
MotifList = list(MotifList+tmp[1].tolist())
tmp1 = pd.DataFrame(data=np.zeros((len(tmp[1].tolist()),2),dtype = np.string_ ),columns=['Motif_ID','TF_Name'])
tmp1['Motif_ID'] = tmp[1].tolist()
tmp1['TF_Name'] = tmp[1].tolist()
Motif2Protein = tmp1.append(Motif2Protein, ignore_index=True)



Background = pd.read_csv(FiMOFolder+"/Summary_CIS_BP_Factorbook_Histone_Methy_Motif.txt",sep='\t');

MotifList = Background['Motif_ID'].tolist()



Region = pd.read_csv(RegionFile,sep='\t',header=None)
Region = Region[Region[0]==Chr]
Region[3] = Region[0]+':'+Region[1].map(str)+'-'+Region[2].map(str)
Region.index = Region[3]

Sample = pd.DataFrame(data=np.zeros((len(Region),len(MotifList)+2)),columns=list(['Region_ID','Region_length']+MotifList),index=Region[3])


#---Loading Motif -----------------------------------
#---TF motif----
FIMO_TF = pd.read_csv(FiMOFolder+"/"+Chr+".fa_fimo_TF.txt/fimo.txt",sep="\t")
print 'I have read '+Chr+' Chromosome FIMO file'
FIMO_TF = FIMO_TF[FIMO_TF['p-value']<FimoPvalue]
Index = FIMO_TF['# motif_id'].isin(TF2Motif.index)
FIMO_TF = FIMO_TF[Index]
#---Motifs from Factorbook------------
tmp = pd.read_csv(FiMOFolder+"/"+Chr+".fa_fimo_Remodeller.txt/fimo.txt",sep="\t")
tmp = tmp[tmp['p-value']<FimoPvalue]
Index = tmp['# motif_id'].isin(GeneGM12878)
tmp = tmp[Index]
FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
tmp = pd.read_csv(FiMOFolder+"/"+Chr+".fa_fimo_RemodellerForK562.txt/fimo.txt",sep="\t")
tmp = tmp[tmp['p-value']<FimoPvalue]
FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
#---histone motif----
tmp = pd.read_csv(FiMOFolder+"/"+Chr+".fa_fimo_HistoneModification.txt/fimo.txt",sep="\t")
tmp = tmp[tmp['p-value']<FimoPvalue]
FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
#---Methylation Motif-----
tmp = pd.read_csv(FiMOFolder+"/methylationMotif/PublishVersion313Motif/"+Chr+".fa_fimo_MethylationMotif313/fimo.txt",sep="\t")
tmp = tmp[tmp['p-value']<FimoPvalue]
FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)


FIMO_TF.loc[:,'start'] = FIMO_TF['start']-1
FIMO_TF.loc[:,'stop'] = FIMO_TF['stop']-1


Index = FIMO_TF['# motif_id'].isin(MotifList)
FIMO_TF = FIMO_TF[Index]

#--Loading Motif end------------------------------
for i in Sample.index:
	tmp = FIMO_TF[(FIMO_TF['start']>=Region.loc[i,1]) & (FIMO_TF['stop']<=Region.loc[i,2])]
	Sample.loc[i,'Region_length'] = Region.loc[i,2]-Region.loc[i,1]+1
	if len(tmp)>0:
		ME = tmp.groupby('# motif_id')['p-value'].median()
		ME = -np.log10(ME)
		Sample.loc[i,ME.index] = ME[0:len(ME)]
	print(i)


Sample['Region_ID'] = Region[3]

Sample.to_csv(Out+'/Motif_MinusLog10pvalueMedian_'+Chr+'.txt',sep="\t",header=True,index=False)

# Background.index = Background['Motif_ID']

# SampleNorm = Sample
# for motif in MotifList:
# 	L = SampleNorm['Region_length']-Background.loc[motif,'Motif_length']+1
# 	Index = L<1
# 	if sum(Index)>0:
# 		L[Index] =1
# 	N = Background.loc[motif,'TotalNum']/100000.0
# 	SampleNorm[motif] = SampleNorm[motif]/L/N
# 	print(motif)

# SampleNorm.to_csv(Out+'/Motif_NormalizedCount_'+Chr+'.txt',sep="\t",header=True,index=False)

del FIMO_TF,tmp
gc.collect()




