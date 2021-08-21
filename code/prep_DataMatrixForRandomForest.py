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
import os
from random import randint

SNVfile = sys.argv[1] 
RegionFile = sys.argv[2]
Out = sys.argv[3]
OutName = sys.argv[4]
FiMOFolder = sys.argv[5]
h = float(sys.argv[6])
BedtoolFolder = sys.argv[7]
FimoPvalue = float(sys.argv[8])


Value = randint(0, 100000)

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


#--Novel motif--------
tmp = pd.read_csv(FiMOFolder+'/MotifFinding/NovelMotif2Chip.txt',sep='\t')
tmp['Motif'] = [x.split(' ')[1] for x in tmp['Motif']]

MotifList = list(MotifList+tmp['Motif'].tolist())

tmp1 = pd.DataFrame(data=np.zeros((len(tmp),2),dtype = np.string_ ),columns=['Motif_ID','TF_Name'])
tmp1['Motif_ID'] = tmp['Motif'].tolist()
tmp1['TF_Name'] = tmp['Protein'].tolist()
Motif2Protein = tmp1.append(Motif2Protein, ignore_index=True)

if h==0:
	Background = pd.read_csv(FiMOFolder+"/TF_Histone_Methy_Novel_Motif_summary.txt",sep='\t')
	MotifList = Background['Motif_ID'].tolist()
else:
	Background = pd.read_csv(FiMOFolder+"/Cluster2431Motif_Info_h"+str(h)+".txt",sep='\t')
	MotifList = Background['ID'].tolist()




SNV = pd.read_csv(SNVfile,sep='\t',header=None)
DonorID = SNV[3].unique()

Region = pd.read_csv(RegionFile,sep='\t',header=None)
Region[3] = Region[0]+':'+Region[1].map(str)+'-'+Region[2].map(str)
CHR = Region[0].unique()

Sample = pd.DataFrame(data=np.zeros((len(DonorID),len(MotifList)+len(Region)+1)),columns=list(['DonorID']+Region[3].tolist()+MotifList),index=DonorID)



# Region = Region[Region[0]==Chr]
# Region[3] = Region[0]+':'+Region[1].map(str)+'-'+Region[2].map(str)
# Region.index = Region[3]

for Chr in CHR:
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
	#---Novel Motif----------
	tmp = pd.read_csv(FiMOFolder+"/MotifFinding/"+Chr+".fa_fimo_NovelMotif_batch1/fimo.txt",sep="\t")
	tmp = tmp[tmp['p-value']<FimoPvalue]
	FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
	tmp = pd.read_csv(FiMOFolder+"/MotifFinding/"+Chr+".fa_fimo_NovelMotif_batch2/fimo.txt",sep="\t")
	tmp = tmp[tmp['p-value']<FimoPvalue]
	FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
	tmp = pd.read_csv(FiMOFolder+"/MotifFinding/"+Chr+".fa_fimo_NovelMotif_batch3/fimo.txt",sep="\t")
	tmp = tmp[tmp['p-value']<FimoPvalue]
	FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
	tmp = pd.read_csv(FiMOFolder+"/MotifFinding/"+Chr+".fa_fimo_NovelMotif_batch4/fimo.txt",sep="\t")
	tmp = tmp[tmp['p-value']<FimoPvalue]
	FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
	tmp = pd.read_csv(FiMOFolder+"/MotifFinding/"+Chr+".fa_fimo_NovelMotif_batch5/fimo.txt",sep="\t")
	tmp = tmp[tmp['p-value']<FimoPvalue]
	FIMO_TF = tmp.append(FIMO_TF, ignore_index=True)
	FIMO_TF.loc[:,'start'] = FIMO_TF['start']-1
	FIMO_TF.loc[:,'stop'] = FIMO_TF['stop']-1
	Index = FIMO_TF['# motif_id'].isin(MotifList)
	FIMO_TF = FIMO_TF[Index]
	FIMO_TF[['sequence_name','start','stop','# motif_id']].to_csv(Out+'/MotifBind_'+Chr+'_'+str(Value)+'.bed',sep="\t",header=False,index=False)
	Commond = 'sort -k1,1 -k2,2n '+'MotifBind_'+Chr+'_'+str(Value)+'.bed > '+Chr+'_'+str(Value)+'.sorted.bed'
	os.system(Commond)
	Commond = BedtoolFolder+'/bedtools intersect -a '+Chr+'_'+str(Value)+'.sorted.bed -b '+SNVfile+' -wo > '+'Motif_SNV_intersect_'+Chr+'_'+str(Value)+'.bed'
	os.system(Commond)
	Region_tmp = Region[Region[0]==Chr]
	Region_tmp.index = range(0,len(Region_tmp))
	SNV_tmp = SNV[SNV[0]==Chr]
	for i in Region_tmp.index:
		tmp = SNV_tmp[(SNV_tmp[1]>=Region_tmp.loc[i,1]) & (SNV_tmp[2]<=Region_tmp.loc[i,2])]
		Count = tmp[3].value_counts()
		Sample.loc[Count.index,Region_tmp.loc[i,3]] = Count
	MotifSNV = pd.read_csv(Out+"/Motif_SNV_intersect_"+Chr+'_'+str(Value)+".bed",sep="\t",header=None)
	Motif_tmpID = MotifSNV[3].unique()
	Donor_tmpID = MotifSNV[7].unique()
	for i in Donor_tmpID:
		tmp = MotifSNV[MotifSNV[7]==i]
		Count = tmp[3].value_counts()
		Sample.loc[i,Count.index] = Count
	del FIMO_TF,Region_tmp,SNV_tmp,MotifSNV,tmp
	gc.collect()
	print(Chr)
	Commond = 'rm '+'MotifBind_'+Chr+'_'+str(Value)+'.bed'
	os.system(Commond)
	Commond = 'rm '+Chr+'_'+str(Value)+'.sorted.bed'
	os.system(Commond)
	Commond = 'rm '+'Motif_SNV_intersect_'+Chr+'_'+str(Value)+'.bed'
	os.system(Commond)



Sample['DonorID'] = DonorID
Sample.to_csv(Out+'/'+OutName,sep="\t",header=True,index=False)



