#---get hyper hypo regions----
#---20200312------------------

options(scipen = 300);

args = commandArgs(TRUE);
Out = as.character(args[1]);
Resultfile = as.character(args[2]);
DataSet = as.character(args[3]);


Top = c(100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000);

#Out = '/home/zmwang/SomaticMutationPrediction/LUSC-KR/Falcon_chromHMM/MergeModel_ChromHMM_apply'
setwd(Out);

#FDR = 0.01
#Resultfile = 'Result.txt'


Result = read.table(Resultfile,header=T,sep='\t',stringsAsFactors =F);
str(Result)

Result[,'hyper_fdr'] = p.adjust(Result[,'hyper_pvalue'],method='fdr')
Result[,'hypo_fdr'] = p.adjust(Result[,'hypo_pvalue'],method='fdr')


FDR = 0.05
Index = which(Result[,'hyper_fdr']<FDR)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hyper_regions_FDR',FDR,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}

FDR = 0.05
Index = which(Result[,'hypo_fdr']<FDR)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hypo_regions_FDR',FDR,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}




FDR = 0.01
Index = which(Result[,'hyper_fdr']<FDR)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hyper_regions_FDR',FDR,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}



FDR = 0.01
Index = which(Result[,'hypo_fdr']<FDR)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hypo_regions_FDR',FDR,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}



Pvalue = 0.05
Index = which(Result[,'hyper_pvalue']<Pvalue)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hyper_regions_pvalue',Pvalue,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}

Pvalue = 0.05
Index = which(Result[,'hypo_pvalue']<Pvalue)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hypo_regions_pvalue',Pvalue,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}

Pvalue = 0.01
Index = which(Result[,'hyper_pvalue']<Pvalue)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hyper_regions_pvalue',Pvalue,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}

Pvalue = 0.01
Index = which(Result[,'hypo_pvalue']<Pvalue)
if(length(Index)>0){
	tmp = Result[Index,'ID'];
	tmp = strsplit(tmp,":");
	tmp = do.call("rbind",tmp);
	Start = tmp[,2]
	Start = strsplit(Start,"-");
	Start = do.call("rbind",Start);
	tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index,'ID']);
	tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
	tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
	tmp = tmp[order(tmp[,1],tmp[,2]),];
	write.table(tmp,file=paste0(DataSet,'_Hypo_regions_pvalue',Pvalue,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');	
}







Result = Result[order(Result[,'hyper_fdr'],decreasing=F),]

for(j in 1:length(Top)){
	Index1 = seq(1,Top[j])
	if(length(Index1)>0){
		tmp = Result[Index1,'ID'];
		tmp = strsplit(tmp,":");
		tmp = do.call("rbind",tmp);
		Start = tmp[,2]
		Start = strsplit(Start,"-");
		Start = do.call("rbind",Start);
		tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index1,'ID']);
		tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
		tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
		tmp = tmp[order(tmp[,1],tmp[,2]),];
		write.table(tmp,file=paste0(DataSet,'_Hyper_regions_Top',Top[j],'.bed'),col.names=F,row.names=F,quote=F,sep='\t');
	}

}



Result = Result[order(Result[,'hypo_fdr'],decreasing=F),]

for(j in 1:length(Top)){
	Index2 = seq(1,Top[j])
	if(length(Index2)>0){
		tmp = Result[Index2,'ID'];
		tmp = strsplit(tmp,":");
		tmp = do.call("rbind",tmp);
		Start = tmp[,2]
		Start = strsplit(Start,"-");
		Start = do.call("rbind",Start);
		tmp = data.frame(chr=tmp[,1],s=Start[,1],e=Start[,2],ID=Result[Index2,'ID']);
		tmp[,c(1,4)] = as.character(as.matrix(tmp[,c(1,4)]));
		tmp[,c(2,3)] = as.numeric(as.matrix(tmp[,c(2,3)]));
		tmp = tmp[order(tmp[,1],tmp[,2]),];
		write.table(tmp,file=paste0(DataSet,'_Hypo_regions_Top',Top[j],'.bed'),col.names=F,row.names=F,quote=F,sep='\t');
	}


}



