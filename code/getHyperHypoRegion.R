#---get hyper hypo regions----
#---20200312------------------

options(scipen = 300);

args = commandArgs(TRUE);
Out = as.character(args[1]);
FDR = as.numeric(args[2]);
Resultfile = as.character(args[3]);


Top = 1000;

#Out = '/home/zmwang/SomaticMutationPrediction/LUSC-KR/Falcon_chromHMM/MergeModel_ChromHMM_apply'
setwd(Out);

#FDR = 0.01
#Resultfile = 'Result.txt'


Result = read.table(Resultfile,header=T,sep='\t',stringsAsFactors =F);
str(Result)

#Result[,'hyper_fdr'] = p.adjust(Result[,'hyper_pvalue'],method='fdr')
#Result[,'hypo_fdr'] = p.adjust(Result[,'hypo_pvalue'],method='fdr')

Index1 = which(Result[,'hyper_fdr']<FDR);
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
	write.table(tmp,file=paste0('Hyper_regions_FDR',FDR,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');
}




Index2 = which(Result[,'hypo_fdr']<FDR);
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
	write.table(tmp,file=paste0('Hypo_regions_FDR',FDR,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');
}





Hyper = Result[Index1,c("NormalizedMutationRate","PredictedNormalizedRated")]

Hypo = Result[Index2,c("NormalizedMutationRate","PredictedNormalizedRated")]

Random = Result[-c(Index1,Index2),c("NormalizedMutationRate","PredictedNormalizedRated")]

MIN = min(c(Result[,"NormalizedMutationRate"],Result[,'PredictedNormalizedRated']));
MAX = max(c(Result[,"NormalizedMutationRate"],Result[,'PredictedNormalizedRated']));



pdf("Scatterplot_hyper_hypo.pdf");
plot(x=Random[,'PredictedNormalizedRated'],y=Random[,'NormalizedMutationRate'],col='#bababa',xlim=c(MIN,MAX),ylim=c(MIN,MAX),
	panel.first = grid(lty = 1, lwd = 2,col='white'),pch=19,type='n',xlab='Prediction',ylab='Observation');
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f7f7f7")
panel.first = grid(lty = 1, lwd = 2,col='white')
points(x=Random[,'PredictedNormalizedRated'],y=Random[,'NormalizedMutationRate'],col='#bababa',pch=20);
points(x=Hyper[,'PredictedNormalizedRated'],y=Hyper[,'NormalizedMutationRate'],col='red',pch=20);
points(x=Hypo[,'PredictedNormalizedRated'],y=Hypo[,'NormalizedMutationRate'],col='blue',pch=20);
lines(x=seq(MIN,MAX,0.1),y=seq(MIN,MAX,0.1));
dev.off();



Result = Result[order(Result[,'hyper_fdr'],decreasing=F),]

Index1 = seq(1,Top)
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
	write.table(tmp,file=paste0('Hyper_regions_Top',Top,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');
}


Result = Result[order(Result[,'hypo_fdr'],decreasing=F),]

Index2 = seq(1,Top)
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
	write.table(tmp,file=paste0('Hypo_regions_Top',Top,'.bed'),col.names=F,row.names=F,quote=F,sep='\t');
}



