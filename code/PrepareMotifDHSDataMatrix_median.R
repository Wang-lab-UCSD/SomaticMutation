#--prepare motif+DHS data matrix---
#--wangzengmiao@gmail.com----------
#--20190522------------------------

options(scipen = 300);
CHR = paste0("chr",seq(1,22));
CHR = c(CHR,"chrX");
args = commandArgs(TRUE);
Dir = as.character(args[1]);
TotalCount = as.numeric(args[2]);


setwd(paste0(Dir,"MotifDHS_hg19"));

DHS = read.table("OpenChromatin_hg19_RPKM.bed",header=F,sep='\t');
DHS[,1] = as.character(as.matrix(DHS[,1]));
L = DHS[,3]-DHS[,2]+1;
DHS[,6] = DHS[,6]/(L/1000*TotalCount/1000000);
ID = paste0(DHS[,1],":",DHS[,2],"-",DHS[,3]);
rownames(DHS) = ID;

i = 'chr1';
MotifCount = read.table(paste0("Motif_MinusLog10pvalueMedian_",i,'.txt'),header=T,sep='\t');


for(i in 2:length(CHR)){
	tmp = read.table(paste0("Motif_MinusLog10pvalueMedian_",CHR[i],'.txt'),header=T,sep='\t');
	MotifCount = rbind(MotifCount,tmp);
	print(i);
}

MotifCount[,"Region_ID"] = as.character(as.matrix(MotifCount[,"Region_ID"]));
dim(MotifCount);
dim(DHS);
Index = match(MotifCount[,"Region_ID"],rownames(DHS));
DHS = DHS[Index,];
sum(rownames(DHS)==MotifCount[,"Region_ID"]);
colnames(MotifCount)[1:2] = c("ID","L");
MotifCount = cbind(MotifCount,DHS=DHS[,6]);


Node = read.table(paste0(Dir,"SomaticRegion_Motif_hg19.bed"),header=F,sep='\t',stringsAsFactors=F);
ID = paste0(Node[,1],":",Node[,2],"-",Node[,3]);
rownames(Node) = ID;
Index = match(MotifCount[,"ID"],ID);
Node = Node[Index,];
sum(rownames(Node)==MotifCount[,"ID"])

MotifCount = cbind(MotifCount[,1:2],Num_mutation=Node[,5],MotifCount[,3:dim(MotifCount)[2]]);


write.table(MotifCount,file="MotifMedian_DHS_Data.txt",row.names=F,col.names=T,sep='\t',quote=F);

SUM = apply(MotifCount[,4:dim(MotifCount)[2]],1,sum);
Index = which(SUM==0);
length(Index);

if(length(Index)>0){


tmp = MotifCount[-Index,];
rm(MotifCount)
gc();
write.table(tmp,file="MotifMedian_DHS_Data_FilterZero.txt",row.names=F,col.names=T,sep='\t',quote=F);
}




