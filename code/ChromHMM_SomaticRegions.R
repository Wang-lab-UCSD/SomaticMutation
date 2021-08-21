#---feature distribution of somatic regions--------------
#---20190228---------------------------------------------
#---wangzengmiao@gmail.com-------------------------------
options(scipen = 300);



CHR = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11',
	'chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX');



args = commandArgs(TRUE);
Out = as.character(args[1]);

setwd(Out);



SNV = read.table("ChromHMM_SomaticMutation.bed",header=F,sep='\t',stringsAsFactors=F)

ID = paste0(SNV[,1],":",SNV[,2],"-",SNV[,3]);
Freq = as.data.frame(table(ID));
Freq[,1] = as.character(as.matrix(Freq[,1]));

ID = Freq[,'ID']
ID = strsplit(ID,":");
ID = do.call('rbind',ID);

tmp = strsplit(ID[,2],"-");
tmp = do.call('rbind',tmp);


Region = data.frame(chr=ID[,1],s=tmp[,1],e=tmp[,2],ID=Freq[,1],Num_mutation=Freq[,2]);
Region[,c(1,4)] = as.character(as.matrix(Region[,c(1,4)]));
Region[,c(2,3)] = as.numeric(as.matrix(Region[,c(2,3)]));
Region = Region[order(Region[,1],Region[,2]),]

write.table(Region,file=paste0("SomaticRegion_Motif_hg19",".bed"),quote=F,sep="\t",col.names=F,row.names=F);

L = Region[,3]-Region[,2];
print("L range is")
range(L);
print("L mean is");
mean(L);
print("L median is");
median(L);
print("L dimension is");
dim(Region)
print("number of mutation is");
sum(Region[,5])

pdf(paste0("Hist_of_L_SomaticRegion_hg19",".pdf"));
hist(L, breaks=50,main="Hist for the length of regions");
legend("topright",legend=c(paste0('range(L): ',round(min(L),2),"-",round(max(L),2)),paste0("mean: ",round(mean(L),2)),
	paste0("median: ",round(median(L),2)),paste0("Num of Regions:",dim(Region)[1]),paste0("Num of mutation:",sum(Region[,5]))));
dev.off();





pdf("MutationCount_L.pdf");
Cor = round(cor(log2(L),log2(Region[,5])),5);
plot(x=log2(L),y=log2(Region[,5]),xlab='log2(length)',ylab='log2(MutationCount)',main=paste0("cor:",Cor))
dev.off()


Rate = Region[,5]/L*1000/sum(Region[,5])*1000000
Rate = log2(Rate+1)
pdf("Rate_L.pdf")
plot(x=log2(L),y=Rate);
dev.off();


pdf("Hist_rate.pdf")
hist(Rate,breaks=100)
dev.off()













