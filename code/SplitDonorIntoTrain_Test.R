#----split samples into train and test------
#----2020-03-11------------------

options(scipen = 300);
args = commandArgs(TRUE);
SNVfile = as.character(args[1]);
Out = as.character(args[2]);

setwd(Out);


SNV = read.table(SNVfile,header=F,sep='\t',stringsAsFactors=F);

Donor = unique(SNV[,4]);


print("The total donor size is ");
length(Donor)


TrainDonor = sample(Donor,length(Donor)*0.7);
print('The donor size in train set is ');
length(TrainDonor);

TestDonor = setdiff(Donor,TrainDonor);
print('The donor size in test set is ');
length(TestDonor);


Index = match(SNV[,4],TrainDonor,nomatch=-1);

Index = which(Index>0);
Train = SNV[Index,];
print('The mutation number in train set is ');
dim(Train)[1];



Test = SNV[-Index,]

print('The mutation number in test set is ');
dim(Test)[1];


write.table(Train,file='SomaticMutation_Train.bed',col.names=F,sep='\t',row.names=F,quote=F)

write.table(Test,file='SomaticMutation_Test.bed',col.names=F,sep='\t',row.names=F,quote=F)



















