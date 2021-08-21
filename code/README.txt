#---------------------------Introduction--------------------------------

Note that only somatic mutation subsitition is focused here.

#------------------------------Input------------------------------------

Two input files are required. 

1. The somatic mutation file from a set of donors. The file should be
   in bed format without header. Each row represents a somatic mutation.
   The first column is the chromosome number(for example, chr1,chr2,...).
   The second column is the start point and the third column is the end
   point. The fourth column is the donorID for this mutation.These columns
   are seperated by tab. If multiple donors have somatic mutation in the 
   same loci, please list it multiple times. For example, three donors 
   have a somatic mutation in chr2:10000, the lines will be: 
   chr2 \tab 10000 \tab 10000 \tab DonorID
   chr2 \tab 10000 \tab 10000 \tab DonorID
   chr2 \tab 10000 \tab 10000 \tab DonorID

   Note that If Split=YES(see how to use section), DonorID is required. 
   Otherwise, DonorID column is option.

2. The open chromatin data in the bam format. The unmapped and duplicated
   multi-aligned reads should be removed from bam file. 

#------------------------Software and Packages---------------------------

1. R. Required packages are: pracma and mixtools 

2. Python. Required packages are: pandas, matplotlib, numpy, scipy, keras,
   sklearn, tensorflow,statsmodels.

3. Bedtools. Please use the versions earlier than 2.19.0. Suggested version
   is 2.18.0. Since there are some functions removed from the latest version.
   I will re-write this part of code to remove the version limitation.

4. Samtools. 


#-----------------------------How to use-----------------------------------------

For each run, modify the Argument_File.txt file and then run the following command 
on the head node:

nohup ./falcon.job Argument_File.txt > run_log.out &

falcon.job is executable file and Argument_File.txt is the parameter file. run_log.out
is the log file for this run. You can give parameter file and log file different names 
for different runs. This will help you track different runs.


Open Argument_File.txt file and modify "Input folders and files", "job parameter"
and "Software bin folder" section. Note that there is NO space before and after "=".

1. "Input folders and files"

   a. ScriptFolder. The absolute path of scripts.

   b. SNVFile. The somatic mutation file with absolute path.

   c. OpenChromFile. The open chromatin file with absolute path.

   d. RegionFile. The chromHMM region file with absolute path for corresponding tissue.

   e. Split. whether or not split the the dataset into train and test. If YES, the ratio
      for train versus test is 7:3. The train data will be used to train CR model.

   f. OutDir. The output directory.

   g. FimoPvalue. The p-value cutoff for FIMO. By default, it is 1e-05.

   h. FiMOFolder. The Motif scanncing folder. I have uploaded the data to the share
      folder. No need to change it.

2. "job parameter"

   a. SLEEP. This is used to control when the job will start. If you have only 1 dataset 
      to run, you can set to 1s. It means that the job will start in 1s. If you have two
      datasets to run, you can set SLEEP to 1s in the first parameter file and set SLEEP
      to 12h for the second parameter files. This means that the job for the first dataset
      start immediately, and the job for the second dataset will start in 12h. 

   b. GPU. Set the gpu to run Contextual Regression Model. Use nvidia-smi command
      to check the gpu usage on the GPU Node. In Node gpu-1-0 of wanglab's server, 
      there are four gpus, which are 0,1,2,3. By default, it is set to 1. Use 
      salloc -p gpu -w gpu-1-0 to login gpu node. CR model uses gpu-1-0. 


3. "Software bin folder"

   a. RFolder. The absolute path of the bin folder for R.

   b. PythonFolder. The absolute path of the bin folder for Python.

   c. BedtoolFolder. The absolute path of the bin folder for bedtools.

   d. SamtoolFolder. The absolute path of the bin folder for samtools.



#---------------The Result folder and files-------------------------------

This script will create a sub-folder named Falcon_chromHMM_pvalue_1e-05 
folder(if p-value is set to 1e-05) in OutDir. In this sub-folder:

"SomaticRegion_Motif_hg19.bed" is the somatic mutation regions.

"Hist_of_L_SomaticRegion_hg19.pdf" is the histgram of length of regions.

"MotifDHS_hg19" folder is the data matrix for training Contextual Regression Model.

"Normalize0-1_cluterMotif_Zscore" folder is results of Contextual Regression Model. 
In total, 2431 motifs and open chromatin are the features. 
Considering the similarity among motifs, these motifs are clustered and one motif 
will be selected to represent the corresponding cluster. By setting cutoff h, 
different number of clusters will be obtained. When h=0.1, the cluster number is 
1210. When h=0.15, the cluster number is 884. When h=0.2, the cluster number is 
640. When h=0, we don't cluster these motifs. Under each h and each Epoch(5,10,15), 
we perform 10-fold cross validation for CR model. "CR_model_NodeSelf_MotifDHS_Epoch5_h0_1.h5"...
"CR_model_NodeSelf_MotifDHS_Epoch5_h0_10.h5" are the 10 CR models under h=0 and Epoch=5. 
"CR_model_NodeSelf_MotifDHS_ContributionModel_Epoch5_h0_1.h5"...
"CR_model_NodeSelf_MotifDHS_ContributionModel_Epoch5_h0_10.h5" are the 10 CR weight
model under h=0 and Epoch=5. "predictionCor_Summary_CV10Epoch5_h0.txt" is the pearson correlation
for 10-fold cv under h=0 and Epoch=5. Scatter plots for training set and testing set under 
different h are also included in this folder.


#------How to apply the trained CR model to other dataset----------------

Open ApplyCR_submit.job file and modify it. And then use the following 
command on the head node:

sbatch ApplyCR_submit.job

After this, you will get Result.txt: this file contains hyper and hypo p-value and FDR for each somatic region.

getHyperHypoRegion.R is script that is used to get the hyper and hyper regions.

Rscript getHyperHypoRegion.R $Out $FDR $Result

$Out is the output path.
$FDR is the cutoff for hyper and hypo regions
$Result is the Result.txt after ApplyCR_submit.job.




