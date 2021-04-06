#!/bin/bash

DEP_THRES=50    ##最低深度
MIN_DEP_THRES=5   ##支持alt的最小reads数
VALID_SIZE_THRES=100
FREQ_THRES=0.02   ###alt的频率最小值
STRANDED_RATIO_THRES=0.1     ## reads链偏性检测
# #download sickle and hammer/spades, install the executables in your $PATH
#sickleBinPath="~/bin/sickle-master"
#spadesBinPath="~/bin/SPAdes-3.7.1-Linux/"/public1/home/liuhj/project_2
fastqPath="./sample_reads"
cleanReadPath=/public1/home/liuhj2/process/Trimmomatic/${sample}
noThreadStep1=11
samtoolsThread=4
bowtieThread=24


customStep=$1   ##3,4,5
#ref_sample=$2

#refFasta=/public1/home/liuhj2/process/pilon/${ref_samp}/${ref_samp}.piloned.fasta     #"./sample_ref/Ebola_genome1.fa"
bowtie2indexPath=/public1/home/liuhj2/process/bwa/index/BGI_bwa2All_ONT_NECAT #"./sample_ref/Ebola_genome"
mpileupPath=/shared/liuhj/coffee/tonglv/Liangci/ntfreq/Same_sample/0147_0148             #"./mpileup_and_ntfreq"

noThreadStep3=6
tablePath=/shared/liuhj/tonglv/process/raw_ntfreq/table_0147_0148
mkdir $tablePath

excludeRegionF=/shared/liuhj/tonglv/process/ntfreq/exclude_region.txt


refFasta="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"     #"./sample_ref/Ebola_genome1.fa"



#########################################################################3
echo 
echo "********************************************************"
echo "*       Intrahost SNV calling for amplicon-seq         *"
echo "********************************************************"
echo "nim, 201709 v2"
echo

if [ $customStep == 0 ]; then
	#customStep="1,2,3,4,5"
	customStep="2,3"
else
	customStep=$1
	echo
	echo ~~~~~~~ custom step selection = $1 ~~~~~~~
	echo 
fi

if [[ $customStep =~ "1" ]];then
bash step1_fastqQC.sh $fastqPath $cleanReadPath  $noThreadStep1
fi



if [[ $customStep =~ "2" ]];then
bash step2_generate_mpileup.sh $cleanReadPath $samtoolsThread $bowtieThread $refFasta $bowtie2indexPath $mpileupPath


echo $cleanReadPath 
echo $samtoolsThread 
echo $bowtieThread 
echo $refFasta 
echo $bowtie2indexPath 
echo $mpileupPath 
fi

if [[ $customStep =~ "3" ]];then
bash step3_mpileup2ntfreq.sh $noThreadStep3 $mpileupPath
fi



#Usage: perl f_generate_summary_tables.pl -in <ntfreq dir> -ref <ref fasta> [-out output_path ]
if [[ $customStep =~ "4" ]];then
perl step4_generate_summary_tables.pl  -in $mpileupPath -ref $refFasta  -out $tablePath
fi

#Usage:	perl filtering_bigtables.pl -in <summary table dir> [options]
#Options:
#	-in <dir>          the path where *ntfreq files are.
# Thresholds:
#	-d <int>           the minimum depth of a genome position. Default = 100
#	-c <int>           the minimum depth of minor allele. Default = 5
#	-s <int>           the minimum number of valid (pass -d and -c) genome positions of a sample. Dafault = 3000
#	-f <fload|0-1>     the minimum minor allele frequency to identify a intrahost single nucleotide variations (iSNVs).If f > 1 - cut, output as SNP. Default = 0.02
#	-st-t <int>        the 4 types to measure stranded bias. 0: total. 1: max allele. 2: minor allele. 3: both max and minor alleles. 4: max or minor alleles. 
#					   Only type 2 available now. Default = 2
#	-st-c <fload|0-1>  the stranded bias ratio cutoff. Default = 0.1
#	-pf <file>         the file of genome region to exclude in iSNV identification, such as the primer regions. Format: stard-end\n... . Default = NULL
#	-pf-d <int>        the upstream and downstream adjcent region of -pf are also excluded, usually for -pf is primer regions.Default = 0
if [[ $customStep =~ "5" ]];then
#perl step5_filtering_bigtables.pl -in $tablePath -d 100 -c 5 -s 3000 -f 0.02 -st-t 2 -st-c 0.1 

echo perl step5_filtering_bigtables.pl -in $tablePath -d $DEP_THRES -c $MIN_DEP_THRES -s $VALID_SIZE_THRES -f $FREQ_THRES -st-t 2 -st-c $STRANDED_RATIO_THRES -output_sites_all_with_snp 1  # -pf  $excludeRegionF

perl step5_filtering_bigtables.pl -in $tablePath -d $DEP_THRES -c $MIN_DEP_THRES -s $VALID_SIZE_THRES -f $FREQ_THRES -st-t 2 -st-c $STRANDED_RATIO_THRES -output_sites_all_with_snp 1   # -pf  $excludeRegionF





mv  ./SNP_info.txt  $tablePath/
mv  ./iSNV_info.txt  $tablePath/
mv  ./samples.statistics.txt  $tablePath/
mv  ./iSNV.all.txt  $tablePath/
mv  ./iSNV.all.formatlab  $tablePath/
mv   ./iSNV_with_SNP.all.txt  $tablePath/




fi

echo ================ done ====================


