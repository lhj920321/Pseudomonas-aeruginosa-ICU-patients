#服务器
bwa_path=/public1/home/liuhj2/nanopore_HP_1/software/bwa-0.7.17
canu_path=/public1/home/liuhj2/nanopore_HP_1/software/canu-1.8/Linux-amd64/bin
java_1_8_path=/public1/home/liuhj2/software/jre1.8.0_191/bin
pilon_path=/public1/home/liuhj2/nanopore_HP_1/software
FastQC_path=/public1/home/liuhj2/nanopore_HP_1/software/FastQC
Trimmomatic_path=/public1/home/liuhj2/nanopore_HP_1/software/Trimmomatic-0.38
samtools_path=/public1/home/liuhj2/nanopore_HP_1/software/samtools-1.9
Prodigal_path=/public1/home/liuhj2/nanopore_HP_1/software/Prodigal-2.6.3
spades_path=/public1/home/liuhj2/nanopore_HP_1/software/SPAdes-3.7.1-Linux/bin
#虚拟机
Porechop_path=/home/liuhj/Documents/software/Porechop  
Unicycler_path=/home/liuhj/Documents/software/Unicycler

# python3_path=/public1/home/liuhj2/software/Python-3.7.2/python



sample_dir=$1
Type=$2
step=$3

samtools_threads=11


ref_genome_path=/public1/home/liuhj2/nanopore_HP_1/ref_genome
bwa_index_output=$ref_genome_path
ref_genome=$bwa_index_output/UM066.fasta
Trimmomatic_outpath=/public1/home/liuhj2/nanopore_HP_1/process/Trimmomatic/${sample_dir}


function bwa_run(){
####################
if [[ $Type == 'BGI' ]]; then
bwa_output=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/BGI/${sample_dir}	
Trimmomatic_outpath=/public1/home/liuhj2/nanopore_HP_1/process/Trimmomatic/${sample_dir}

mkdir $bwa_output
	echo $ref_genome
	echo $Trimmomatic_outpath
$bwa_path/bwa  mem  $bwa_index_output/UM066.fasta \
$Trimmomatic_outpath/${sample_dir}.R1.paired.fastq   $Trimmomatic_outpath/${sample_dir}.R2.paired.fastq  \
>$bwa_output/BGI_${sample_dir}.sam

$samtools_path/samtools view  --threads $samtools_threads  -bS  $bwa_output/BGI_${sample_dir}.sam  > $bwa_output/BGI_${sample_dir}.bam
$samtools_path/samtools sort  --threads $samtools_threads  $bwa_output/BGI_${sample_dir}.bam  -o  $bwa_output/BGI_${sample_dir}.sorted.bam
$samtools_path/samtools index $bwa_output/BGI_${sample_dir}.sorted.bam


rm $bwa_output/BGI_${sample_dir}.sam
rm $bwa_output/BGI_${sample_dir}.bam


fi


##################

if [[ $Type == 'ONT' ]]; then
bwa_output=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/ONT/ONT_${sample_dir}
ONT_data_P=/public1/home/liuhj2/nanopore_HP_1/data/BJXWZ-201810005Y-1/barcode/1.merge
ONT_data=$ONT_data_P/${sample_dir}.fastq 

mkdir $bwa_output

	
echo $ref_genome
echo $ONT_data
$bwa_path/bwa  mem  $ref_genome  $ONT_data >$bwa_output/ONT_${sample_dir}.sam


$samtools_path/samtools view  --threads $samtools_threads  -bS  $bwa_output/ONT_${sample_dir}.sam  > $bwa_output/ONT_${sample_dir}.bam
$samtools_path/samtools sort  --threads $samtools_threads  $bwa_output/ONT_${sample_dir}.bam  -o  $bwa_output/ONT_${sample_dir}.sorted.bam
$samtools_path/samtools index $bwa_output/ONT_${sample_dir}.sorted.bam


#rm $bwa_output/ONT_${sample_dir}.sam
#rm $bwa_output/ONT_${sample_dir}.bam

fi


}



function IGV_2k(){

start_point=50000
end_point=52000


if [[ $Type == 'BGI' ]]; then
bwa_output=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/BGI/${sample_dir}	
raw_bam_F=$bwa_output/BGI_${sample_dir}.sorted.bam  ##改

#out_P=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/IGV/BGI
out_P=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/mpileup


tiqu_sam_F=$out_P/${sample_dir}.${start_point}.${end_point}.sam
tiqu_sam_Header_F=$out_P/${sample_dir}.${start_point}.${end_point}.header.sam
tiqu_bam_F=$out_P/${sample_dir}.${start_point}.${end_point}.header.bam
tiqu_bam_sort=$out_P/BGI_${sample_dir}.${start_point}.${end_point}.sorted.bam

fi


if [[ $Type == 'ONT' ]]; then
bwa_output=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/ONT/ONT_${sample_dir}
raw_bam_F=$bwa_output/ONT_${sample_dir}.sorted.bam  ##改

#out_P=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/IGV/BGI
out_P=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/mpileup


tiqu_sam_F=$out_P/${sample_dir}.${start_point}.${end_point}.sam
tiqu_sam_Header_F=$out_P/${sample_dir}.${start_point}.${end_point}.header.sam
tiqu_bam_F=$out_P/${sample_dir}.${start_point}.${end_point}.header.bam
tiqu_bam_sort=$out_P/ONT_${sample_dir}.${start_point}.${end_point}.sorted.bam

fi


ref_ID_1=`samtools  view -H  $raw_bam_F | grep SQ  | awk '{print $2}'`
ref_ID=`echo ${ref_ID_1#*:}`
echo $ref_ID
cp $ref_genome  $out_P 
$samtools_path/samtools view  --threads  $samtools_threads  $raw_bam_F   ${ref_ID}:${start_point}-${end_point} >$tiqu_sam_F
$samtools_path/samtools view  --threads  $samtools_threads  -T $ref_genome  -h  $tiqu_sam_F  >$tiqu_sam_Header_F 
$samtools_path/samtools view  --threads  $samtools_threads -bS  $tiqu_sam_Header_F  >$tiqu_bam_F
$samtools_path/samtools sort  --threads  $samtools_threads  $tiqu_bam_F  -o  $tiqu_bam_sort
$samtools_path/samtools index  $tiqu_bam_sort


rm $tiqu_sam_F
rm $tiqu_sam_Header_F
rm $tiqu_bam_F

}



###只需输入一个样本号，就可将所有的
function Mpileup_ntfreq(){

ref_ID=UM066
echo $ref_ID

mpileupPath=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/mpileup


script_P=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/iSNV_calling
#bash $script_P/step2_generate_mpileup.sh  $Trimmomatic_outpath $samtools_threads  $samtools_threads  $ref_genome  $bwa_index_output  $mpileupPath


####mpileup2ntfreq
bash  $script_P/step3_mpileup2ntfreq.sh  $samtools_threads $mpileupPath

}

function snv_summary(){
DEP_THRES=20
MIN_DEP_THRES=2
VALID_SIZE_THRES=1000
FREQ_THRES=0.02
STRANDED_RATIO_THRES=0.1 


mpileupPath=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/mpileup
tablePath=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/mpileup/out

mkdir $tablePath
script_P=/public1/home/liuhj2/nanopore_HP_1/process/bwa/BGI_ONT_bwa2UM066/iSNV_calling

perl $script_P/step4_generate_summary_tables.pl  -in $mpileupPath -ref $ref_genome  -out $tablePath
perl $script_P/step5_filtering_bigtables.pl -in $tablePath -d $DEP_THRES -c $MIN_DEP_THRES -s $VALID_SIZE_THRES -f $FREQ_THRES -st-t 2 -st-c $STRANDED_RATIO_THRES -output_sites_all_with_snp 1

}



##1
if [[ $step == 'bwa' ]]; then
	bwa_run
fi

##2
if [[ $step == 'IGV' ]]; then
	IGV_2k
fi


###只需要输入一个样本即可
if [[ $step == 'Mpileup' ]]; then
	Mpileup_ntfreq
fi

if [[ $step == 'summary' ]]; then
	snv_summary
fi




