#!/bin/bash

cleanPath=$1 #"./clean_fastq"
samtoolsThread=$2 #8
bowtieThread=$3 #16
reffasta=$4 #"./ref/Ebola_genome1.fa"
bowtie2indexpath=$5 #"./ref/Ebola_genome"
outputpath=$6 #"./mpileup_and_ntfreq"

#bowtie2binpath=""
#samtoolsbinpath=""

echo 
echo ============== step 2: generate mpileup "(samtools 1.3 required)" ==============
echo

if [ "${cleanpath:0-1}"x = "/"x ]; then
	cleanpath=${cleanpath%?}
fi

echo "clean fastq path: "$cleanPath
echo "reference fasta file: "$reffasta

echo "output path: "$outputpath

# ----------------------------------
echo -------------- step 2.1: bowtie2-align and sort bam -------------


echo
echo --------------- step 2.2: generate mpileup file -----------------
for file in $outputpath/*sorted.bam
do
	echo $file
	while [ $(ps -Af|grep "samtools mpileup -f"|wc -l) -gt $samtoolsThread ]
	do
		sleep 2
	done
	if [ ! -e ${file/sorted.bam/mpileup} ];then
		samtools mpileup -f $reffasta -d 1000000 $file  > ${file/sorted.bam/mpileup} &
	else
		echo ${file/sorted.bam/mpileup} already exists, skip samtools mpileup.
	fi
done
	
while [ $(ps -Af|grep "samtools mpileup -f"|wc -l) -gt 1 ]; do
	sleep 2
done

echo


