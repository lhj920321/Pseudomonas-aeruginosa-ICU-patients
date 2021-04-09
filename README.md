# title

## summary

## Dependencies

- python 3.7
- perl 

- RAxML
- R version 3.5.1
- bwa
- samtools
- Spades
- 

## Repo Contents

- meta-data: sample data and preprocessed data used for analysis.

  PA ref genome:

  PA ref gff file :

  dir: genomes of the other bacteria 

- results: 

- scripts:

## SNPs calling

1. map to reference genome :

2. SNPs calling 

   iSNV_calling.sh 3

   iSNV_calling.sh 4



## repeat region 

- cut the ref genome to 150bp simulation sequences:

  python  cut 

- map the  simulation reads to the other bacteria genomes

- 



## Phylogenetic analysis

- prepare SNPs sequences:

​     ` python XXX`

- Phylogenetic analysis was 

`raxmlHPC-PTHREADS -s $align_trim_file -n raxmltree_result -m GTRGAMMAI -f a -x 12345 -N 1000 -p 123456 -T 24 -k`



## Genome assembly

`genome assembly : Spade XXX`



## variations annotation

`snpeff XXX`

1. To convert iSNV table into vcf format(VCFv4.1) that can be recognized  by SnpEff software: `python ./scripts/Variantions_annotion/iSNVTable_2_vcf_allsample.py  -i  ./data/   -o $allsamplesVcfPath   -r   MN908947.3`
2. Annoting the mutation by using SnpEff software ,for example: `java -jar snpEff.jar ann  MN908947   $$allsamplesVcfPath/${sample}.vcf    >$outSnpeff/{sample}.snpeffAnno.vcf`

## Citation

If you use data, results or conclusion from this work, please cite:



## Acknowledgement



