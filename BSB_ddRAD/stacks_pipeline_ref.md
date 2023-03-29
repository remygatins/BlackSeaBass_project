# BSB RADseq STACKS pipeline with a Reference Genome




Lotterhos lab partitions

There are 2 nodes each with 36 cores, so we can run parallel programs with up to 36 cores, and serial on up to 72 cores.

There is 5GB of memory on each core, with 36 cores/node and 72 cores total across both. So it should be possible to submit an array with 140 jobs at a time, each with 2.5 GB

Paths:

Working directory: /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref
Reference genome:   /work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_01.fasta
Trimmed sequences: /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags 

## Demultiplex, trim, and quality check
Samples have already been demultiplexed previously (see Thais'[pipeline.md](https://github.com/thais-neu/BlackSeaBass_project/blob/master/BSB_ddRAD/pipeline.md))

For trimming and quality check see [stacks_pipeline.md](https://github.com/remygatins/BlackSeaBass_project/edit/master/BSB_ddRAD/stacks_pipeline.md)

## Map files to the genome


### BWA  v 0.7.17

map trimmed Illumina reads to our draft genome

We first create an index assembly file of the draft genome C_striata_01.fasta


```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=bwa              # Name your job something useful for easy tracking
#SBATCH --output=out/bwa.out
#SBATCH --error=out/bwa.err
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
##SBATCH --array=0-117%10		#there are 118 samples and it will run a maximum of 10 jobs at a time

#--------------MODULES---------------

module load bwa

#--------------COMMAND----------------

DIR=/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome

#index genome
bwa index $DIR/C_striata_01.fasta
```

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa



### Create bwa job array ###
Now, we will create a bwa job with tasks, or a job array. 

Let's work from the `../stacks_ref/jobs` folder

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=bwa              # Name your job something useful for easy tracking
#SBATCH --output=out/bwa_%A_%a.out
#SBATCH --error=out/bwa_%A_%a.err
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=25G                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --array=0-117%50		#there are 118 samples and it will run a maximum of 10 jobs at a time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

echo “using $SLURM_CPUS_ON_NODE CPUs”
echo “Start Run”
echo “start time is `date`”

#--------------MODULES---------------

module load bwa

#--------------COMMAND----------------

echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

REF=/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_01.fasta
SEQ_DIR=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags
SAMPLE_LIST=($(</work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/BSB_sample_list_uniq))
FILENAME=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

echo "My input file is ${FILENAME}"

bwa mem -t16 $REF $SEQ_DIR/${FILENAME}.1.fq.gz $SEQ_DIR/${FILENAME}.2.fq.gz > $SEQ_DIR/${FILENAME}_aligned.sam


echo = `date` job $JOB_NAME done



#--------- Diagnostic/Logging Information---------------

echo “using $NSLOTS CPUs”
echo `date`

```

convert  `.sam` file to `.bam` and sort using samtools 

```
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=samtools             # Name your job something useful for easy tracking
#SBATCH --output=out/samtools_%A_%a.out
#SBATCH --error=out/samtools_%A_%a.err
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=25G                          # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=1-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
#SBATCH --array=0-117%10		#there are 118 samples and it will run a maximum of 10 jobs at a time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

echo “using $SLURM_CPUS_ON_NODE CPUs”
echo “Start Run”
echo “start time is `date`”

#--------------MODULES---------------

module load samtools

#--------------COMMAND----------------
#Set variables/paths
REF=/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_01.fasta
SEQ_DIR=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags
SAMPLE_LIST=($(</work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref/samples/BSB_sample_list_uniq))
FILENAME=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}


echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "My input file is ${FILENAME}"

# convert .sam to .bam 
samtools view -Sb -@ 16 -O BAM -o $SEQ_DIR/${FILENAME}_aligned.bam $SEQ_DIR/${FILENAME}_aligned.sam

#sort output
samtools sort -o $SEQ_DIR/${FILENAME}_aligned_sorted.bam -O BAM -@ 16 $SEQ_DIR/${FILENAME}_aligned.bam


#--------- Diagnostic/Logging Information---------------
echo = `date` job $JOB_NAME done
echo `date`

```


`-Sb`		input format bam\
`-@`		threads\
`-o` 		FILE  output file name\
`-O` 		output format (SAM, BAM, CRAM)\
`-b` 		output BAM\

`.sam` files are much larger than `.bam` files. 
Once we have converted sam to bam files we can delete the sam files to save space

`rm *sam`

#calculate percentage of reads mapped to reference 

```bash
module load samtools

samtools flagstat SN_191_aligned_sorted.bam > SN_191_aligned_sorted_stats.out
```

```
3966250 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
108834 + 0 supplementary
0 + 0 duplicates
3948393 + 0 mapped (99.55% : N/A)
3857416 + 0 paired in sequencing
1928708 + 0 read1
1928708 + 0 read2
3616750 + 0 properly paired (93.76% : N/A)
3830346 + 0 with itself and mate mapped
9213 + 0 singletons (0.24% : N/A)
189774 + 0 with mate mapped to a different chr
92016 + 0 with mate mapped to a different chr (mapQ>=5)
```

run interactive mode

    srun -p lotterhos -N 1 --pty /bin/bash
    
run stats for all samples using a for loop

```bash
module load samtools

for file in `cat BSB_sample_list_uniq`;
do
    samtools flagstat ${file}_aligned_sorted.bam > ${file}_aligned_sorted_stats.out
done
```

print line 5 from all output files into a summary file

    awk "FNR==5" *.out > flagstat_summary.txt

## Run Stacks

I'm going to move all bam files from my `stacks` directory to my `stacks_ref` directory

        mv *.bam /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref/samples/

### Reference map

Make a popmap file with all your sample `.bam` file names


```bash
module load lotterhos
module load stacks

DIR=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref

ref_map.pl -T 10 -o $DIR/stacks --popmap $DIR/popmap/BSB_all --samples $DIR/samples --rm-pcr-duplicates -X "populations: --fstats --vcf --genepop
```

<img width="472" alt="image" src="https://user-images.githubusercontent.com/26288352/178889797-cd5e9522-9934-4f5a-9146-986732716cdc.png">


### Populations r 0.80

```bash
populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.80 --vcf --genepop --structure --fstats --hwe -t 10
```

### Populations r 0.80 min-maf 0.01
```bash
populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.80 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 30
```

### Populations r 0.80 min-maf 0.05

```bash
populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.80 --min-maf 0.05 --vcf --genepop --structure --fstats --hwe -t 30
```

### Populations r 0.60 min-maf 0.01

```bash
populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.60 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10
```


## Summary table

|popmap |populations        | Loci   |% missing data|
|:------|:-----------------:|:------:|-------------:|
|BSB_all| none          |213,852  | 47.3%  |
|BSB_all|-r 0.80        |101,818    |30.14% |
|BSB_all|-r 0.80 --min-maf 0.01 | 90,571  | 31.99%|
|BSB_all|-r 0.80 --min-maf 0.05 | 61,135  | 33.53%|
|BSB_all|-r 0.60 --min-maf 0.01 |   | |
|BSB_all|-r 0.80 --min-maf 0.01 --write-single-snp|24,896|21.22 |


I'm going to use the vcf file where I filtered loci to show up in 80% of individuals per population `-r 0.80` and had a minimum allele frequecy of 0.01 `--min-maf 0.01`

We first need to convert our VCF file (.vcf) into a PLINK format (.bed/.bim/.fam)

```bash
module load miniconda3
conda activate plink

plink --vcf populations.snps.vcf --out BSB_plink --double-id --allow-extra-chr
```

`--double-id` My ind ID's had a "_" 
`--allow-extra-chr` chromosome name was too long and needs this to accept

This should have created 3 files:
```bash
BSB_plink.bed
BSB_plink.bim
BSB_plink.fam
```

Now open OpenOnDemand to work with R interactively

Open `bigsnpr.R`


------------------- [11/03/2022] -----------

run interactive mode

    srun -p lotterhos -N 1 --pty /bin/bash

1. Create chromosome map

```bash
module load miniconda3
conda activate plink
conda activate bcftools

bcftools view -H filename.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > filename.chrom-map.txt
```

2. Make ped file using this chromosome map

```bash
module load vcftools

vcftools --vcf populations.snps.vcf --out BSB_r0.6_maf0.01_v2 --plink --chrom-map filename.chrom-map.txt
```
Now open OpenOnDemand to work with R interactively

https://ood.discovery.neu.edu/pun/sys/dashboard

An R course has already been created by Katie with most packages already installed. So go to `Courses` > `EEMB3465_RStudio`

Select the time and memory you need


# March 28, 2023

Adegenet_BSB_ref.R

```r
rm(list = ls())
#install.packages("adegenet", dep=TRUE)

library("adegenet")
library("ape") 
library("RColorBrewer")
library("poppr") # create trees and treat missing data

#show Adegenet version
packageDescription("adegenet", fields = "Version") 

#set working directory
setwd("/Users/remygatins/GoogleDrive_gmail/Work/Projects/2021_Black\ Sea\ Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.01")

#Load .vcf data output from populations (STACKS)
#install.packages("vcfR")
library(vcfR)
vcf <- read.vcfR(file = "populations.snps.vcf")
#vcf <- read.vcfR(file = "filtered_ind_2.recode.vcf")
#vcf <- read.vcfR(file = "filtered.recode.vcf")
data <- vcfR2genlight(vcf)
data
```
<img width="436" alt="image" src="https://user-images.githubusercontent.com/26288352/228348369-8e2892b2-67ac-485f-bebf-77b5dc3de066.png">


```r
#__________________
#Label populations 
#------------------

#popmap1
pop(data)<- c(rep("NM",20), rep("MD",20), rep("ME",18), rep("NC",13), rep("NJ",17), rep("SN",29)) 


Dgen <- dist.gene(data)
Dgeo <- dist(data)
ploidy(data) <- 2

##############
# Check data #
##############
#plot missing data
glPlot(data, posi="topleft")
glPlot(data, col=bluepal(6))
```

<img width="543" alt="image" src="https://user-images.githubusercontent.com/26288352/228347855-aaae9e91-e96d-4b93-a5fb-138e0a98b444.png">




