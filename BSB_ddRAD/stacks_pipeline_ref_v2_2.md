# BSB RADseq STACKS pipeline with an updated reference genome `C_striata_v2.fasta`

Paths:

Working directory: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2`
Reference genome:   `/projects/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`  
Trimmed sequences: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags` 


## Run Ref_map
`cat stacks_2.sh`

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short
##SBATCH --array=0-117%50		#there are 118 samples and it will run a maximum of 10 jobs at a time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

echo “using $SLURM_CPUS_ON_NODE CPUs”
echo “Start Run”
echo “start time is `date`”

#--------------MODULES---------------

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

#--------------START diagnostics----------------

echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

#--------------VARIABLES----------------

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/

#--------------COMMAND----------------

ref_map.pl -T 10 -o $WOR_DIR/stacks/ref_map_1_stacks2.6 --popmap $WOR_DIR/popmap/BSB_all --samples $WOR_DIR/samples --rm-pcr-duplicates -X "populations: -r 0.80 --min-maf 0.01 --fstats --vcf --genepop --hwe --structure"

#--------- END Diagnostics/Logging Information---------------
echo = `date` job $JOB_NAME done
echo “using $NSLOTS CPUs”
```


## Now run Populations
Within the `populations` directory, make the following directories

`mkdir p1_maf_0.01  p1_maf_0.05  p6_maf_0.01  p6_maf_0.05`

Now run the different populations

```bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2


populations \
  -P $WOR_DIR/stacks/ref_map_1_stacks2.6 \
  --popmap --popmap $WOR_DIR/popmap/BSB_all \
  -r 0.80 \
  -p 1 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O $WOR_DIR/populations/p1_maf_0.01 \
  -t 30
```

To obtain number of loci from each populations run
`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|populations|No. Loci|
|-----------|--------|
|maf 0.01|27664| 
|p1_maf_0.01|28707|
|p1_maf_0.05|28509|
|p6_maf_0.01|12693|
|p6_maf_0.05|12692|

I'm going with p1 maf 0.01 since it gives us the most loci and we will be doing more filtering with vcftools which may remove those extra loci anyway

## Filter VCF with VCFTOOLS


```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=vcftools              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short

#load program
module load vcftools # this is just globally available on explorer

# paths
INPUT_VCF="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01/populations.snps.vcf"
OUTDIR="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01/"

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10
echo "done minimum depth"

# filter for sites present in >= 80% of individuals (refiltered any SNPs found in ≤80% of individuals)
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8
echo "done max maxmiss"

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.8.recode.vcf \
         --remove ${OUTDIR}/remove_individuals.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8_filtInd
echo "done filtering piepline"
```


output
```bash

# filter by minimum depth per genotype (minDP = 10)

After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 24993 out of a possible 24993 Sites
Run Time = 4.00 seconds

# filter for sites present in >= 80% of individuals (refiltered any SNPs found in ≤80% of individuals)
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 0 out of a possible 24993 Sites
No data left for analysis!
Run Time = 0.00 seconds
```
Filtering by 0.8 max missingness seems to stringent, so I tried running it with 0.5 instead and keep having the same issue.
I believe it has to do with the format of the file that may be missing something from the stacks output. I will filter this directly within the populations command using -R 0.8.

```bash
export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2

populations \
  -P $WOR_DIR/stacks/ref_map_1_stacks2.6 \
  --popmap $WOR_DIR/popmap/BSB_all \
  -r 0.80 \
  -p 1 \
  --min-maf 0.01 \
  --R 0.80 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O $WOR_DIR/populations/p1_maf_0.01_R0.8 \
  -t 30
```

`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|populations|No. Loci|
|-----------|--------|
|p1_maf_0.01_R0.8|13692|

ok now I will filter again with vcftools but i'll skip the `max_missingness` as we did this already 


```bash
#load program
module load vcftools # this is just globally available on explorer

# paths
INPUT_VCF="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01_R0.8/populations.snps.vcf"
OUTDIR="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01_R0.8/"

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10
echo "done minimum depth"

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --remove ${OUTDIR}/remove_individuals.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_filtInd
echo "done filtering piepline"
```

