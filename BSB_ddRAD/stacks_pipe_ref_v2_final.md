# BSB RADseq STACKS pipeline with a updated reference genome `C_striata_v2.fasta`

Paths:

Working directory: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2`

Reference genome:   `/projects/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`  
Trimmed sequences: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags` 
Mapped sequences: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/samples`


Run Stacks

```bash
(base) [r.gatins@d3037 jobs]$ cat stacks_2.sh

#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
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

## Run populations

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
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

#--------------VARIABLES----------------

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/

#--------------COMMAND----------------

populations \
  -P $WOR_DIR/stacks/ref_map_1_stacks2.6 \
  --popmap $WOR_DIR/popmap/BSB_all \
  -r 0.80 \
  -p 1 \
  --min-maf 0.01 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O $WOR_DIR/populations/p1_maf_0.01 \
  -t 30
```

Run for `--min-maf` with `0.01` and `0.05` as well as `-p 1` and `-p 5`

|populations | No. Loci |
|------------|----------|
|p1_maf_0.01||
|p1_maf_0.05||
|p6_maf_0.01||
|p6_maf_0.05||
|maf_0.01||


```bash
#!/bin/bash

#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=vcftools
#SBATCH --output=vcftools_filter_%j.log
#SBATCH --error=vcftools_filter_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eppley.m@northeastern.edu

module load vcftools # this is just globally available on explorer, slay

# input VCF from stacks population run
INPUT_VCF="/home/r.gatins/BSB_ddRAD/stacks_ref_v2/populations.snps.vcf"

OUTDIR="/projects/gatins/2025_Mobulid/tarapacana/vcftools_filtered/n2_p1"
mkdir -p ${OUTDIR}

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10

# filter for sites present in >= 80% of individuals
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8

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
```


