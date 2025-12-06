When analysing all samples together, a few individuals in NC seem to be the most different, potentially as a spillover from the southern cluster that was believed to have stopped in Cape Hatteras. 


<img width="1274" height="459" alt="image" src="https://github.com/user-attachments/assets/09ab77fc-a19c-4142-b542-b178a48eca22" />


So, let's now rerun the analysis but without the NC population.

First create a new popmap without the NC samples

`cp BSB_all BSB_noNC`
 I manually went in and just deleted the NC samples but I'm sure there is a fancy way to do this with code =)

now create a new directory 

`mkdir /projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/stacks/final_noNC`

now edit the stacks and vcftools sbatch file to update the directories and make sure to change the `-p` to `5` since we removed a population. 

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short


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

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2

#--------------COMMAND----------------

ref_map.pl \
  -T 10 \
  -o $WOR_DIR/stacks/final_noNC \
  --popmap $WOR_DIR/popmap/BSB_all \
  --samples $WOR_DIR/samples \
  -X "populations: -r 0.80 --min-maf 0.01 -p 5 --write-single-snp --fstats --vcf --genepop --hwe --structure"

#--------- END Diagnostics/Logging Information---------------
echo = `date` job $JOB_NAME done
echo “using $NSLOTS CPUs”
```

output
```bash
Genotyped 492299 loci:
  effective per-sample coverage: mean=31.2x, stdev=9.5x, min=3.2x, max=48.5x
  mean number of sites per locus: 170.5
  a consistent phasing was found for 649400 of out 725405 (89.5%) diploid loci needing phasing
```

Now filter the vcf

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=vcftools              # Name your job something useful for easy tracking
#SBATCH --output=out/vcftools_%A.out
#SBATCH --error=out/vcftools_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short

#load program
module load vcftools # this is just globally available on explorer

# paths
INPUT_VCF="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/stacks/final_noNC/populations.snps.vcf"
OUTDIR="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/stacks/finalnoNC/"

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10
echo "done minimum depth"

# filter for sites present in >= 70% of individuals (sites need to be present in >=70% of individuals)
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --max-missing 0.7 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.7
echo "done maxmiss"

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.7.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10_maxmiss0.7.recode.vcf \
         --remove ${OUTDIR}/remove_individuals.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.7_filtInd
echo "done filtering pipeline"
```

```bash
After filtering, kept 96 out of 104 Individuals
Outputting VCF file...
After filtering, kept 11454 out of a possible 11454 Sites
Run Time = 2.00 seconds
```

individuals removed

```bash
(miniconda3) [r.gatins@explorer-01 final_noNC]$ cat remove_individuals.txt
INDV
MA_303_aligned_sorted
MA_304_aligned_sorted
MA_314_aligned_sorted
MA_316_aligned_sorted
MA_318_aligned_sorted
MA_324_aligned_sorted
ME_253_aligned_sorted
ME_254_aligned_sorted
```

move vcf file and popmap to my computer and run the analysis in R

```bash
/// GENLIGHT OBJECT /////////

 // 96 genotypes,  11,454 binary SNPs, size: 1.6 Mb
 88818 (8.08 %) missing data

 // Basic content
   @gen: list of 96 SNPbin

 // Optional content
   @ind.names:  96 individual labels
   @loc.names:  11454 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @other: a list containing: elements without names 
```

<img width="542" height="470" alt="image" src="https://github.com/user-attachments/assets/b619ca3e-4cc0-4f0c-845f-cccb1765c427" />
Structure with SNMF
<img width="1024" height="453" alt="image" src="https://github.com/user-attachments/assets/8f03e5fa-7762-4626-8228-1f253c250a76" />
PCA
<img width="678" height="563" alt="image" src="https://github.com/user-attachments/assets/3aef48e3-5fe6-42e5-9159-2ea4c4200de2" />
DAPC
<img width="826" height="571" alt="image" src="https://github.com/user-attachments/assets/6ab3e2eb-5e87-46d6-b6f5-67ed1e056e56" />

