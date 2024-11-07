# BSB RADseq STACKS pipeline with a updated reference genome `C_striata_v2.fasta`

Paths:

Working directory: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2`
or using the shortcut: `/home/r.gatins/BSB_ddRAD/stacks_ref_v2/`

Reference genome:   `/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`
Trimmed sequences: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags` 

## Demultiplex, trim, and quality check
Samples have already been demultiplexed previously (see Thais'[pipeline.md](https://github.com/thais-neu/BlackSeaBass_project/blob/master/BSB_ddRAD/pipeline.md))

For trimming and quality check see [stacks_pipeline.md](https://github.com/remygatins/BlackSeaBass_project/edit/master/BSB_ddRAD/stacks_pipeline.md)

## 1. Map files to the genome

### 1.1 Index Assembly
***BWA  v 0.7.17***

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
bwa index $DIR/C_striata_v2.fasta
```
bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa

```bash
(base) [r.gatins@login-00 final_genome]$ ls C_striata_v2.fasta*

C_striata_v2.fasta      C_striata_v2.fasta.ann  C_striata_v2.fasta.pac
C_striata_v2.fasta.amb  C_striata_v2.fasta.bwt  C_striata_v2.fasta.sa
```
### Run BWA array

Now, we will create a bwa job with tasks, or a job array. 

Let's work from the `/home/r.gatins/BSB_ddRAD/stacks_ref_v2/jobs` folder

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
#SBATCH --array=0-117%50		            #there are 118 samples and it will run a maximum of 10 jobs at a time

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

WOR_DIR=/home/r.gatins/BSB_ddRAD/stacks_ref_v2/
REF=/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta
SEQ_DIR=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags
SAMPLE_LIST=($(</work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/BSB_sample_list_uniq))
FILENAME=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

echo "My input file is ${FILENAME}"

bwa mem -t16 $REF $SEQ_DIR/${FILENAME}.1.fq.gz $SEQ_DIR/${FILENAME}.2.fq.gz > $WOR_DIR/samples/${FILENAME}_aligned.sam


echo = `date` job $JOB_NAME done

#--------- Diagnostic/Logging Information---------------

echo “using $NSLOTS CPUs”
echo `date`
```

<details>
<summary>Optimize your code</summary>
<br>

Use `seff job_ID` to get a summary of your run and adjust your parameters in your job script accordingly

```bash
(base) [r.gatins@login-01 samples]$ seff 45016272_23

Job ID: 45016357
Array Job ID: 45016272_23
Cluster: discovery
User/Group: r.gatins/users
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 5
CPU Utilized: 00:27:10
CPU Efficiency: 84.90% of 00:32:00 core-walltime
Job Wall-clock time: 00:06:24
Memory Utilized: 4.62 GB
Memory Efficiency: 18.47% of 25.00 GB
```
So in this case I will change my array to only call on 10G instead of 25G to run more efficiently 

</details>




### 1.2 convert  `.sam` file to `.bam` and sort using samtools 
***samtools v 1.9***

```
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=samtools             # Name your job something useful for easy tracking
#SBATCH --output=out/samtools_%A_%a.out
#SBATCH --error=out/samtools_%A_%a.err
#SBATCH --partition=lotterhos
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=10G                          # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=1-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
#SBATCH --array=0-117%50		#there are 118 samples and it will run a maximum of 50 jobs at a time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

echo “using $SLURM_CPUS_ON_NODE CPUs”
echo “Start Run”
echo “start time is `date`”

#--------------MODULES---------------

module load samtools

#--------------VARIABLES----------------
#Set variables/paths
WOR_DIR=/home/r.gatins/BSB_ddRAD/stacks_ref_v2/
REF=/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta
SEQ_DIR=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags
SAMPLE_LIST=($(</work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/BSB_sample_list_uniq))
FILENAME=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

#--------------COMMAND----------------
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "My input file is ${FILENAME}"

# convert .sam to .bam 
samtools view -Sb -@ 16 -O BAM -o $WOR_DIR/samples/${FILENAME}_aligned.bam $WOR_DIR/samples/${FILENAME}_aligned.sam

#sort output
samtools sort -o $WOR_DIR/samples/${FILENAME}_aligned_sorted.bam -O BAM -@ 16 $WOR_DIR/samples/${FILENAME}_aligned.bam

#--------- Diagnostic/Logging Information---------------
echo = `date` job $JOB_NAME done
echo `date`

```
Run time: 


PARAMETERS:
`-Sb`		input format bam\
`-@`		threads\
`-o` 		FILE  output file name\
`-O` 		output format (SAM, BAM, CRAM)\
`-b` 		output BAM\

### remove sam files

`.sam` files are much larger than `.bam` files. So, once we have converted sam to bam files we can delete the sam files to save space

`rm *sam`

### calculate percentage of reads mapped to reference 

```bash
module load samtools

samtools flagstat SN_191_aligned_sorted.bam > SN_191_aligned_sorted_stats.out
```
check output file
```
(base) [r.gatins@d3037 samples]$ cat SN_191_aligned_sorted_stats.out

3966287 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
108871 + 0 supplementary
0 + 0 duplicates
3948427 + 0 mapped (99.55% : N/A)
3857416 + 0 paired in sequencing
1928708 + 0 read1
1928708 + 0 read2
3616798 + 0 properly paired (93.76% : N/A)
3830340 + 0 with itself and mate mapped
9216 + 0 singletons (0.24% : N/A)
189760 + 0 with mate mapped to a different chr
92089 + 0 with mate mapped to a different chr (mapQ>=5)
```

run interactive mode 

    srun -p lotterhos -N 1 --pty /bin/bash
    
run stats for all samples using a for-loop. However, first create a file with all sample names `BSB_sample_list_uniq`

```bash
module load samtools

for file in `cat BSB_sample_list_uniq`;
do
    samtools flagstat ${file}_aligned_sorted.bam > ${file}_aligned_sorted_stats.out
done
```

print line 5 from all output files into a summary file to show percent mapped sequences 

    awk "FNR==5" *.out > flagstat_summary.txt

<details>
<summary>see flagstat_summary.txt</summary>
<br>

```bash
(base) [r.gatins@d3037 samples]$ cat flagstat_summary.txt

3019933 + 0 mapped (99.56% : N/A)
3700023 + 0 mapped (99.35% : N/A)
0 + 0 mapped (N/A : N/A)
1591197 + 0 mapped (99.48% : N/A)
159307 + 0 mapped (91.82% : N/A)
1132389 + 0 mapped (99.20% : N/A)
3220325 + 0 mapped (99.61% : N/A)
2653597 + 0 mapped (99.55% : N/A)
2614079 + 0 mapped (99.52% : N/A)
2754520 + 0 mapped (99.55% : N/A)
1913548 + 0 mapped (99.41% : N/A)
216286 + 0 mapped (97.83% : N/A)
2439143 + 0 mapped (99.46% : N/A)
601359 + 0 mapped (99.10% : N/A)
521182 + 0 mapped (98.17% : N/A)
2955987 + 0 mapped (99.55% : N/A)
3344079 + 0 mapped (99.56% : N/A)
2892574 + 0 mapped (99.62% : N/A)
154004 + 0 mapped (96.67% : N/A)
1828337 + 0 mapped (99.25% : N/A)
1176721 + 0 mapped (98.97% : N/A)
3072036 + 0 mapped (99.59% : N/A)
2626289 + 0 mapped (99.74% : N/A)
3008290 + 0 mapped (99.68% : N/A)
3393021 + 0 mapped (99.56% : N/A)
1840265 + 0 mapped (99.60% : N/A)
2262295 + 0 mapped (99.70% : N/A)
2436646 + 0 mapped (99.62% : N/A)
1898891 + 0 mapped (99.74% : N/A)
3318220 + 0 mapped (99.63% : N/A)
2959262 + 0 mapped (99.12% : N/A)
3454767 + 0 mapped (99.64% : N/A)
3274770 + 0 mapped (99.60% : N/A)
2922702 + 0 mapped (99.66% : N/A)
1717195 + 0 mapped (99.71% : N/A)
1970811 + 0 mapped (99.61% : N/A)
1829953 + 0 mapped (98.42% : N/A)
2088208 + 0 mapped (99.63% : N/A)
2601217 + 0 mapped (99.66% : N/A)
2962264 + 0 mapped (98.00% : N/A)
1894622 + 0 mapped (97.15% : N/A)
2311953 + 0 mapped (99.48% : N/A)
2208838 + 0 mapped (99.66% : N/A)
3184820 + 0 mapped (99.51% : N/A)
2727634 + 0 mapped (99.50% : N/A)
2960233 + 0 mapped (99.58% : N/A)
3470495 + 0 mapped (99.42% : N/A)
4106189 + 0 mapped (99.49% : N/A)
2386576 + 0 mapped (92.46% : N/A)
1635133 + 0 mapped (98.01% : N/A)
2650901 + 0 mapped (98.09% : N/A)
2347653 + 0 mapped (67.61% : N/A)
2387934 + 0 mapped (70.24% : N/A)
3879468 + 0 mapped (99.52% : N/A)
3266308 + 0 mapped (99.65% : N/A)
4181329 + 0 mapped (99.54% : N/A)
3415092 + 0 mapped (99.60% : N/A)
3668197 + 0 mapped (99.43% : N/A)
3477542 + 0 mapped (99.51% : N/A)
2981084 + 0 mapped (91.42% : N/A)
3033797 + 0 mapped (99.57% : N/A)
3488455 + 0 mapped (99.58% : N/A)
3390097 + 0 mapped (99.29% : N/A)
2659236 + 0 mapped (99.50% : N/A)
1510839 + 0 mapped (59.73% : N/A)
1858077 + 0 mapped (49.65% : N/A)
3374232 + 0 mapped (99.11% : N/A)
2214836 + 0 mapped (59.56% : N/A)
1999444 + 0 mapped (99.49% : N/A)
1875251 + 0 mapped (74.61% : N/A)
2176212 + 0 mapped (99.43% : N/A)
1154297 + 0 mapped (89.81% : N/A)
1844101 + 0 mapped (99.68% : N/A)
1789786 + 0 mapped (99.45% : N/A)
2164564 + 0 mapped (99.60% : N/A)
1778311 + 0 mapped (99.55% : N/A)
2818595 + 0 mapped (99.58% : N/A)
3265478 + 0 mapped (99.57% : N/A)
3034932 + 0 mapped (99.58% : N/A)
1716473 + 0 mapped (99.59% : N/A)
1511968 + 0 mapped (99.60% : N/A)
2638481 + 0 mapped (99.56% : N/A)
2339741 + 0 mapped (99.67% : N/A)
2589376 + 0 mapped (99.56% : N/A)
3203067 + 0 mapped (99.57% : N/A)
3233846 + 0 mapped (99.60% : N/A)
3245552 + 0 mapped (99.59% : N/A)
2972584 + 0 mapped (99.59% : N/A)
2583460 + 0 mapped (99.68% : N/A)
2990187 + 0 mapped (99.52% : N/A)
3458137 + 0 mapped (99.52% : N/A)
2139039 + 0 mapped (99.52% : N/A)
3417092 + 0 mapped (99.65% : N/A)
3723592 + 0 mapped (99.62% : N/A)
3086659 + 0 mapped (99.65% : N/A)
2436045 + 0 mapped (99.56% : N/A)
2035721 + 0 mapped (99.68% : N/A)
3087931 + 0 mapped (99.55% : N/A)
2767226 + 0 mapped (99.63% : N/A)
2464771 + 0 mapped (99.58% : N/A)
3824511 + 0 mapped (99.63% : N/A)
3946485 + 0 mapped (99.59% : N/A)
4062153 + 0 mapped (99.55% : N/A)
3460483 + 0 mapped (99.60% : N/A)
2455671 + 0 mapped (99.62% : N/A)
3361716 + 0 mapped (99.51% : N/A)
2973346 + 0 mapped (99.69% : N/A)
3515476 + 0 mapped (99.59% : N/A)
3123325 + 0 mapped (99.61% : N/A)
3517921 + 0 mapped (99.61% : N/A)
2263580 + 0 mapped (99.60% : N/A)
3280117 + 0 mapped (99.58% : N/A)
2502499 + 0 mapped (99.62% : N/A)
3176334 + 0 mapped (99.53% : N/A)
2790539 + 0 mapped (99.61% : N/A)
2531255 + 0 mapped (99.46% : N/A)
3722023 + 0 mapped (99.55% : N/A)
3948427 + 0 mapped (99.55% : N/A)

```
</details>

The third sample seems to have failed, from MA_300. 

`samtools flagstat MA_300_aligned_sorted.bam > MA_300_aligned_sorted_stats.out`

Looking back at the original .fq.gz file, this file contains no sequences so we will drop it from the analysis. 

## Run Stacks

First make a popmap file with all your sample `.bam` file names

This can be done in excel and copied into a file in our `/popmap` directory

<details>
<summary>BSB_all popmap</summary>
<br>
    
```bash
(base) [r.gatins@login-00 popmap]$ cat BSB_all

MA_298_aligned_sorted	NM
MA_299_aligned_sorted	NM
MA_302_aligned_sorted	NM
MA_303_aligned_sorted	NM
MA_304_aligned_sorted	NM
MA_306_aligned_sorted	NM
MA_307_aligned_sorted	NM
MA_310_aligned_sorted	NM
MA_311_aligned_sorted	NM
MA_313_aligned_sorted	NM
MA_314_aligned_sorted	NM
MA_315_aligned_sorted	NM
MA_316_aligned_sorted	NM
MA_318_aligned_sorted	NM
MA_320_aligned_sorted	NM
MA_321_aligned_sorted	NM
MA_323_aligned_sorted	NM
MA_324_aligned_sorted	NM
MA_325_aligned_sorted	NM
MA_327_aligned_sorted	NM
MD_136_aligned_sorted	MD
MD_137_aligned_sorted	MD
MD_138_aligned_sorted	MD
MD_139_aligned_sorted	MD
MD_140_aligned_sorted	MD
MD_141_aligned_sorted	MD
MD_142_aligned_sorted	MD
MD_143_aligned_sorted	MD
MD_145_aligned_sorted	MD
MD_149_aligned_sorted	MD
MD_150_aligned_sorted	MD
MD_151_aligned_sorted	MD
MD_152_aligned_sorted	MD
MD_154_aligned_sorted	MD
MD_158_aligned_sorted	MD
MD_159_aligned_sorted	MD
MD_160_aligned_sorted	MD
MD_161_aligned_sorted	MD
MD_162_aligned_sorted	MD
MD_163_aligned_sorted	MD
ME_164_aligned_sorted	ME
ME_165_aligned_sorted	ME
ME_166_aligned_sorted	ME
ME_167_aligned_sorted	ME
ME_176_aligned_sorted	ME
ME_248_aligned_sorted	ME
ME_249_aligned_sorted	ME
ME_250_aligned_sorted	ME
ME_251_aligned_sorted	ME
ME_252_aligned_sorted	ME
ME_253_aligned_sorted	ME
ME_254_aligned_sorted	ME
ME_255_aligned_sorted	ME
ME_256_aligned_sorted	ME
ME_257_aligned_sorted	ME
ME_258_aligned_sorted	ME
ME_261_aligned_sorted	ME
ME_262_aligned_sorted	ME
NC_233_aligned_sorted	NC
NC_234_aligned_sorted	NC
NC_235_aligned_sorted	NC
NC_237_aligned_sorted	NC
NC_238_aligned_sorted	NC
NC_239_aligned_sorted	NC
NC_240_aligned_sorted	NC
NC_241_aligned_sorted	NC
NC_242_aligned_sorted	NC
NC_243_aligned_sorted	NC
NC_244_aligned_sorted	NC
NC_245_aligned_sorted	NC
NC_246_aligned_sorted	NC
NJ_106_aligned_sorted	NJ
NJ_108_aligned_sorted	NJ
NJ_109_aligned_sorted	NJ
NJ_112_aligned_sorted	NJ
NJ_113_aligned_sorted	NJ
NJ_114_aligned_sorted	NJ
NJ_118_aligned_sorted	NJ
NJ_119_aligned_sorted	NJ
NJ_121_aligned_sorted	NJ
NJ_122_aligned_sorted	NJ
NJ_124_aligned_sorted	NJ
NJ_128_aligned_sorted	NJ
NJ_129_aligned_sorted	NJ
NJ_130_aligned_sorted	NJ
NJ_131_aligned_sorted	NJ
NJ_132_aligned_sorted	NJ
NJ_133_aligned_sorted	NJ
RI_328_aligned_sorted	SN
RI_329_aligned_sorted	SN
RI_330_aligned_sorted	SN
RI_331_aligned_sorted	SN
RI_332_aligned_sorted	SN
RI_333_aligned_sorted	SN
RI_334_aligned_sorted	SN
RI_335_aligned_sorted	SN
RI_336_aligned_sorted	SN
RI_337_aligned_sorted	SN
RI_338_aligned_sorted	SN
RI_339_aligned_sorted	SN
RI_340_aligned_sorted	SN
RI_341_aligned_sorted	SN
RI_342_aligned_sorted	SN
RI_343_aligned_sorted	SN
RI_344_aligned_sorted	SN
RI_345_aligned_sorted	SN
RI_346_aligned_sorted	SN
RI_347_aligned_sorted	SN
RI_348_aligned_sorted	SN
RI_349_aligned_sorted	SN
SN_009_aligned_sorted	SN
SN_179_aligned_sorted	SN
SN_182_aligned_sorted	SN
SN_185_aligned_sorted	SN
SN_189_aligned_sorted	SN
SN_190_aligned_sorted	SN
SN_191_aligned_sorted	SN
```
</details>

### run stacks

ref_map.pl --samples [path] --popmap [path] -o [path] [--rm-pcr-duplicates] [-X prog:"opts" ...]

ref_map.pl -T 10 -o $WOR_DIR/stacks --popmap $WOR_DIR/popmap/BSB_all --samples $WOR_DIR/samples --rm-pcr-duplicates -X "populations: -r 0.80 --min-maf 0.01 --fstats --vcf --genepop --hwe --structure"

<details>
<summary>See Stacks job script</summary>
<br>
    
```bash
(base) [r.gatins@d3037 jobs]$ cat stacks.sh

#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
##SBATCH --array=0-117%50		#there are 118 samples and it will run a maximum of 10 jobs at a time

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

echo “using $SLURM_CPUS_ON_NODE CPUs”
echo “Start Run”
echo “start time is `date`”

#--------------MODULES---------------

module load lotterhos
module load stacks

#--------------START diagnostics----------------

echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

#--------------VARIABLES----------------

WOR_DIR=/home/r.gatins/BSB_ddRAD/stacks_ref_v2/

#--------------COMMAND----------------

ref_map.pl -T 10 -o $WOR_DIR/stacks --popmap $WOR_DIR/popmap/BSB_all --samples $WOR_DIR/samples --rm-pcr-duplicates -X "populations: -r 0.80 --min-maf 0.01 --fstats --vcf --genepop --hwe --structure"

#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.80 --vcf --genepop --structure --fstats --hwe -t 10
#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.80 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10
#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.80 --min-maf 0.05 --vcf --genepop --structure --fstats --hwe -t 10
#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all -r 0.60 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10

#write-single-snp
#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_all --write-single-snp -r 0.80 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10
#filter out NC samples
#populations -P $DIR/stacks/ -M $DIR/popmap/BSB_noNC --write-single-snp -r 0.80 --min-maf 0.01 --vcf --genepop --structure --fstats --hwe -t 10

#--------- END Diagnostics/Logging Information---------------
echo = `date` job $JOB_NAME done
echo “using $NSLOTS CPUs”
```
</details>

 








