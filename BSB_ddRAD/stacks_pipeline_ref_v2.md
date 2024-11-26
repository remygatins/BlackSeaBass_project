# BSB RADseq STACKS pipeline with a updated reference genome `C_striata_v2.fasta`

Paths:

Working directory: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2`
or using the shortcut: `/home/r.gatins/BSB_ddRAD/stacks_ref_v2/`

Reference genome:   `/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`  
Trimmed sequences: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags` 

## Demultiplex, trim, and quality check
Samples have already been demultiplexed previously (see Thais'[pipeline.md](https://github.com/thais-neu/BlackSeaBass_project/blob/master/BSB_ddRAD/pipeline.md))

For trimming and quality check see [stacks_pipeline.md](https://github.com/remygatins/BlackSeaBass_project/edit/master/BSB_ddRAD/stacks_pipeline.md)

In summary:
Raw sequences were quality-checked using FastQC v0.12.1 and low-quality bases and adapter sequences were removed with TRIMGALORE. We used STACKS v 2.41 (Catchen et al., 2013; Rochette et al., 2019) to further quality filter, trim sequences to 130 bp, and remove PCR clones using process_radtags and clone_filter.

```bash

module load lotterhos
module load stacks

#--------------COMMAND----------------

for file in `cat ../samples/BSB_sample_list`;
do
   process_radtags -1 ../samples/no_adapter/${file}.F_val_1.fq.gz \
   -2 ../samples/no_adapter/${file}.R_val_2.fq.gz \
   -o ../samples/no_adapter/process_radtags \
   --renz-1 mspI --renz-2 bamHI -c -q -t 130
done

# Rename sequences

for f in *.F_val_1.1.fq.gz;
    do mv "$f" "${f%.F_val_1.1.fq.gz}.1.fq.gz";
    done

for f in *.R_val_2.2.fq.gz; do mv "$f" "${f%.R_val_2.2.fq.gz}.2.fq.gz"; done
```

Let's clone filter using stacks

```bash
#--------------MODULES---------------
module load lotterhos
module load stacks

#--------------COMMAND----------------
WOR_DIR=/home/r.gatins/BSB_ddRAD/stacks/

for file in `cat $WOR_DIR/samples/BSB_sample_list`;
do
    clone_filter -1 $WOR_DIR/samples/no_adapter/process_radtags/${file}.1.fq.gz \
    -2 $WOR_DIR/samples/no_adapter/process_radtags/${file}.2.fq.gz \
    -o $WOR_DIR/samples/no_adapter/process_radtags/clone_filter \
    -i gzfastq
done
```


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

#--------- END Diagnostics/Logging Information---------------
echo = `date` job $JOB_NAME done
echo “using $NSLOTS CPUs”
```
</details>

Job statistics

(base) [r.gatins@login-00 out]$ seff 45020667
Job ID: 45020667
Cluster: discovery
User/Group: r.gatins/users
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 10
CPU Utilized: 03:44:04
CPU Efficiency: 77.35% of 04:49:40 core-walltime
Job Wall-clock time: 00:28:58
Memory Utilized: 1010.54 MB
Memory Efficiency: 4.93% of 20.00 GB

## Filter VCF file
***vcftools v0.1.17***

from our stacks directory lets now make a new filtering directory

```bash
mkdir filtering
cd filtering
```
Run interactive mode
```bash
srun -p lotterhos -N 1 --pty /bin/bash
module load vcftools
```

To make this file more manageable, let’s start by applying three step filter. We are going to only keep variants that have been successfully genotyped in 50% of individuals, a minimum quality score of 30, and a minor allele count of 3.

`vcftools --gzvcf ../populations.snps.vcf --max-missing 0.5 --mac 1 --minQ 30 --recode --recode-INFO-all --out snps.g5mac3`

In this code, we call vcftools, feed it a vcf file after the `--vcf` flag, `--max-missing 0.5` tells it to filter genotypes called below 50% (across all individuals) the `--mac 3` flag tells it to filter SNPs that have a minor allele count less than 1.

The --recode flag tells the program to write a new vcf file with the filters, `--recode-INFO-all` keeps all the INFO flags from the old vcf file in the new one. Lastly, `--out` designates the name of the output. The output will scroll through a lot of lines, but should end like:

```bash
After filtering, kept 40 out of 40 Individuals
After filtering, kept 78434 out of a possible 147540 Sites
Outputting VCF file... Done
Run Time = 40.00 seconds 
Those two simple filters got rid of 50% of the data and will make the next filtering steps run much faster.
```
Missing 
vcftools --vcf ../populations.snps.vcf --max-missing 0.5 --recode --recode-INFO-all --out snps.g5

```bash
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 64570 out of a possible 90578 Sites
Run Time = 10.00 seconds
```
check missing data per ind
vcftools --vcf snps.g5dp10.recode.vcf --missing-indv

```bash 
(base) [r.gatins@d3037 filtering]$ cat out.imiss.g5
INDV	N_DATA	N_GENOTYPES_FILTERED	N_MISS	F_MISS
MA_298_aligned_sorted	64570	0	17250	0.267152
MA_299_aligned_sorted	64570	0	16799	0.260167
MA_302_aligned_sorted	64570	0	19123	0.296159
MA_303_aligned_sorted	64570	0	34830	0.539415
MA_304_aligned_sorted	64570	0	19639	0.304151
MA_306_aligned_sorted	64570	0	17206	0.26647
MA_307_aligned_sorted	64570	0	17078	0.264488
MA_310_aligned_sorted	64570	0	17975	0.27838
MA_311_aligned_sorted	64570	0	17659	0.273486
MA_313_aligned_sorted	64570	0	17364	0.268917
MA_314_aligned_sorted	64570	0	24860	0.385009
MA_315_aligned_sorted	64570	0	17435	0.270017
MA_316_aligned_sorted	64570	0	27116	0.419947
MA_318_aligned_sorted	64570	0	36102	0.559114
MA_320_aligned_sorted	64570	0	17244	0.267059
MA_321_aligned_sorted	64570	0	16925	0.262119
MA_323_aligned_sorted	64570	0	16995	0.263203
MA_324_aligned_sorted	64570	0	29137	0.451247
MA_325_aligned_sorted	64570	0	17158	0.265727
MA_327_aligned_sorted	64570	0	17630	0.273037
MD_136_aligned_sorted	64570	0	15308	0.237076
MD_137_aligned_sorted	64570	0	16321	0.252764
MD_138_aligned_sorted	64570	0	17560	0.271953
MD_139_aligned_sorted	64570	0	15312	0.237138
MD_140_aligned_sorted	64570	0	15804	0.244758
MD_141_aligned_sorted	64570	0	15308	0.237076
MD_142_aligned_sorted	64570	0	16724	0.259006
MD_143_aligned_sorted	64570	0	19190	0.297197
MD_145_aligned_sorted	64570	0	15100	0.233855
MD_149_aligned_sorted	64570	0	15429	0.23895
MD_150_aligned_sorted	64570	0	15228	0.235837
MD_151_aligned_sorted	64570	0	14888	0.230571
MD_152_aligned_sorted	64570	0	15153	0.234676
MD_154_aligned_sorted	64570	0	18784	0.290909
MD_158_aligned_sorted	64570	0	18033	0.279278
MD_159_aligned_sorted	64570	0	18085	0.280084
MD_160_aligned_sorted	64570	0	17006	0.263373
MD_161_aligned_sorted	64570	0	15491	0.23991
MD_162_aligned_sorted	64570	0	15493	0.239941
MD_163_aligned_sorted	64570	0	18277	0.283057
ME_164_aligned_sorted	64570	0	8094	0.125352
ME_165_aligned_sorted	64570	0	9423	0.145935
ME_166_aligned_sorted	64570	0	7921	0.122673
ME_167_aligned_sorted	64570	0	7934	0.122874
ME_176_aligned_sorted	64570	0	8117	0.125709
ME_248_aligned_sorted	64570	0	7834	0.121326
ME_249_aligned_sorted	64570	0	7771	0.12035
ME_250_aligned_sorted	64570	0	9029	0.139833
ME_251_aligned_sorted	64570	0	8451	0.130881
ME_252_aligned_sorted	64570	0	7958	0.123246
ME_253_aligned_sorted	64570	0	59148	0.916029
ME_254_aligned_sorted	64570	0	55543	0.860198
ME_255_aligned_sorted	64570	0	7731	0.119731
ME_256_aligned_sorted	64570	0	7905	0.122425
ME_257_aligned_sorted	64570	0	7634	0.118228
ME_258_aligned_sorted	64570	0	7870	0.121883
ME_261_aligned_sorted	64570	0	7785	0.120567
ME_262_aligned_sorted	64570	0	7825	0.121186
NC_233_aligned_sorted	64570	0	4823	0.0746941
NC_234_aligned_sorted	64570	0	5115	0.0792164
NC_235_aligned_sorted	64570	0	5119	0.0792783
NC_237_aligned_sorted	64570	0	5294	0.0819885
NC_238_aligned_sorted	64570	0	5302	0.0821124
NC_239_aligned_sorted	64570	0	5953	0.0921945
NC_240_aligned_sorted	64570	0	4818	0.0746167
NC_241_aligned_sorted	64570	0	4751	0.0735791
NC_242_aligned_sorted	64570	0	5343	0.0827474
NC_243_aligned_sorted	64570	0	5810	0.0899799
NC_244_aligned_sorted	64570	0	5153	0.0798049
NC_245_aligned_sorted	64570	0	5363	0.0830571
NC_246_aligned_sorted	64570	0	5208	0.0806567
NJ_106_aligned_sorted	64570	0	5034	0.0779619
NJ_108_aligned_sorted	64570	0	4535	0.0702339
NJ_109_aligned_sorted	64570	0	4721	0.0731144
NJ_112_aligned_sorted	64570	0	5398	0.0835992
NJ_113_aligned_sorted	64570	0	4475	0.0693046
NJ_114_aligned_sorted	64570	0	4430	0.0686077
NJ_118_aligned_sorted	64570	0	4194	0.0649528
NJ_119_aligned_sorted	64570	0	4783	0.0740746
NJ_121_aligned_sorted	64570	0	4506	0.0697847
NJ_122_aligned_sorted	64570	0	4368	0.0676475
NJ_124_aligned_sorted	64570	0	4561	0.0706365
NJ_128_aligned_sorted	64570	0	4424	0.0685148
NJ_129_aligned_sorted	64570	0	4550	0.0704662
NJ_130_aligned_sorted	64570	0	4176	0.064674
NJ_131_aligned_sorted	64570	0	4435	0.0686851
NJ_132_aligned_sorted	64570	0	4043	0.0626142
NJ_133_aligned_sorted	64570	0	4580	0.0709308
RI_328_aligned_sorted	64570	0	2538	0.0393062
RI_329_aligned_sorted	64570	0	2354	0.0364566
RI_330_aligned_sorted	64570	0	2522	0.0390584
RI_331_aligned_sorted	64570	0	2610	0.0404212
RI_332_aligned_sorted	64570	0	2338	0.0362088
RI_333_aligned_sorted	64570	0	2914	0.0451293
RI_334_aligned_sorted	64570	0	3491	0.0540654
RI_335_aligned_sorted	64570	0	3063	0.0474369
RI_336_aligned_sorted	64570	0	3091	0.0478705
RI_337_aligned_sorted	64570	0	2754	0.0426514
RI_338_aligned_sorted	64570	0	2800	0.0433638
RI_339_aligned_sorted	64570	0	2512	0.0389035
RI_340_aligned_sorted	64570	0	2211	0.0342419
RI_341_aligned_sorted	64570	0	2395	0.0370915
RI_342_aligned_sorted	64570	0	2346	0.0363327
RI_343_aligned_sorted	64570	0	2557	0.0396004
RI_344_aligned_sorted	64570	0	2691	0.0416757
RI_345_aligned_sorted	64570	0	3000	0.0464612
RI_346_aligned_sorted	64570	0	2446	0.0378814
RI_347_aligned_sorted	64570	0	2716	0.0420629
RI_348_aligned_sorted	64570	0	2861	0.0443085
RI_349_aligned_sorted	64570	0	10834	0.167787
SN_009_aligned_sorted	64570	0	2777	0.0430076
SN_179_aligned_sorted	64570	0	3365	0.052114
SN_182_aligned_sorted	64570	0	2780	0.043054
SN_185_aligned_sorted	64570	0	4778	0.0739972
SN_189_aligned_sorted	64570	0	2497	0.0386712
SN_190_aligned_sorted	64570	0	2755	0.0426669
SN_191_aligned_sorted	64570	0	2813	0.0435651
```


with min depth 10

`vcftools --vcf snps.g5.recode.vcf --minDP 10 --recode --recode-INFO-all --out snps.g5dp10`

```bash
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 64570 out of a possible 64570 Sites
Run Time = 9.00 seconds
```
check missing data per individual

`vcftools --vcf snps.g5dp10.recode.vcf --missing-indv`

```bash
(base) [r.gatins@d3037 filtering]$ cat out.imiss.g5dp10
INDV	N_DATA	N_GENOTYPES_FILTERED	N_MISS	F_MISS
MA_298_aligned_sorted	64570	0	64570	1
MA_299_aligned_sorted	64570	0	64570	1
MA_302_aligned_sorted	64570	0	64570	1
MA_303_aligned_sorted	64570	0	64570	1
MA_304_aligned_sorted	64570	0	64570	1
MA_306_aligned_sorted	64570	0	64570	1
MA_307_aligned_sorted	64570	0	64564	0.999907
MA_310_aligned_sorted	64570	0	64569	0.999985
MA_311_aligned_sorted	64570	0	64570	1
MA_313_aligned_sorted	64570	0	64570	1
MA_314_aligned_sorted	64570	0	64570	1
MA_315_aligned_sorted	64570	0	64570	1
MA_316_aligned_sorted	64570	0	64570	1
MA_318_aligned_sorted	64570	0	64570	1
MA_320_aligned_sorted	64570	0	64570	1
MA_321_aligned_sorted	64570	0	64570	1
MA_323_aligned_sorted	64570	0	64570	1
MA_324_aligned_sorted	64570	0	64570	1
MA_325_aligned_sorted	64570	0	64570	1
MA_327_aligned_sorted	64570	0	64570	1
MD_136_aligned_sorted	64570	0	64565	0.999923
MD_137_aligned_sorted	64570	0	64570	1
MD_138_aligned_sorted	64570	0	64570	1
MD_139_aligned_sorted	64570	0	64570	1
MD_140_aligned_sorted	64570	0	64570	1
MD_141_aligned_sorted	64570	0	64570	1
MD_142_aligned_sorted	64570	0	64570	1
MD_143_aligned_sorted	64570	0	64570	1
MD_145_aligned_sorted	64570	0	64570	1
MD_149_aligned_sorted	64570	0	64570	1
MD_150_aligned_sorted	64570	0	64567	0.999954
MD_151_aligned_sorted	64570	0	64570	1
MD_152_aligned_sorted	64570	0	64570	1
MD_154_aligned_sorted	64570	0	64570	1
MD_158_aligned_sorted	64570	0	64570	1
MD_159_aligned_sorted	64570	0	64570	1
MD_160_aligned_sorted	64570	0	64570	1
MD_161_aligned_sorted	64570	0	64566	0.999938
MD_162_aligned_sorted	64570	0	64570	1
MD_163_aligned_sorted	64570	0	64566	0.999938
ME_164_aligned_sorted	64570	0	64570	1
ME_165_aligned_sorted	64570	0	64570	1
ME_166_aligned_sorted	64570	0	64568	0.999969
ME_167_aligned_sorted	64570	0	64570	1
ME_176_aligned_sorted	64570	0	64570	1
ME_248_aligned_sorted	64570	0	64570	1
ME_249_aligned_sorted	64570	0	64570	1
ME_250_aligned_sorted	64570	0	64570	1
ME_251_aligned_sorted	64570	0	64570	1
ME_252_aligned_sorted	64570	0	64570	1
ME_253_aligned_sorted	64570	0	64570	1
ME_254_aligned_sorted	64570	0	64570	1
ME_255_aligned_sorted	64570	0	64570	1
ME_256_aligned_sorted	64570	0	64570	1
ME_257_aligned_sorted	64570	0	64570	1
ME_258_aligned_sorted	64570	0	64570	1
ME_261_aligned_sorted	64570	0	64570	1
ME_262_aligned_sorted	64570	0	64570	1
NC_233_aligned_sorted	64570	0	64570	1
NC_234_aligned_sorted	64570	0	64570	1
NC_235_aligned_sorted	64570	0	64570	1
NC_237_aligned_sorted	64570	0	64570	1
NC_238_aligned_sorted	64570	0	64570	1
NC_239_aligned_sorted	64570	0	64570	1
NC_240_aligned_sorted	64570	0	64570	1
NC_241_aligned_sorted	64570	0	64570	1
NC_242_aligned_sorted	64570	0	64570	1
NC_243_aligned_sorted	64570	0	64570	1
NC_244_aligned_sorted	64570	0	64570	1
NC_245_aligned_sorted	64570	0	64570	1
NC_246_aligned_sorted	64570	0	64570	1
NJ_106_aligned_sorted	64570	0	64570	1
NJ_108_aligned_sorted	64570	0	64570	1
NJ_109_aligned_sorted	64570	0	64570	1
NJ_112_aligned_sorted	64570	0	64570	1
NJ_113_aligned_sorted	64570	0	64570	1
NJ_114_aligned_sorted	64570	0	64570	1
NJ_118_aligned_sorted	64570	0	64570	1
NJ_119_aligned_sorted	64570	0	64570	1
NJ_121_aligned_sorted	64570	0	64570	1
NJ_122_aligned_sorted	64570	0	64570	1
NJ_124_aligned_sorted	64570	0	64570	1
NJ_128_aligned_sorted	64570	0	64570	1
NJ_129_aligned_sorted	64570	0	64570	1
NJ_130_aligned_sorted	64570	0	64570	1
NJ_131_aligned_sorted	64570	0	64570	1
NJ_132_aligned_sorted	64570	0	64570	1
NJ_133_aligned_sorted	64570	0	64570	1
RI_328_aligned_sorted	64570	0	64570	1
RI_329_aligned_sorted	64570	0	64570	1
RI_330_aligned_sorted	64570	0	64570	1
RI_331_aligned_sorted	64570	0	64570	1
RI_332_aligned_sorted	64570	0	64563	0.999892
RI_333_aligned_sorted	64570	0	64570	1
RI_334_aligned_sorted	64570	0	64570	1
RI_335_aligned_sorted	64570	0	64570	1
RI_336_aligned_sorted	64570	0	64570	1
RI_337_aligned_sorted	64570	0	64570	1
RI_338_aligned_sorted	64570	0	64570	1
RI_339_aligned_sorted	64570	0	64570	1
RI_340_aligned_sorted	64570	0	64570	1
RI_341_aligned_sorted	64570	0	64570	1
RI_342_aligned_sorted	64570	0	64570	1
RI_343_aligned_sorted	64570	0	64570	1
RI_344_aligned_sorted	64570	0	64570	1
RI_345_aligned_sorted	64570	0	64570	1
RI_346_aligned_sorted	64570	0	64570	1
RI_347_aligned_sorted	64570	0	64570	1
RI_348_aligned_sorted	64570	0	64570	1
RI_349_aligned_sorted	64570	0	64570	1
SN_009_aligned_sorted	64570	0	64570	1
SN_179_aligned_sorted	64570	0	64570	1
SN_182_aligned_sorted	64570	0	64570	1
SN_185_aligned_sorted	64570	0	64568	0.999969
SN_189_aligned_sorted	64570	0	64570	1
SN_190_aligned_sorted	64570	0	64570	1
SN_191_aligned_sorted	64570	0	64570	1
```

Why do we have more missing data after removing? In the manual:
```bash
--minDP <float>
--maxDP <float>
           Includes only genotypes greater than or equal to the "--minDP" value and less than or equal to the "--maxDP" value. This option requires that the "DP" FORMAT tag is specified for all sites.
```

so maybe because some sites don't seem to have a DP value it is not doing a good job?

calculate Depth per individual after removing 0.5 missing data

`vcftools --vcf snps.g5.recode.vcf --depth`

```bash
(base) [r.gatins@d3037 filtering]$ cat out.idepth
INDV	N_SITES	MEAN_DEPTH
MA_298_aligned_sorted	47320	1.06989
MA_299_aligned_sorted	47771	1.10452
MA_302_aligned_sorted	45447	1.04801
MA_303_aligned_sorted	29740	1.01258
MA_304_aligned_sorted	44931	1.03659
MA_306_aligned_sorted	47364	1.06338
MA_307_aligned_sorted	47492	1.07397
MA_310_aligned_sorted	46595	1.07303
MA_311_aligned_sorted	46911	1.0767
MA_313_aligned_sorted	47206	1.06933
MA_314_aligned_sorted	39710	1.02339
MA_315_aligned_sorted	47135	1.08809
MA_316_aligned_sorted	37454	1.02819
MA_318_aligned_sorted	28468	1.03576
MA_320_aligned_sorted	47326	1.08374
MA_321_aligned_sorted	47645	1.09394
MA_323_aligned_sorted	47575	1.0684
...
```

Now calculate mean depth per site

`vcftools --vcf snps.g5.recode.vcf --site-mean-depth`

```bash
CHROM   POS     MEAN_DEPTH      VAR_DEPTH
Scaffold_1      22453   1.02778 0.0273865
Scaffold_1      22350   1.02817 0.0277666
Scaffold_1      22233   1.02778 0.0273865
Scaffold_1      244873  2.03896 1.59057
Scaffold_1      244847  2.09836 1.75683
Scaffold_1      244776  1.96154 1.49201
Scaffold_1      300317  1.01887 0.0186882
Scaffold_1      308317  1.07792 0.0727956
Scaffold_1      308491  1.06897 0.0647676
Scaffold_1      308470  1.06034 0.0571964
Scaffold_1      308470  1       0
Scaffold_1      310050  1       0
Scaffold_1      310274  1.03509 0.0341562
Scaffold_1      310279  1.02609 0.0256293
Scaffold_1      310287  1.03448 0.0335832
Scaffold_1      310339  1.03419 0.0333039
Scaffold_1      351152  1.08696 0.0800915
Scaffold_1      351227  1.08036 0.0745656
Scaffold_1      351321  1.08696 0.0800915
Scaffold_1      351330  1.00943 0.00943396
...
```

thin 

`vcftools --vcf snps.g5.recode.vcf --thin 5000 --recode --recode-INFO-all --out snps.g5t5000`

```bash
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 15146 out of a possible 64570 Sites
Run Time = 3.00 seconds
```

The next step is to get rid of individuals that did not sequence well. We can do this by assessing individual levels of missing data.

vcftools --vcf snps.g5t5000.recode.vcf --missing-indv



#----------------------------------------------------------------
running parameters from https://researchdata.gla.ac.uk/1059/2/biogeography_whitefish.Rmd

I created a new directory within samples

```bash
mkdir q20
cd 20
```

comvert sam to bam files but only retained if mapping quality was >20 with samtools v.1.7 (Li et al., 2009)
```bash
# convert .sam to .bam
samtools view -Sbq 20 -@ 16 -O BAM -o $WOR_DIR/samples/q20/${FILENAME}_aligned.bam $WOR_DIR/samples/q20/${FILENAME}_aligned.sam

#sort output
samtools sort -o $WOR_DIR/samples/q20/${FILENAME}_aligned_sorted.bam -O BAM -@ 16 $WOR_DIR/samples/q20/${FILENAME}_aligned.bam
```
## Build Catalog loci

```{bash, eval = FALSE}
WOR_DIR=/home/r.gatins/BSB_ddRAD/stacks_ref_v2/

ref_map.pl -T 10 --popmap $WOR_DIR/popmap/BSB_all --samples $WOR_DIR/samples/q20 -o $WOR_DIR/stacks/ref_map_q20 --rm-pcr-duplicates
```
### Population genomics and phylogenetics analyses 

From this point onward, we are filtering and generating datasets for different analyses.

#### Generate a vcf file using [populations](http://catchenlab.life.illinois.edu/stacks/comp/populations.php)

Here we are telling `populations` to retain loci only if present in at least 80% of individuals per population in at least 9 populations. We are also removing loci with heterozygosity higher than 0.6 and minor allele frequency below 0.05. Only one snp per locus is retained to reduce the effect of linkage.
re
```{bash, eval = FALSE}
#run with same parameters as before to see if anything changes
populations -P $WOR_DIR/stacks/ref_map_q20 -M $WOR_DIR/popmap/BSB_all -t 10 -r 0.8 --min_maf 0.01 --vcf -O $WOR_DIR/stacks/ref_map_q20/r0.8_maf0.01

populations -P $WOR_DIR/stacks/ref_map_q20 -M $WOR_DIR/popmap/BSB_all -t 4 -p 9 -r 0.8 --max_obs_het 0.6 --min_maf 0.01 --write_single_snp --vcf

```

`‐p 9` (minimum number of populations a locus must be present in to be retained)
`‐r 0.8` (minimum proportion of samples in a population required to have a locus)
`‐‐max_obs_het 0.6` (maximum observed heterozygosity for a nucleotide at a locus)
`‐‐min_maf 0.05` (minor allele frequency across populations required for a SNP to be included in the dataset)
`‐‐write_single_snp` (retain one SNP per locus)


```bash
srun -p lotterhos -N 1 --pty /bin/bash
module load vcftools

vcftools --vcf populations.snps.vcf --max-missing 0.5 --recode --recode-INFO-all --out snps.g5
```

```bash
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 63596 out of a possible 88920 Sites
Run Time = 10.00 seconds
```
```bash
vcftools --vcf populations.snps.vcf --max-missing 0.5 --minDP 5 --recode --recode-INFO-all --out snps.g5dp5
```

After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 6 out of a possible 88920 Sites
Run Time = 2.00 seconds

populations -P $WOR_DIR/stacks/ref_map_q20 -M $WOR_DIR/popmap/BSB_all -t 10 -r 0.8 --min_maf 0.01 --write-single-snp -r 0.7 -R 0.5 --vcf -O $WOR_DIR/stacks/ref_map_q20/r0.8_maf0.01_snp1_r0.7_R0.5


# Denovo Map comparison
with a denovo map

```bash
WOR_DIR=/home/r.gatins/BSB_ddRAD/stacks_ref_v2/

denovo_map.pl --samples $WOR_DIR/samples/ --popmap $WOR_DIR/popmap/BSB_all -o $WOR_DIR/stacks/denovo --paired --rm-pcr-duplicates -m 4 -M 3 -n 2 -T 12

```

depth from denovo map output

```bash
(base) [r.gatins@d3037 BSB_pop_r0.8_minmaf0.01]$ cat out.idepth
INDV	N_SITES	MEAN_DEPTH
MA_298	62616	54.0363
MA_299	61749	73.1809
MA_302	47564	35.1175
MA_304	34745	28.1168
MA_306	62226	58.2528
MA_307	62676	48.8396
MA_310	60930	42.8115
MA_311	61439	46.9788
MA_313	49141	36.1576
MA_315	62238	47.4712
MA_320	63071	52.3185
MA_321	63003	60.8062
MA_323	62284	53.7729
MA_325	59572	36.8039
MA_327	57648	21.134
MD_136	48918	61.0627
MD_137	47729	53.8987
MD_138	47812	60.3542
MD_139	47772	68.271
MD_140	40344	38.8521
MD_141	45828	50.2015
MD_142	48791	50.5635
MD_143	42058	37.8581
MD_145	49044	72.1611
MD_149	49306	64.0227
MD_150	49133	73.8174
MD_151	48888	71.2082
MD_152	49154	61.558
MD_154	43581	36.5593
MD_158	45075	38.4997
MD_159	46067	36.7049
MD_160	47589	45.9335
MD_161	47391	50.5915
MD_162	49125	63.0058
MD_163	42479	34.8322
ME_164	74214	37.4779
ME_165	71654	37.9465
ME_166	92126	53.8507
ME_167	91340	43.3377
ME_176	90745	45.209
ME_248	91924	54.7131
ME_249	91864	64.1217
ME_250	80392	33.849
ME_251	59154	30.1621
ME_252	89897	40.287
ME_255	92105	60.388
ME_256	91184	50.6543
ME_257	92091	69.0608
ME_258	92311	55.5227
ME_261	92108	58.2058
ME_262	92403	55.6806
NC_233	77945	48.5908
NC_234	74516	50.6483
NC_235	78006	60.2312
NC_237	77658	58.3382
NC_238	66772	44.5515
NC_239	67245	25.4385
NC_240	78624	30.1026
NC_241	78305	56.8471
NC_242	76286	35.9531
NC_243	77329	35.8159
NC_244	75383	30.9245
NC_245	78565	39.2575
NC_246	59876	21.37
NJ_106	98732	31.0582
NJ_108	101969	30.7783
NJ_109	102094	35.7574
NJ_112	99384	29.5756
NJ_113	102019	46.2554
NJ_114	102543	53.5239
NJ_118	102421	50.0896
NJ_119	81574	28.1472
NJ_121	68971	25.7374
NJ_122	102894	44.5737
NJ_124	98384	36.2231
NJ_128	102405	42.5261
NJ_129	102602	52.4054
NJ_130	102667	53.1598
NJ_131	102432	52.5133
NJ_132	101978	45.9503
NJ_133	99791	40.8591
RI_328	105334	49.6823
RI_329	103405	51.1087
RI_330	89261	33.6694
RI_331	105636	55.8545
RI_332	105496	60.3625
RI_333	104343	48.6725
RI_334	82481	37.9824
RI_335	99900	32.8473
RI_336	101483	42.9803
RI_337	104644	44.9268
RI_338	105042	40.6128
RI_339	104597	58.2468
RI_340	105513	62.3818
RI_341	105654	66.0332
RI_342	105079	52.7224
RI_343	99340	37.6658
RI_344	105468	56.8856
RI_345	104633	48.5511
RI_346	105965	56.426
RI_347	105479	51.378
RI_348	105625	58.5115
RI_349	73507	37.7108
SN_009	99427	49.4208
SN_179	94002	37.3859
SN_182	101164	50.6852
SN_185	93425	40.0842
SN_189	100626	39.9464
SN_190	100929	58.3985
SN_191	100906	58.8745
```










