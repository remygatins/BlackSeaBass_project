GenomeScope & PSMC
================
Remy Gatins
last edited: June 27, 2022

  - [GenomeScope](#genomescope)
  - [PSMC](#psmc)
      - [1. Prepare your genome assembly
        data](#1-prepare-your-genome-assembly-data)
      - [2. Call diploid genome and covert to psmc
        file](#2-call-diploid-genome-and-covert-to-psmc-file)
      - [3. Run PSMC analysis](#3-run-psmc-analysis)
      - [4. Plot PSMC results](#4-plot-psmc-results)
  - [PSMC with Bootstrap](#psmc-with-bootstrap)
      - [1. Prepare your genome assembly
        data](#1-prepare-your-genome-assembly-data-1)
      - [2. Call diploid genome](#2-call-diploid-genome)
      - [3. Run PSMC analysis](#3-run-psmc-analysis-1)
      - [4. Generate split file for
        bootstrap](#4-generate-split-file-for-bootstrap)
      - [5. Run PSMC in bootstrap mode using a job
        array](#5-run-psmc-in-bootstrap-mode-using-a-job-array)
      - [6. Plot PSMC results](#6-plot-psmc-results)

## Check Raw PacBio sequences

``` bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=fastqc              # Name your job something useful for easy tracking
#SBATCH --output=out/fastqc.out
#SBATCH --error=out/fastqc.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
##SBATCH --array=0-117%10		#there are 118 samples and it will run a maximum of 10 jobs at a time

#--------------MODULES---------------

module load oracle_java/jdk1.8.0_181
module load fastqc/0.11.8

#--------------COMMAND----------------

fastqc /work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/PacBio_Denovo/raw_sequences/DTG-DNA-1126.r64296e173242G01.subreads_ccs.fastq.gz -o fastqc

```

<img width="1055" alt="image" src="https://user-images.githubusercontent.com/26288352/176452603-67306086-5f57-48ea-a192-bfdc7f00c1ad.png">


## GenomeScope

A program used to estimate genome heterozigosity, repeat content, and
size from sequencing reads using a kmer-based statistical approach

1.  You will first need to download and install Jellyfish:
    <http://www.genome.umd.edu/jellyfish.html#Release>

2.  Once you have jellyfish on your server you now need to count K-mers
    from your raw Illumina fastq data

<!-- end list -->

Run an intereactive node and load modules
```bash
srun -p lotterhos -N 1 --pty /bin/bash
module load jellyfish
```

``` bash
jellyfish count -C -m 21 -s 60000000000 -t 12 /work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/PacBio_Denovo/raw_sequences/DTG-DNA-1126.r64296e173242G01.subreads_ccs.fastq.gz -o BSB_PacBio.jf
```

Jellyfish can’t read gzipped files so if you have fq.gz files use the
following command:

``` bash
jellyfish count -C -m 21 -s 60000000000 -t 12 <(zcat /work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/PacBio_Denovo/raw_sequences/DTG-DNA-1126.r64296e173242G01.subreads_ccs.fastq.gz) -o BSB_PacBio.jf
```
  Run time 07:34:07

*Parameter key:*  
`-C` indicates to count canonical kmers (**don’t change this**) 
`-m` kmer length (default -m 21)  
`-s` allocate memory  
`-t` threads  
`-o` output name (default is mer\_counts.jf)

3.  Now export the kmer count histogram

<!-- end list -->

``` bash
jellyfish histo -t 12 BSB_PacBio.jf > BSB_PacBio.histo
```
  Run time 00:00:59

4.  Upload reads.histo to GenomeScope: <http://qb.cshl.edu/genomescope/>

<img width="1221" alt="image" src="https://user-images.githubusercontent.com/26288352/177917337-76494f83-6b99-4eca-9e0e-d1b41ff4f03d.png">

[GenomeScope analysis](http://genomescope.org/analysis.php?code=FJQq5xxcIAPcCIvhw0dz)

After looking at your GenomeScope results and have a good idea of your
average coverage, continue on to the PSMC analysis.


## PSMC

The PSMC analysis can infer historic population size from a diploid
sequence using the Pairwise Sequentially Markovian Coalescent (PSMC)
model.

Specific details and information on this package can be found:
<https://github.com/lh3/psmc>

#### 1. Prepare your genome assembly data

First we need to index our genome assembly (should be a fasta file)


``` bash
#load programs
module load bwa

#set the working directory
DIR=/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome

#run bwa
bwa index $DIR/final_genome/C_striata_01.fasta 
```
  Run time 00:16:41

Now map adapter trimmed Illumina reads to your genome:

``` bash
bwa mem -t32 $DIR/final_genome/C_striata_01.fasta $DIR/PacBio_Denovo/raw_sequences/DTG-DNA-1126.r64296e173242G01.subreads_ccs.fastq.gz > $DIR/PSMC/C_striata_01.sam
```
  Run time 06:26:02

Use samtools to convert your sam file to a bam file:

``` bash
module load samtools

samtools view -Sb -@ 30 -O BAM -o $DIR/PSMC/C_striata_01.bam $DIR/PSMC/C_striata_01.sam 
```
  Run time 00:06:15

Sort your bam file:

``` bash
samtools sort -o $DIR/PSMC/C_striata_01_sorted.bam -O BAM -@ 20 $DIR/PSMC/C_striata_01.bam 
```
  Run time 00:07:51

Finally, index your sorted bam file:

``` bash
samtools index -b -@ 20 $DIR/PSMC/C_striata_01_sorted.bam  
```
  Run time 00:00:40

#### 2\. Call diploid genome and covert to psmc file

Once you have mapped and sorted your genome the next step is to call a
diploid genome and convert it to the required format to input into the
PSMC analysis.

So, first we call a diploid genome using our reference genome and mapped
reads using bcftools. 

NOTE: I installed bcftools and PSMC within the PSMC conda environment. To activate:
 `conda activate /work/lotterhos/programs/PSMC`

``` bash
bcftools mpileup -C50 -Ou --threads 12 -f $DIR/final_genome/C_striata_01.fasta $DIR/PSMC/C_striata_01_sorted.bam | bcftools call -c --threads 12 | vcfutils.pl vcf2fq -d 8 -D 50 | gzip> $DIR/PSMC/diploid_C_striata_01_8_50.fq.gz
```
Run time 2-08:28:18

*Parameters:*  
`mpileup` multi-way pileup producing genotype likelihoods  
`-C50` One may consider to add -C50 to mpileup if mapping quality is
overestimated for reads containing excessive mismatches. Applying this
option usually helps for BWA-backtrack alignments, but may not other
aligners.  
`-f`, `--fasta-ref` FILE  
`-c`, `--consensus-caller`  
`-d` sets and minimum read depth  
`-D` sets the maximum read depth  
\*It is recommended to set `-d` to a third of your average depth and
`-D` to twice the average depth (which can be seen from the GenomeScope
results)

#### 3\. Run PSMC analysis

Here we are running our PSMC analysis on the original diploid genome you
have.

First, convert your diploid.fastq file into a psmcfa file

``` bash
fq2psmcfa -q20 $DIR/PSMC/diploid_C_striata_01_8_50.fq.gz > $DIR/PSMC/diploid_C_striata_01_8_50.psmcfa
```
Run time 00:00:16

Now run the PSMC

``` bash
psmc -N30 -t30 -r5 -p "4+30*2+4+6+10" -o $DIR/PSMC/diploid_C_striata_01_8_50.psmc $DIR/PSMC/diploid_C_striata_01_8_50.psmcfa 
psmc -N30 -t30 -r5 -p "4*4+13*2+4*4+6" -o $DIR/PSMC/diploid_C_striata_01_8_50.psmc $DIR/PSMC/diploid_C_striata_01_8_50.psmcfa
```
Run time 01:33:26

*PSMC parameters:*  
`-p` STR pattern of parameters \[4+5\*3+4\]  
`-t` FLOAT maximum 2N0 coalescent time \[15\]  
`-N` INT maximum number of iterations \[30\]  
`-r` FLOAT initial theta/rho ratio \[4\]  
`-o` FILE output file \[stdout\]

#### 4\. Plot PSMC results

Before plotting your results you will need to know the mutation rate
`-u` and generation time `-g` for your specific species.  
Generation time is defined as the age at which 50% of the population has
reproduced. *Centropristis striata* is a protogynous species that has been
estimated to reproduce as a female around 4 years of age while males
tend to reproduce when they are 6 years old. Therefore we chose to set
generation time `g` at 5 years old.  
Additionally, mutation rate `-u` is a hard parameter to precisely
predict, thus we decided it would be best to show a range between
1x10<sup>-8</sup> and 1x10<sup>-9</sup>. Which is why we went ahead and
plotted both scenarios.

considering a mutation rate`-u` of 1x10<sup>-9</sup>:

``` bash
psmc_plot.pl -u 1e-09 -g 5 C_striata_01_8_50_t30r5_plot_u1-9g5 $DIR/PSMC/diploid_C_striata_01_8_50.psmc  
```
Run time 00:00:01

considering a mutation rate`-u` of 1x10<sup>-9</sup>:

``` bash
psmc_plot.pl -u 1e-08 -g 5 C_striata_01_8_50_t30r5_plot_u1-8g5 $DIR/PSMC/diploid_C_striata_01_8_50.psmc  
```
Run time 00:00:01

This last command will output two files .eps and .par. Copy both files
to your local computer. Open the .eps file to view your PSMC plot I like
to use Cyberduck to transfer (sftp) documents between the server and my
computer. You might notice that your PSMC plot will start at 10,000
years. This is because the PSMC analysis is better to look at population
size at older time scales. If you want to look at more recent time
scales you might consider running a MSMC analysis instead.

**mutation rate = 1x10<sup>-8</sup>**

<img width="635" alt="image" src="https://user-images.githubusercontent.com/26288352/178328690-e536c78e-dde9-4ea3-bb80-1129c7f00789.png">

**mutation rate = 1x10<sup>-9</sup>**

<img width="640" alt="image" src="https://user-images.githubusercontent.com/26288352/178328177-22d71a1b-b805-4d37-ac05-d3b310b18fcc.png">


## Run PSMC with using Barth et al 2017 parameters
#### 3\. Run PSMC analysis

Convert your diploid.fastq file into a psmcfa file

``` bash
fq2psmcfa -q20 -g10000 -s20 $DIR/PSMC/diploid_C_striata_01_8_50.fq.gz > $DIR/PSMC/diploid_C_striata_01_8_50_v2.psmcfa
```

`-q` minimum quality threshold
`-g` scaffold-good-size
`-s` window size

Run time 00:00:16

Now run the PSMC

``` bash
psmc -N25 -t30 -r5 -p "4*4+13*2+4*4+6" -o $DIR/PSMC/diploid_C_striata_01_8_50_v2.psmc $DIR/PSMC/diploid_C_striata_01_8_50_v2.psmcfa
```
Run time 01:33:26

*PSMC parameters:*  
`-p` STR pattern of parameters \[4+5\*3+4\]  
`-t` FLOAT maximum 2N0 coalescent time \[15\]  
`-N` INT maximum number of iterations \[30\]  
`-r` FLOAT initial theta/rho ratio \[4\]  
`-o` FILE output file \[stdout\]

#### 4\. Plot PSMC results

considering a mutation rate`-u` of 1x10<sup>-9</sup>:

``` bash
psmc_plot.pl -u 1e-09 -g 5 C_striata_01_8_50_t30r5_plot_u1-9g5_v2 $DIR/PSMC/diploid_C_striata_01_8_50_v2.psmc  
```
Run time 00:00:01

considering a mutation rate`-u` of 1x10<sup>-9</sup>:

``` bash
psmc_plot.pl -u 1e-08 -g 5 C_striata_01_8_50_t30r5_plot_u1-8g5_v2 $DIR/PSMC/diploid_C_striata_01_8_50_v2.psmc 
```
Run time 00:00:01

**mutation rate = 1x10<sup>-8</sup>**

![image](https://user-images.githubusercontent.com/26288352/178598810-c6e4d180-02c5-4bc1-a1a5-3d2ae7da1e6e.png)


**mutation rate = 1x10<sup>-9</sup>**

![image](https://user-images.githubusercontent.com/26288352/178598863-f685b1d8-2f93-45f3-a0a0-2564dbd7515a.png)


## PSMC with Bootstrap

Overall the steps to run a PSMC with bootstrap are the very similar to
those above, with the difference that you will need to run a prior step
to generate a split file and will need to run a job array to
simultaneously run the bootstrap PSMC analyses. Running a job array can
be different depending on your server’s task management. The steps below
illustrate the steps to run a job array on SLURM.

**Steps 1, 2, and 3 are exactly the same as the steps explained above. If
you have already done these, skip to step 4.**

#### 1\. Prepare your genome assembly data

this step you might already have from your first PSMC analysis howe

First we need to index our genome assembly (should be a fasta file)

``` bash
bwa index HPA_assembly.fa 
```

Now map adapter trimmed Illumina reads to your genome:

``` bash
bwa mem -t32 HPA_assembly.fasta HPA_HiSeq_R1.fq.gz HPA_HiSeq_R2.fq.gz > HPA_bwa_aligned.sam
```

Use samtools to convert your sam file to a bam file:

``` bash
samtools view -Sb -@ 30 -O BAM -o HPA_bwa_aligned.bam HPA_bwa_aligned.sam 
```

Sort your bam file:

``` bash
samtools sort -o HPA_bwa_aligned_sorted.bam -O BAM -@ 20 HPA_bwa_aligned.bam 
```

Finally, index your sorted bam file:

``` bash
samtools index -b -@ 20 HPA_bwa_aligned_sorted.bam  
```

#### 2\. Call diploid genome

Once you have mapped and sorted your genome the next step is to call a
diploid genome and convert it to the required format to input into the
PSMC analysis.

So, first we call a diploid genome using our reference genome and mapped
reads using bcftools. \*DeLeonLab– NOTE: I installed bcftools within the
samtools conda environment on the chimera server

``` bash
bcftools mpileup -C50 -Ou --threads 12 -f HPA_pilon.fasta HPA_bwa_aligned_sorted.bam | bcftools call -c --threads 10 | vcfutils.pl vcf2fq -d 35 -D 220 | gzip> diploid_HPA_35_220.fq.gz
```

*Parameters:*  
`mpileup` multi-way pileup producing genotype likelihoods  
`-C50` One may consider to add -C50 to mpileup if mapping quality is
overestimated for reads containing excessive mismatches. Applying this
option usually helps for BWA-backtrack alignments, but may not other
aligners.  
`-f`, `--fasta-ref` FILE  
`-c`, `--consensus-caller`  
`-d` sets and minimum read depth  
`-D` sets the maximum read depth  
\*It is recommended to set `-d` to a third of your average depth and
`-D` to twice the average depth (which can be seen from the GenomeScope
results)

#### 3\. Run PSMC analysis

Here we are running our PSMC analysis on the original diploid genome you
have.

First, convert your diploid.fastq file into a psmcfa file

``` bash
/hpcstor4/data01/DeLeonLab/apps/psmc/utils/fq2psmcfa -q20 diploid_HPA_35_220.fq.gz > diploid_HPA_35_220.psmcfa
```

Now run the PSMC

``` bash
/hpcstor4/data01/DeLeonLab/apps/psmc/psmc -N30 -t30 -r5 -p "4+30*2+4+6+10" -o diploid_HPA_35_220_t30r5.psmca.psmc diploid_HPA_35_220.psmcfa 
```

*PSMC options:*  
`-p` STR pattern of parameters \[4+5\*3+4\]  
`-t` FLOAT maximum 2N0 coalescent time \[15\]  
`-N` INT maximum number of iterations \[30\]  
`-r` FLOAT initial theta/rho ratio \[4\]  
`-o` FILE output file \[stdout\]

#### 4\. Generate split file for bootstrap

We sill use the splitfa command to split long chromosome sequences found
in your diploid.psmcfa file to shorter segments for bootstrapping.

``` bash
splitfa diploid_C_striata_01_8_50.psmcfa > diploid_C_striata_01_8_50_split.psmcfa
```

Once you have your diploid\_split.psmcfa file you will need to copy this
file into 100 independent files. I personally like to do this in a
separate directory, so:

``` bash
mkdir bootstrap
```

now copy(or move) your diploid\_split.psmcfa and your original psmc run
outfile into your new bootstrap directory. The psmc file will be used
after you run the bootstrap to concatenate with the other output files.

``` bash
cp diploid_C_striata_01_8_50_split.psmcfa bootstrap
cp diploid_C_striata_01_8_50.psmc bootstrap
```

Now move into your bootstrap directory

``` bash
cd bootstrap
```

In order to run a job array we need to generate input files numbered
from 1-100. I like to do this using an interactive session.

To run an interactive session in a SLURM environment (e.g. Chimera
server):

``` bash
srun -p short -N 1 --pty /bin/bash
```

Once in the interactive session, copy your diploid\_split.psmcfa file
into 100 different fileas labelled from 001-100.

``` bash
echo split_CST_{001..100}.psmcfa| xargs -n 1 cp diploid_C_striata_01_8_50_split.psmcfa
```

Exit interactive mode

``` bash
exit
```

#### 5\. Run PSMC in bootstrap mode using a job array

Now you need to create a job array that will submit multiple jobs at
once. For the SLURM manager on the Chimera server you can just create a
new document called psmc\_array.sh using nano or the text editor of your
choice:

``` bash
nano psmc_array.sh
```

Now copy the following SLURM script into your psmc\_array.sh \*You will
need to edit your job submission parameters to correspond to your
specific server as well as your specific input file name of your 100
input files.

``` bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=psmc_array              # Name your job something useful for easy tracking
#SBATCH --output=out/array_%A_%a.out
#SBATCH --error=out/array_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=6000MB
#SBATCH --time=00-04:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=short
#SBATCH --array=1-101%20		#run a maximum of 20 jobs at a time


# ----------------Modules------------------------- #
module load miniconda3
source activate /work/lotterhos/programs/PSMC

# ----------------Your Commands------------------- #
#
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# select our filename
N=${SLURM_ARRAY_TASK_ID}
# Comment one of the following two lines, depending on if the file names have leading zeros
#FILENAME=run-${N} # without leading zeros
 FILENAME=split_CST_$(printf "%03d" ${N}) # with leading zeros
      # adjust "%03d" to as many digits as are in the numeric part of the file name

echo "My input file is ${FILENAME}"

#
echo $P
#
psmc -N30 -t30 -r5 -b -p "4+30*2+4+6+10" -o ${FILENAME}.psmc ${FILENAME}.psmcfa
#

echo "Job finished" `date`
echo "My input file is ${FILENAME}"
```

The `-b` found in your psmc command indicates it is in bootstrap mode.

Now submit your job array as you would submit any other job.

``` bash
sbatch psmc_array.sh
```

Once your job array is finished you should have 100 different .psmc
files in addition to your original psmc analysis you ran in step 3.

We need to concatenate all .psmc files:

``` bash
cat *.psmc > HPA_35_220_combined.psmc
```

#### 6\. Plot PSMC results

The final step is to plot your combined results. Remember you will need
to change the mutation rate and generation time accordingly to your
species.

For *Holacanthus passer* I used a generation time `-g` of 5 and plotted
using a mutation rate `-u` for 1x10<sup>-8</sup> and 1x10<sup>-9</sup>

``` bash
/hpcstor4/data01/DeLeonLab/apps/psmc/utils/psmc_plot.pl -u 1e-08 -g 5 HPA_35_220_t30r5_plot_u1-8g5 diploid_HPA_35_220_t30r5.psmc  
```

``` bash
/hpcstor4/data01/DeLeonLab/apps/psmc/utils/psmc_plot.pl -u 1e-09 -g 5 HPA_35_220_t30r5_plot_u1-9g5 diploid_HPA_35_220_t30r5.psmc  
```

mutation rate = 1x10<sup>-8</sup>
<p align="center">
<img src="images/HPA_35_220_t30r5_plot_u1-8g5_boot.png" width="600"/>
</p>

mutation rate = 1x10<sup>-9</sup>
<p align="center">
<img src="images/HPA_35_220_t30r5_plot_u1-9g5_boot.png" width="600"/>
</p>
