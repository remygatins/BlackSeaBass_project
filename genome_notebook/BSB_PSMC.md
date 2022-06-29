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

``` bash
jellyfish count -C -m 21 -s 60000000000 -t 12 /work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/PacBio_Denovo/raw_sequences/DTG-DNA-1126.r64296e173242G01.subreads_ccs.fastq.gz -o BSB_PacBio.jf
```

Jellyfish can’t read gzipped files so if you have fq.gz files use the
following command:

``` bash
jellyfish count -C -m 21 -s 60000000000 -t 12 <(zcat /work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/PacBio_Denovo/raw_sequences/DTG-DNA-1126.r64296e173242G01.subreads_ccs.fastq.gz) -o BSB_PacBio.jf
```

*Parameter key:*  
`-C` indicates to count canonical kmers (**don’t change this**) 
`-m` kmer length (default -m 21)  
`-s` allocate memory  
`-t` threads  
`-o` output name (default is mer\_counts.jf)

3.  Now export the kmer count histogram

<!-- end list -->

``` bash
jellyfish histo -t 12 HPA_HiSeq.jf > HPA_HiSeq.histo
```

4.  Upload reads.histo to GenomeScope: <http://qb.cshl.edu/genomescope/>
<p align="center">
<img src="images/GenomeScope_profile.png" width="500"/>
</p>

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

#### 2\. Call diploid genome and covert to psmc file

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
reproduced. *Holacanthus passer* is a protogynous species that has been
estimated to reproduce as a female around 4 years of age while males
tend to reproduce when they are 6 years old. Therefore we chose to set
generation time `g` at 5 years old.  
Additionally, mutation rate `-u` is a hard parameter to precisely
predict, thus we decided it would be best to show a range between
1x10<sup>-8</sup> and 1x10<sup>-9</sup>. Which is why we went ahead and
plotted both scenarios.

considering a mutation rate`-u` of 1x10<sup>-8</sup>:

``` bash
/hpcstor4/data01/DeLeonLab/apps/psmc/utils/psmc_plot.pl -u 1e-09 -g 5 HPA_35_220_t30r5_plot_u1-9g5 diploid_HPA_35_220_t30r5.psmcfa.psmc  
```

considering a mutation rate`-u` of 1x10<sup>-9</sup>:

``` bash
/hpcstor4/data01/DeLeonLab/apps/psmc/utils/psmc_plot.pl -u 1e-09 -g 5 HPA_35_220_t30r5_plot_u1-9g5 diploid_HPA_35_220_t30r5.psmcfa.psmc  
```

This last command will output two files .eps and .par. Copy both files
to your local computer. Open the .eps file to view your PSMC plot I like
to use Cyberduck to transfer (sftp) documents between the server and my
computer. You might notice that your PSMC plot will start at 10,000
years. This is because the PSMC analysis is better to look at population
size at older time scales. If you want to look at more recent time
scales you might consider running a MSMC analysis instead.

mutation rate = 1x10<sup>-8</sup>
<p align="center">
<img src="images/HPA_35_220_t30r5_plot_u1-8g5.png" width="600"/>
</p>

mutation rate = 1x10<sup>-9</sup>
<p align="center">
<img src="images/HPA_35_220_t30r5_plot_u1-9g5.png" width="600"/>
</p>


## PSMC with Bootstrap

Overall the steps to run a PSMC with bootstrap are the very similar to
those above, with the difference that you will need to run a prior step
to generate a split file and will need to run a job array to
simultaneously run the bootstrap PSMC analyses. Running a job array can
be different depending on your server’s task management. The steps below
illustrate the steps to run a job array on SLURM.

Steps 1, 2, and 3 are exactly the same as the steps explained above. If
you have already done these, skip to step 4.

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
/hpcstor4/data01/DeLeonLab/apps/psmc/utils/splitfa diploid_HPA_35_220.psmcfa > diploid_HPA_35_220_split.psmcfa
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
cp diploid_HPA_35_220_split.psmcfa bootstrap
cp diploid_HPA_35_220_t30r5.psmca.psmc bootstrap
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
srun -N 1 -n 12 -p AMD6128 --mem=50000MB -t 04:00:00 --pty bash
```

Once in the interactive session, copy your diploid\_split.psmcfa file
into 100 different fileas labelled from 001-100.

``` bash
echo split_HPA_{001..100}.psmcfa| xargs -n 1 cp diploid_HPA_35_220_split.psmcfa
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
#
#
#SBATCH -p AMD6276   # Partition name
#SBATCH -J psmc_array  # Job name
#SBATCH --mail-user=user_email@umb.edu  #change this to your email
#SBATCH --mail-type=ALL
#SBATCH -o array_%A_%a.out    # Name of stdout output file
#SBATCH -e array_%A_%a.err    # Name of stdout output file
#SBATCH --array=1-100
#SBATCH -N 1        # Total number of nodes requested
#SBATCH -n 2        # Total number of mpi tasks requested per node
#SBATCH -t 03-24:00:00  # Run Time (DD-HH:MM:SS) - 1.5 hours (optional)
#SBATCH --mem=6000MB  # Memory to be allocated PER NODE (default: 1gb)


echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


# ----------------Modules------------------------- #
source ~/.bashrc
conda activate samtools
#
# ----------------Your Commands------------------- #
#
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# select our filename
N=${SLURM_ARRAY_TASK_ID}
# Comment one of the following two lines, depending on if the file names have leading zeros
#FILENAME=run-${N}.inp # without leading zeros
 FILENAME=split_HPA_$(printf "%03d" ${N}).psmcfa # with leading zeros
# adjust "%03d" to as many digits as are in the numeric part of the file name
echo "My input file is ${FILENAME}"

#
echo $P
#
/hpcstor4/data01/DeLeonLab/apps/psmc/psmc -N30 -t30 -r5 -b -p "4+30*2+4+6+10" -o /hpcstor4/data01/DeLeonLab/remy/HPA_genome/PSMC/bootstrap_35_220/${FILENAME}.psmc /hpcstor4/data01/DeLeonLab/remy/HPA_genome/PSMC/bootstrap_35_220/${FILENAME}
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
