# STACKS

## file sequence check
Samples have already been demultiplexed previously (see Thais' [pipeline.md](https://github.com/thais-neu/BlackSeaBass_project/blob/master/BSB_ddRAD/pipeline.md)) so I am going to go ahead and check `.fq.gz` files first

I have installed `multiqc` using conda from within the lotterhos module

```bash
module load lotterhos/2020-08-24
conda create -n multiqc
source activate multiqc
conda install -c bioconda multiqc
```
I ran `multiqc` on trimmed and synced files from `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed`

Thais has a separate directory where she has the same files but did not trim the adapter sequence to use for dDocent:  
`/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed`

for now I will use the original trimmed sequences and check the length. Stacks needs sequences to all be trimmed to the same length so I will check there is no contamination as well as sequence length to trim if needed. 

To run interactive jobs on the lotterhos partition:
lotterhos - 2 nodes - 36 cores 

`srun -p lotterhos -N 1 --pty /bin/bash`

OR For 20 minutes on the debug partition:

`srun -p debug -N 1 --pty /bin/bash`

We first need to run fastqc on all files and then we use multiqc to view all reports in one. 


```bash
srun -p lotterhos -N 1 --pty /bin/bash
module load oracle_java
module load fastqc
```
now run fastqc for all files

    fastqc /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/*fq.gz -o .

now load multiqc

      module load lotterhos
      source activate multiqc
run

      mulitqc .

now download the html file to your computer to open. I am interested in sequence length to be able to trim them all equal

![sequence length summmary](/img/multiqc_synced_trimmed_sequence_length.jpg)
![adapter presence](/img/multiqc_synced_trimmed_adapter.png)

The shortest sequence is 138 so I will trim to 138 and remove any leftover nexterra adapters.

## trimgalore test
```bash
module load lotterhos
source activate trimgalore

trim_galore --phred33 --fastqc --nextera -o ../samples --paired --cores 2 /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.F.fq.gz /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.R.fq.gz
```
this gave me an error that may be due to overlapping programs between trimgalore and salmon
*I also re-installed trim-galore as a conda environment because trimgalore from the lotterhos module could not find cutadapt*

The following worked:
```bash
module load miniconda3
source activate trimgalore

trim_galore --phred33 --fastqc --nextera -o ../samples/130bp --paired --cores 2 --length 100 --hardtrim5 130 --basename trim130 /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.F.fq.gz /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.R.fq.gz
```
when I use `--hardtrim` this takes precedence over everything else and only trims all sequences. So, I will need to remove the nexterra adapter and then run hardtrim first and filter out any smaller sequences 

```bash
module load miniconda3
source activate trimgalore

trim_galore --phred33 --fastqc --nextera -o ../samples/130bp --paired --cores 2 --length 100 /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.F.fq.gz /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.R.fq.gz
```

hardtrim to 135
```bash
trim_galore --phred33 -o . --paired --cores 2 --hardtrim5 135 /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.F.fq.gz /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/ME_165.R.fq.gz
```
Then filter smaller fragements out (min length 135)
```bash
trim_galore --phred33 --fastqc -o . --paired --cores 2 --length 135 --basename ME_156_trim135 /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/130bp/ME_165.F.135bp_5prime.fq.gz /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/130bp/ME_165.R.135bp_5prime.fq.gz
```


I plan to use an array so I first made a list of the unique sample names
```bash
head BSB_sample_list_uniq 
MA_298
MA_299
MA_300
MA_302
MA_303
MA_304
MA_306
MA_307
MA_310
...
```

## trim galore array

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=trim_galore              # Name your job something useful for easy tracking
#SBATCH --output=out/trim_galore.out
#SBATCH --error=out/trim_galore.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
#SBATCH --array=0-117%10		#there are 118 samples and it will run a maximum of 10 jobs at a time

#--------------MODULES---------------

module load miniconda3
source activate trimgalore

#--------------COMMAND----------------

sample=($(<../samples/BSB_sample_list_uniq))

echo "My input file is ${sample}"

trim_galore --phred33 -o ../samples/hardtrim --paired --cores 2 --hardtrim5 135 \
/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${sample}.F.fq.gz \
/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${sample}.R.fq.gz
```
This did not work. It kept repeating th first file over and over....

## trim galore for loop to trim to 135bp

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=trimgal              # Name your job something useful for easy tracking
#SBATCH --output=out/trimgal.out
#SBATCH --error=out/trim_gal.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos

#--------------MODULES---------------

module load miniconda3
source activate trimgalore

#--------------COMMAND----------------

for file in `cat BSB_sample_list_uniq`;
do
    trim_galore --phred33 -o ../samples/hardtrim --paired --cores 2 --hardtrim5 135 \
    /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${file}.F.fq.gz \
    /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${file}.R.fq.gz
done
```
or

```bash
samples=("MA_298
    MA_299
    MA_300
    MA_302
    MA_303
    MA_304
    MA_306
    MA_307
    MA_310
    MA_311")
     
for file in $samples;
do
    trim_galore --phred33 -o ../samples/hardtrim --paired --cores 2 --hardtrim5 135 \
    /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${file}.F.fq.gz \
    /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${file}.R.fq.gz
done
```
# Check fastqc stats

```bash
module load oracle_java
module load fastqc
```
now run fastqc for all files

    fastqc ../samples/hardtrim/*fq.gz -o ../samples/hardtrim/fastqc

now load multiqc

fist open an interactive node

       srun -p lotterhos -N 1 --pty /bin/bash

load multiqc

      module load lotterhos
      source activate multiqc
run

      mulitqc .

now download the html file to your computer to open. I am interested in sequence length to be able to trim them all equal

![multiqc_hardtrim](/img/multiqc_hardtrim_135bp.png)

samples still had adapters and some sequence lengths varied so we will first need to remove the adapter sequences and drop sequences under 130 bp. Then I will hardtrim all sequences.



## remove nexterra adapter and sequences smaller than 120bp for loop

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=trimgal              # Name your job something useful for easy tracking
#SBATCH --output=out/trimgal.out
#SBATCH --error=out/trim_gal.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos

#--------------MODULES---------------

module load miniconda3
source activate trimgalore

#--------------COMMAND----------------

for file in `cat ../samples/BSB_sample_list_uniq`;
do
    trim_galore --phred33 --fastqc --nextera -o ../samples/no_adapter --paired --cores 2 --length 130  \
    /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${file}.F.fq.gz \
    /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/${file}.R.fq.gz
done
```
Run time 05:48:11 (118 individuals)

Output notes:
The `*_trimmed.fq.gz` are produced as intermediate output (as R1 and R2 are trimmed individually in the first instance). Once the trimming has completed, Trim Galore will launch a round of 'validation' (which is is where the files get the val in their names from), which primarily performs length-cutoff filtering (and a few more optional things I believe). Once the validation is complete, the trimmed files will be deleted, and you are left with only the files `N1_1_val_1.fq.gz` and `N1_2_val_2.fq.gz`.

now load multiqc

fist open an interactive node

       srun -p lotterhos -N 1 --pty /bin/bash

load multiqc

      module load lotterhos
      source activate multiqc
run

      mulitqc .


## trim galore for loop to trim to 135bp

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=trimgal              # Name your job something useful for easy tracking
#SBATCH --output=out/trimgal.out
#SBATCH --error=out/trim_gal.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos

#--------------MODULES---------------

module load miniconda3
source activate trimgalore

#--------------COMMAND----------------

for file in `cat ../samples/BSB_sample_list_uniq`;
do
    trim_galore --phred33 -o ../samples/no_adapter/hardtrim --paired --cores 2 --hardtrim5 130 \
    ../samples/no_adapter/${file}.F_val_1.fq.gz \
    ../samples/no_adapter/${file}.R_val_2.fq.gz
done
```

