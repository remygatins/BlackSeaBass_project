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

```bash
multiqc -d /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/trim_padding/synced_renamed/*.fq.gz -l ../multiqc/BSB_synced_trimmed_list -i BSB_synced_trimmed -o .
```


