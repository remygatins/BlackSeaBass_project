# BSB RADseq STACKS pipeline with a updated reference genome `C_striata_v2.fasta`

Paths:

Working directory: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2`
Reference genome:   `/work/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`
Trimmed sequences: `/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags` 

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
bwa index $DIR/C_striata_v2.fasta
```

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa



