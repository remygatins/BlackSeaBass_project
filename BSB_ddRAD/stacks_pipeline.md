# BSB RADseq STACKS pipeline

- [check sequences](#file-sequence-check)
- [Filter adapters](#trim-adapters)
- [STACKS](#stacks)
    - [Process_radtags](#process-radtags) 
    - [Denovo_map test run](#denovo-map)
    - [optimize parameters](#optimize-parameters)
    - [Populations](#populations)
        - [BSB_r0.8_R0.8_Filt](#BSB_r0.8_R0.8_Filt) 
        - [BSB_r0.8_R0.8_Filt_noNC](#BSB_r0.8_R0.8_Filt_noNC)   




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

## Trim adapters
### trimgalore test
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
fastqc

    fastqc ../samples/no_adapter/hardtrim/*fq.gz -o ../samples/no_adapter/hardtrim/fastqc
    
multiqc

     srun -p lotterhos -N 1 --pty /bin/bash

load multiqc

      module load lotterhos
      source activate multiqc
run

      mulitqc .
      
All sequences are trimmed to a total length of 130bp starting from 5' end and have no adapter. (There may be a possibility that these may no longer be paired perfectly since the 3' end may have been cut and the reciprocal end may not. will check this.

# STACKS

## Process radtags

```bash
module load lotterhos
module load stacks

process_radtags -p ../samples/no_adapter --paired -o ../samples/no_adapter/process_radtags --renz-1 bamHI --renz-2 mspI -c -q -t 130 

```
error

```bash

process_radtags -1 MA_299.F_val_1.fq.gz -2 MA_299.R_val_2.fq.gz \
-o /work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags \
--renz-1 bamHI --renz-2 mspI -c -q -t 130
```

output
```bash
3670680 total sequences
      0 barcode not found drops (0.0%)
      0 low quality read drops (0.0%)
3670680 RAD cutsite not found drops (100.0%)
      0 retained reads (0.0%)
 ```
 Corresponding enzyme and read are flipped.  
 
 Read 1 has cutsite mspI while R2 has cutsite bamHI

```bash
 process_radtags -1 MA_299.F_val_1.fq.gz -2 MA_299.R_val_2.fq.gz \
 -o ./process_radtags \
 --renz-1 mspI --renz-2 bamHI \
 -c -q -t 130
 ```
 
 output
 ```bash
3670680 total sequences
      0 barcode not found drops (0.0%)
   1042 low quality read drops (0.0%)
      0 RAD cutsite not found drops (0.0%)
3669638 retained reads (100.0%)
 ```
 
 
 ```bash
#!/bin/bash
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks.out
#SBATCH --error=out/stacks.err
#SBATCH --cpus-per-task=14
#SBATCH --mem=5000                        # Allocate 5GB of RAM.  You must declare --mem in all scripts
#SBATCH --time=2-24:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu	   # replace "cruzid" with your user id
#SBATCH --mail-type=ALL                   # Only send emails when jobs end or fail
#SBATCH --partition=lotterhos
#--------------MODULES---------------

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
```
output 
```bash
SN_191.F_val_1.1.fq.gz
SN_191.F_val_1.rem.1.fq.gz
SN_191.R_val_2.2.fq.gz
SN_191.R_val_2.rem.2.fq.gz
```
`*.rem.1.fq` and `*.rem.2.fq` are remaining sequences that do not have a pair due to being discaded due to quality.

I want to rename my sequences to make them more straightforward. 


`for f in *.F_val_1.1.fq.gz; do mv "$f" "${f%.F_val_1.1.fq.gz}.1.fq.gz"; done`

`for f in *.R_val_2.2.fq.gz; do mv "$f" "${f%.R_val_2.2.fq.gz}.2.fq.gz"; done`

`for f in *.F_val_1.rem.1.fq.gz; do mv "$f" "${f%.F_val_1.rem.1.fq.gz}.rem.1.fq.gz"; done`

`for f in *.R_val_2.rem.2.fq.gz; do mv "$f" "${f%.R_val_2.rem.2.fq.gz}.rem.2.fq.gz"; done`

```bash
...
SN_191.1.fq.gz
SN_191.2.fq.gz
SN_191.rem.1.fq.gz
SN_191.rem.2.fq.gz
```
## Denovo map

```bash
denovo_map.pl -m 4 -M 3 -n 2 -T 12 -d \
-o ../stacks/ \
--samples ../samples/no_adapter/process_radtags \
--popmap ../popmap/BSB_all --paired \
-X "populations:-r 0.8 --genepop --vcf --write_single_snp"
```
I removed MA_300 due to low number of sequences. 

denovo_map.log
```bash
#cstacks
MA_298	MA	42.60x
MA_299	MA	43.40x
MA_302	MA	34.98x
MA_303	MA	21.75x
MA_304	MA	29.09x
MA_306	MA	42.18x
MA_307	MA	36.14x
MA_310	MA	38.08x
MA_311	MA	39.84x
MA_313	MA	33.30x
MA_314	MA	14.56x
MA_315	MA	37.25x
MA_316	MA	23.57x
MA_318	MA	25.48x
MA_320	MA	37.66x
MA_321	MA	42.68x
MA_323	MA	35.42x
MA_324	MA	8.62x
MA_325	MA	23.77x
MA_327	MA	17.87x
MD_136	MD	39.44x
MD_137	MD	44.32x
MD_138	MD	50.90x
MD_139	MD	43.49x
MD_140	MD	35.03x
MD_141	MD	43.29x
MD_142	MD	40.80x
MD_143	MD	37.12x
MD_145	MD	45.55x
MD_149	MD	43.84x
MD_150	MD	47.33x
MD_151	MD	43.61x
MD_152	MD	40.37x
MD_154	MD	37.23x
MD_158	MD	33.72x
MD_159	MD	32.40x
MD_160	MD	38.00x
MD_161	MD	42.58x
MD_162	MD	44.66x
MD_163	MD	31.44x
ME_164	ME	35.45x
ME_165	ME	39.00x
ME_166	ME	39.43x
ME_167	ME	37.77x
ME_176	ME	40.33x
ME_248	ME	45.17x
ME_249	ME	47.66x
ME_250	ME	25.59x
ME_251	ME	30.39x
ME_252	ME	36.61x
ME_253	ME	40.50x
ME_254	ME	42.78x
ME_255	ME	47.88x
ME_256	ME	44.72x
ME_257	ME	51.40x
ME_258	ME	44.83x
ME_261	ME	43.97x
ME_262	ME	44.23x
NC_233	NC	42.17x
NC_234	NC	43.52x
NC_235	NC	46.29x
NC_237	NC	44.17x
NC_238	NC	36.86x
NC_239	NC	31.68x
NC_240	NC	28.49x
NC_241	NC	46.20x
NC_242	NC	33.47x
NC_243	NC	30.66x
NC_244	NC	28.01x
NC_245	NC	31.70x
NC_246	NC	21.07x
NJ_106	NJ	28.75x
NJ_108	NJ	25.45x
NJ_109	NJ	31.31x
NJ_112	NJ	27.14x
NJ_113	NJ	38.94x
NJ_114	NJ	42.53x
NJ_118	NJ	40.87x
NJ_119	NJ	25.85x
NJ_121	NJ	25.47x
NJ_122	NJ	34.16x
NJ_124	NJ	34.01x
NJ_128	NJ	36.18x
NJ_129	NJ	43.95x
NJ_130	NJ	43.80x
NJ_131	NJ	43.46x
NJ_132	NJ	37.54x
NJ_133	NJ	37.62x
RI_328	RI	37.00x
RI_329	RI	36.30x
RI_330	RI	28.80x
RI_331	RI	45.57x
RI_332	RI	47.19x
RI_333	RI	43.55x
RI_334	RI	37.46x
RI_335	RI	30.51x
RI_336	RI	38.68x
RI_337	RI	38.39x
RI_338	RI	34.62x
RI_339	RI	43.61x
RI_340	RI	46.72x
RI_341	RI	49.04x
RI_342	RI	42.16x
RI_343	RI	33.73x
RI_344	RI	41.30x
RI_345	RI	41.89x
RI_346	RI	46.02x
RI_347	RI	43.36x
RI_348	RI	46.99x
RI_349	RI	36.99x
SN_009	SN	43.00x
SN_179	SN	37.80x
SN_182	SN	39.38x
SN_185	SN	40.53x
SN_189	SN	31.60x
SN_190	SN	45.76x
SN_191	SN	46.75x

#gstacks
Genotyped 243815 loci:
  effective per-sample coverage: mean=44.6x, stdev=11.0x, min=8.1x, max=66.1x
  mean number of sites per locus: 229.0
  a consistent phasing was found for 411986 of out 478031 (86.2%) diploid loci needing phasing

#populations
Removed 224945 loci that did not pass sample/population constraints from 243815 loci.
Kept 18870 loci, composed of 4428444 sites; 7984 of those sites were filtered, 18067 variant sites remained.
Number of loci with PE contig: 18870.00 (100.0%);
  Mean length of loci: 224.68bp (stderr 0.36);
Number of loci with SE/PE overlap: 10390.00 (55.1%);
  Mean length of overlapping loci: 200.42bp (stderr 0.30); mean overlap: 28.52bp (stderr 0.05);
Mean genotyped sites per locus: 230.25bp (stderr 0.34).
```
adegenet summary:
```r
 /// GENLIGHT OBJECT /////////

 // 117 genotypes,  18,067 binary SNPs, size: 7.3 Mb
 1087754 (51.46 %) missing data

 // Basic content
   @gen: list of 117 SNPbin

 // Optional content
   @ind.names:  117 individual labels
   @loc.names:  18067 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 7-22)
   @other: a list containing: elements without names 
```


```bash
populations -P ../stacks/ -M ../popmap/BSB_all -r 0.80 -p 7 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30
```
Adegenet summary:
```bash
/// GENLIGHT OBJECT /////////

 // 117 genotypes,  1,008 binary SNPs, size: 382.6 Kb
 8533 (7.24 %) missing data

 // Basic content
   @gen: list of 117 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  117 individual labels
   @loc.names:  1008 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 7-22)
   @other: a list containing: elements without names 

```

## Optimize parameters

`mkdir opt`

```bash
## Denovo map

```bash
src=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks
params="
2
3
4
5
6
7"

for p in $params
do
    denovo_map.pl -m 4 -M ${p} -n ${p} -T 12 -o $src/opt/M${p} --samples $src/samples/no_adapter/process_radtags --popmap $src/popmap/BSB_15x --paired -X "populations:-r 0.8 --write_single_snp"
done

```

**Options:**
`-m` — Minimum depth of coverage required to create a stack (default 3).  
`-M` — number of mismatches allowed between stacks within individuals (for ustacks).  
`-n` — number of mismatches allowed between stacks between individuals (for cstacks).  
`-samples [path]` — specify a path to the directory of samples (samples will be read from population map).  
`--popmap [path]` — path to a population map file (format is "[name] TAB [pop]", one sample per line).  
`-o [path]` — path to write pipeline output files.  
`-p`,`--min-populations` — minimum number of populations a locus must be present in to process a locus (for populations; default: 1). 
`-r`,`--min-samples-per-pop` — minimum percentage of individuals in a population required to process a locus for that population (for populations; default: 0). 

Summary table

`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|Opt    | Loci  |
|:------|------:|
|M2 	|18824  |
|**M3**    |**19018**	|
|M4		|18948  |
|M5		|18898  |
|M6     |18686  |
|M7     |18503  |

M3 obtained the most loci and is the optimum parameter.

Copy catalog and stacks from `opt/M3` to `stacks`

`cp ../opt/M3/* .`

## Populations

```bash
populations -P ../stacks/ -M ../popmap/BSB_15x -r 0.80 -p 7 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30
populations -P ../stacks/ -M ../popmap/BSB_15x -r 0.80 -p 7 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30
populations -P ../stacks/ -M ../popmap/BSB_15x -r 0.80 -p 6 --min-maf 0.05 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30
populations -P ../stacks/ -M ../popmap/BSB_15x -r 0.80 -p 6 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30
```

`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|populations        | Loci   |% missing data|
|:------------------|:------:|-------------:|
|r0.8 	            |18559   |   |*not sure why I obtained less loci after running populations again on original catalog that from denovo
|r0.8_p7 	        |1703    | |
|r0.8_p7_minmaf0.05 |1703    | |
|r0.8_p6 	        |5063    | |
|r0.8_p6_minmaf0.05 |5063    | |
|r0.8_R0.8_filt     |6167    | |

From Longo *et al* 2021:  

"We then used populations to remove loci that failed to meet the following criteria: present in ≥80% of individuals, minor allele frequency ≥1%, and maximum observed heterozygosity of 70%... We exported the resulting SNP dataset from Stacks and further filtered using VCFtools v.0.1.13 (Danecek et al.,   2011). We then dropped all but the first SNP from each RADseq locus (--thin 5000), removed loci in individuals that were below 10x depth of coverage (--minDP 10), refiltered for loci found in ≥80% of individuals (--max-missing 0.8), and then removed individuals with >30% missing loci from the final dataset (--remove), which was exported for downstream analyses"

In order to keep more loci I will use the r0.8 output that has the most loci and filter with vcftools

vcftools --vcf ../stacks/BSB_pop_r0.8/populations.snps.vcf --out ../stacks/BSB_pop_r0.8/vcftools/BSB_min10x_maxmiss0.8 --minDP 10 --max-missing 0.8 --recode --recode-INFO-all

vcftools --vcf populations.snps.vcf --minDP 10 --max-missing 0.8 --recode --recode-INFO-all

vcftools --vcf populations.snps.vcf --missing-indv

`cat remove_ind_0.40`

```bash
ME_253
ME_254
MA_303
MA_314
MA_318
MA_316
```

vcftools --vcf populations.snps.vcf --minDP 10 --max-missing 0.8 --remove remove_ind_0.40 --recode --recode-INFO-all --out filtered_ind

vcftools --vcf filtered_ind.recode.vcf --missing-indv --out filtered_ind

vcftools --vcf out.recode.vcf --remove remove_ind_0.40 --recode --recode-INFO-all --out filtered_ind_2

vcftools --vcf filtered_ind_2.recode.vcf --missing-indv --out filtered_ind_2


### Filtered_ind_2
Stacks: r0.8
vcftools: --minDP 10 --max-missing 0.8 --remove remove_ind_0.40
*removed individuals with more than 40% missing data

```bash
/// GENLIGHT OBJECT /////////

 // 110 genotypes,  27,125 binary SNPs, size: 3.9 Mb
 213457 (7.15 %) missing data

 // Basic content
   @gen: list of 110 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  110 individual labels
   @loc.names:  27125 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 7-22)
   @other: a list containing: elements without names 
```
![PCA](/img/PCA.png)
![dapc](/img/dapc_60pc.png)
![compoplot](/img/compoplot_60pc.png)

COnvert vcf to structure using radiator in R

```r
#--------------------------------
# COnvert vcf file to structure
#--------------------------------
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("thierrygosselin/radiator")
library(radiator)

#To verify your file is detected by radiator as the correct format:
radiator::detect_genomic_format(
  data="filtered_ind_2.recode.vcf")

strata <- radiator::read_strata(strata = "BSB.strata.filt.txt")

tidy.vcf.data <-  tidy_vcf(data="filtered_ind_2.recode.vcf",
                       strata = "BSB.strata.filt.txt")

#rename colummn GT_BIN to GT
library(dplyr)
tidy.vcf.data_2<- tidy.vcf.data %>% 
                      rename(GT= GT_BIN)

write_structure(tidy.vcf.data_2, filename = "str_r0.8_filtered")
```

## Structure

src=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks
iteration="
1
2
3"

for i in $iteration
do
    structure -K 1 -L 1 -N 110 -i $src/stacks/BSB_pop_r0.8/str_r0.8_filtered.str -o $src/structure/results/BSB_r0.8_filt_K1_${i}
done

structure -K 4 -L 27125 -N 110 -i $src/stacks/BSB_pop_r0.8/str_r0.8_filtered.str -o $src/structure/results/BSB_r0.8_filt_K4_1

error
```bash
# Entries:   Line numbers
     17430:   1
     34862:   2--111
----------------------------------
```
when converting from vcf to structure loci got dropped... ?


## BSB_r0.8_R0.8_Filt
Re-run populations
use filtered popmap `BSB_15x_filt`
keep loci found in 80% of all individuals `-R 0.80`

`populations -P ../stacks/ -M ../popmap/BSB_15x_filt -r 0.80 -R 0.80 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30`

```bash
/// GENLIGHT OBJECT /////////

 // 110 genotypes,  6,116 binary SNPs, size: 1.3 Mb
 48092 (7.15 %) missing data

 // Basic content
   @gen: list of 110 SNPbin

 // Optional content
   @ind.names:  110 individual labels
   @loc.names:  6116 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @other: a list containing: elements without names 
```


## Structure

`structure -K 4 -L 6116 -N 110 -i $src/stacks/BSB_r0.8_R0.8_filt/populations_str.gdv -o $src/structure/results/BSB_r0.8_R0.8_filt_K4_1`


```bash
#--------------MODULES---------------
module load lotterhos
source activate structure

#--------------COMMAND----------------
src=/work/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks
iteration="
1
2
3"

for i in $iteration
do
    structure -K 1 -L 6116 -N 110 -i $src/stacks/BSB_r0.8_R0.8_filt/populations_str.gdv -o $src/structure/results/BSB_r0.8_R0.8_filt_K1_${i}
done
```
![pca](/BSB_ddRAD/img/pca.png)

![structure](/BSB_ddRAD/img/K2_K3_K4.png)




## _BSB_r0.8_R0.8_Filt_noNC_
remove NC samples
use filtered popmap `BSB_15x_filt`
keep loci found in 80% of all individuals `-R 0.80`

`populations -P ../stacks/ -M ../popmap/BSB_15x_filt_noNC -r 0.80 -R 0.80 --write-single-snp --vcf --genepop --structure --fstats --hwe -t 30`

```r
 /// GENLIGHT OBJECT /////////

 // 97 genotypes,  5,639 binary SNPs, size: 1.1 Mb
 30828 (5.64 %) missing data

 // Basic content
   @gen: list of 97 SNPbin

 // Optional content
   @ind.names:  97 individual labels
   @loc.names:  5639 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @other: a list containing: elements without names
```

![pca_noNC](/BSB_ddRAD/img/pca_noNC.png)
![dapc_noNC](/BSB_ddRAD/img/dapc_noNC.png)
