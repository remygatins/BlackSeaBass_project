# BSB RADseq STACKS pipeline with an updated reference genome `C_striata_v2.fasta`

Paths:

Working directory: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2`
Reference genome:   `/projects/lotterhos/2021_BlackSeaBass_genomics/BSB_genome/final_genome/C_striata_v2.fasta`  
Trimmed sequences: `/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks/samples/no_adapter/process_radtags` 


## Run Ref_map
`cat stacks_2.sh`

```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
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


## Now run Populations
Within the `populations` directory, make the following directories

`mkdir p1_maf_0.01  p1_maf_0.05  p6_maf_0.01  p6_maf_0.05`

Now run the different populations

```bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=stacks              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short

export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2


populations \
  -P $WOR_DIR/stacks/ref_map_1_stacks2.6 \
  --popmap --popmap $WOR_DIR/popmap/BSB_all \
  -r 0.80 \
  -p 1 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O $WOR_DIR/populations/p1_maf_0.01 \
  -t 30
```

To obtain number of loci from each populations run
`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|populations|No. Loci|
|-----------|--------|
|maf 0.01|27664| 
|p1_maf_0.01|28707|
|p1_maf_0.05|28509|
|p6_maf_0.01|12693|
|p6_maf_0.05|12692|

I'm going with p1 maf 0.01 since it gives us the most loci and we will be doing more filtering with vcftools which may remove those extra loci anyway

## Filter VCF with VCFTOOLS


```bash
#!/bin/bash
#--------------SLURM COMMANDS--------------
#SBATCH --job-name=vcftools              # Name your job something useful for easy tracking
#SBATCH --output=out/stacks_%A.out
#SBATCH --error=out/stacks_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=48:00:00                 # Limit run to N hours max (prevent jobs from wedging in the queues)
#SBATCH --mail-user=r.gatins@northeastern.edu      # replace "cruzid" with your user id
#SBATCH --mail-type=BEGIN,END,FAIL                   # Only send emails when jobs end or fail
#SBATCH --partition=short

#load program
module load vcftools # this is just globally available on explorer

# paths
INPUT_VCF="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01/populations.snps.vcf"
OUTDIR="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01/"

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10
echo "done minimum depth"

# filter for sites present in >= 80% of individuals (refiltered any SNPs found in ≤80% of individuals)
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_maxmiss0.8
echo "done max maxmiss"

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
echo "done filtering piepline"
```


output
```bash

# filter by minimum depth per genotype (minDP = 10)

After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 24993 out of a possible 24993 Sites
Run Time = 4.00 seconds

# filter for sites present in >= 80% of individuals (refiltered any SNPs found in ≤80% of individuals)
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 0 out of a possible 24993 Sites
No data left for analysis!
Run Time = 0.00 seconds
```
Filtering by 0.8 max missingness seems to stringent, so I tried running it with 0.5 instead and keep having the same issue.
I believe it has to do with the format of the file that may be missing something from the stacks output. I will filter this directly within the populations command using -R 0.8.

```bash
export PATH=/projects/gatins/programs_explorer/stacks_2.68/bin:$PATH

WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2

populations \
  -P $WOR_DIR/stacks/ref_map_1_stacks2.6 \
  --popmap $WOR_DIR/popmap/BSB_all \
  -r 0.80 \
  -p 1 \
  --min-maf 0.01 \
  --R 0.80 \
  --write-single-snp \
  --vcf \
  --genepop \
  --structure \
  --fstats \
  --hwe \
  -O $WOR_DIR/populations/p1_maf_0.01_R0.8 \
  -t 30
```

`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|populations|No. Loci|
|-----------|--------|
|p1_maf_0.01_R0.8|13692|

ok now I will filter again with vcftools but i'll skip the `max_missingness` as we did this already 


```bash
#load program
module load vcftools # this is just globally available on explorer

# paths
INPUT_VCF="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01_R0.8/populations.snps.vcf"
OUTDIR="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/populations/p1_maf_0.01_R0.8/"

# filter by minimum depth per genotype (minDP = 10)
vcftools --vcf ${INPUT_VCF} \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10
echo "done minimum depth"

# remove individuals with >40% missing data
# 1: compute missingness per individual
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --missing-indv \
         --out ${OUTDIR}/missingness

# 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' ${OUTDIR}/missingness.imiss > ${OUTDIR}/remove_individuals.txt

# now filter out individuals
vcftools --vcf ${OUTDIR}/minDP10.recode.vcf \
         --remove ${OUTDIR}/remove_individuals.txt \
         --recode --recode-INFO-all \
         --out ${OUTDIR}/minDP10_filtInd
echo "done filtering piepline"
```


## Fixed ref map coverage error!

I figured it out!! 

last run:

```bash
Genotyped 246427 loci:
  effective per-sample coverage: mean=1.1x, stdev=0.0x, min=1.0x, max=1.1x
  mean number of sites per locus: 217.4
  a consistent phasing was found for 366964 of out 367070 (100.0%) diploid loci needing phasing
```

Now, when running `ref_map.pl` with `--remove-duplicates` it was removing 95% of my data. So after rerunning ref_map.pl without removing pcrr duplicates I get:



```bash
Genotyped 509504 loci:
  effective per-sample coverage: mean=31.0x, stdev=9.5x, min=3.2x, max=48.5x
  mean number of sites per locus: 170.8
  a consistent phasing was found for 731471 of out 819475 (89.3%) diploid loci needing phasing
```

Depth per sample
```bash
BEGIN effective_coverages_per_sample
# For mean_cov_ns, the coverage at each locus is weighted by the number of
# samples present at that locus (i.e. coverage at shared loci counts more).
sample	n_loci	n_used_fw_reads	mean_cov	mean_cov_ns
MA_298_aligned_sorted	49196	1227023	24.942	35.874
MA_299_aligned_sorted	54067	1551141	28.689	41.496
MA_302_aligned_sorted	40144	644066	16.044	22.137
MA_303_aligned_sorted	18296	61267	3.349	3.713
MA_304_aligned_sorted	38002	457046	12.027	15.503
MA_306_aligned_sorted	53964	1334626	24.732	37.150
MA_307_aligned_sorted	50746	1100308	21.683	30.625
MA_310_aligned_sorted	49690	1009794	20.322	30.228
MA_311_aligned_sorted	56638	1078240	19.037	31.723
MA_313_aligned_sorted	51383	737043	14.344	19.881
MA_314_aligned_sorted	25789	87833	3.406	3.817
MA_315_aligned_sorted	44818	1005950	22.445	31.110
MA_316_aligned_sorted	26883	245526	9.133	10.890
MA_318_aligned_sorted	20581	209521	10.180	12.011
MA_320_aligned_sorted	52626	1205223	22.902	33.590
MA_321_aligned_sorted	51869	1379906	26.604	38.616
MA_323_aligned_sorted	53844	1206189	22.402	32.545
MA_324_aligned_sorted	21863	64637	2.956	3.180
MA_325_aligned_sorted	45226	779629	17.239	21.796
MA_327_aligned_sorted	41244	482932	11.709	14.425
MD_136_aligned_sorted	58129	1248460	21.477	34.321
MD_137_aligned_sorted	54077	997689	18.449	31.658
MD_138_aligned_sorted	57094	1111874	19.474	37.139
MD_139_aligned_sorted	54791	1404790	25.639	38.812
MD_140_aligned_sorted	54697	693328	12.676	20.002
MD_141_aligned_sorted	47645	878733	18.443	27.431
MD_142_aligned_sorted	54215	906150	16.714	28.855
MD_143_aligned_sorted	54385	674159	12.396	22.667
MD_145_aligned_sorted	56817	1308301	23.027	37.314
MD_149_aligned_sorted	54384	1147749	21.105	34.113
MD_150_aligned_sorted	59233	1326466	22.394	38.216
MD_151_aligned_sorted	56769	1293392	22.783	36.243
MD_152_aligned_sorted	58706	1151654	19.617	31.683
MD_154_aligned_sorted	45224	624826	13.816	22.382
MD_158_aligned_sorted	59246	694715	11.726	22.055
MD_159_aligned_sorted	54288	663968	12.230	21.727
MD_160_aligned_sorted	48211	790004	16.386	26.063
MD_161_aligned_sorted	56178	961180	17.110	28.203
MD_162_aligned_sorted	59466	1129729	18.998	33.836
MD_163_aligned_sorted	58395	657193	11.254	20.816
ME_164_aligned_sorted	52197	926052	17.741	26.295
ME_165_aligned_sorted	45776	885426	19.343	26.991
ME_166_aligned_sorted	47194	1346222	28.525	38.991
ME_167_aligned_sorted	47064	1129732	24.004	33.012
ME_176_aligned_sorted	51404	1208884	23.517	35.017
ME_248_aligned_sorted	55711	1446217	25.959	39.320
ME_249_aligned_sorted	54784	1715533	31.314	45.809
ME_250_aligned_sorted	71228	939663	13.192	22.330
ME_251_aligned_sorted	49829	662724	13.300	19.181
ME_252_aligned_sorted	46885	1095511	23.366	31.394
ME_253_aligned_sorted	10869	233815	21.512	18.551
ME_254_aligned_sorted	12267	236471	19.277	13.258
ME_255_aligned_sorted	51210	1612075	31.480	45.021
ME_256_aligned_sorted	53305	1353937	25.400	37.669
ME_257_aligned_sorted	53246	1768789	33.219	48.539
ME_258_aligned_sorted	48455	1444386	29.809	41.164
ME_261_aligned_sorted	59494	1524936	25.632	40.960
ME_262_aligned_sorted	54460	1453584	26.691	39.657
NC_233_aligned_sorted	52007	1203984	23.150	34.513
NC_234_aligned_sorted	50751	1251502	24.660	35.756
NC_235_aligned_sorted	49520	1464220	29.568	42.158
NC_237_aligned_sorted	52718	1408663	26.721	40.029
NC_238_aligned_sorted	49333	1107210	22.444	31.373
NC_239_aligned_sorted	44373	582902	13.136	17.651
NC_240_aligned_sorted	42641	687721	16.128	20.842
NC_241_aligned_sorted	50306	1397872	27.787	40.221
NC_242_aligned_sorted	43884	851001	19.392	25.326
NC_243_aligned_sorted	43236	841326	19.459	25.547
NC_244_aligned_sorted	46582	742954	15.949	21.623
NC_245_aligned_sorted	43379	923487	21.289	27.780
NC_246_aligned_sorted	40211	479803	11.932	14.572
NJ_106_aligned_sorted	40723	780810	19.174	24.033
NJ_108_aligned_sorted	41315	765054	18.518	23.240
NJ_109_aligned_sorted	44189	913300	20.668	27.287
NJ_112_aligned_sorted	41821	744258	17.796	22.830
NJ_113_aligned_sorted	46397	1188537	25.617	34.762
NJ_114_aligned_sorted	47892	1376453	28.741	39.858
NJ_118_aligned_sorted	46066	1286575	27.929	37.640
NJ_119_aligned_sorted	43104	716838	16.630	21.358
NJ_121_aligned_sorted	43368	624036	14.389	18.589
NJ_122_aligned_sorted	44990	1125148	25.009	32.990
NJ_124_aligned_sorted	46298	971089	20.975	28.434
NJ_128_aligned_sorted	44511	1097851	24.665	32.552
NJ_129_aligned_sorted	46669	1348910	28.904	39.759
NJ_130_aligned_sorted	47526	1367169	28.767	39.978
NJ_131_aligned_sorted	49085	1363568	27.780	39.283
NJ_132_aligned_sorted	52119	1233049	23.658	34.423
NJ_133_aligned_sorted	45557	1083192	23.777	32.028
RI_328_aligned_sorted	48258	1268864	26.293	35.589
RI_329_aligned_sorted	66362	1442899	21.743	35.645
RI_330_aligned_sorted	49928	897514	17.976	24.709
RI_331_aligned_sorted	47856	1454297	30.389	41.186
RI_332_aligned_sorted	52318	1569343	29.996	43.044
RI_333_aligned_sorted	47481	1297417	27.325	37.638
RI_334_aligned_sorted	49653	996724	20.074	28.125
RI_335_aligned_sorted	45383	854764	18.834	24.960
RI_336_aligned_sorted	60732	1214421	19.996	33.044
RI_337_aligned_sorted	47113	1168778	24.808	33.498
RI_338_aligned_sorted	44925	1042212	23.199	30.305
RI_339_aligned_sorted	59088	1580694	26.752	40.639
RI_340_aligned_sorted	57172	1655893	28.963	43.386
RI_341_aligned_sorted	52032	1715979	32.979	46.935
RI_342_aligned_sorted	55160	1451103	26.307	39.037
RI_343_aligned_sorted	49945	1021547	20.453	27.969
RI_344_aligned_sorted	47861	1430601	29.891	40.406
RI_345_aligned_sorted	45597	1253604	27.493	36.970
RI_346_aligned_sorted	50677	1473654	29.079	41.056
RI_347_aligned_sorted	46916	1323515	28.210	38.009
RI_348_aligned_sorted	48625	1495473	30.755	42.780
RI_349_aligned_sorted	52600	872651	16.590	28.173
SN_009_aligned_sorted	53901	1335188	24.771	37.300
SN_179_aligned_sorted	49474	1008188	20.378	29.073
SN_182_aligned_sorted	47850	1329646	27.788	38.010
SN_185_aligned_sorted	50626	1101822	21.764	32.518
SN_189_aligned_sorted	45009	1063334	23.625	30.412
SN_190_aligned_sorted	47927	1549012	32.320	44.197
SN_191_aligned_sorted	58432	1596787	27.327	43.849
END effective_coverages_per_sample
```

So I will use the output from the following ref map command
```bash
WOR_DIR=/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2

ref_map.pl -T 10 -o $WOR_DIR/stacks/ref_map_1_stacks2.6_2 --popmap $WOR_DIR/popmap/BSB_all --samples $WOR_DIR/samples -X "populations:-r 0.80 --min-maf 0.01 --fstats --vcf --genepop --hwe --structure`
```
Check number of Loci
`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

33203

Great! Now let's filter the vcf file


run interactive mode
```bash
srun -p short -N 1 --pty /bin/bash
module load vcftools
```

filter by minimum depth per genotype (minDP = 10)
```bash
vcftools --vcf populations.snps.vcf \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out minDP10
```

output
```bash
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 306227 out of a possible 306227 Sites
Run Time = 51.00 seconds
```
filter for sites present in >= 80% of individuals (refiltered any SNPs found in ≤80% of individuals)
```bash
vcftools --vcf minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out minDP10_maxmiss0.8
```
output
```bash
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 65631 out of a possible 306227 Sites
Run Time = 18.00 seconds
```

## Now remove individuals with >40% missing data

```bash
### 1: compute missingness per individual
vcftools --vcf minDP10_maxmiss0.8.recode.vcf \
         --missing-indv \
         --out missingness

### 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' missingness.imiss > remove_individuals.txt

### now filter out individuals
vcftools --vcf minDP10_maxmiss0.8.recode.vcf \
         --remove remove_individuals.txt \
         --recode --recode-INFO-all \
         --out minDP10_maxmiss0.8_filtInd
```

```bash
After filtering, kept 109 out of 117 Individuals
Outputting VCF file...
After filtering, kept 65631 out of a possible 65631 Sites
Run Time = 13.00 seconds
```
removed individuals
```bash
cat remove_individuals.txt
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

Lets run populations with different parameters to see if we can keep more individuals even if we have fewer loci

`cat populations.hapstats.tsv | grep -v "^#" | cut -f 1 | uniq | wc -l`

|Parameters|Loci|
|----------|----|
|maf 0.01_v2|33203|* this does not add up to the number of sites found in vcftools
|p_1_maf_0.01_v2|35278|
|p_1_maf_0.05_v2|35248|
|p_2_maf_0.01_v2|30565|
|p_3_maf_0.01_v2|27521|
|p_4_maf_0.01_v2|24883|
|p_5_maf_0.01_v2|20562|
|p_6_maf_0.01_v2|15438|
|p_6_maf_0.05_v2|15438|

**When running without the `-p` flag we don't get as many loci as with p1, however after filtering we have way more sites. Why??**

In `p_6_maf_0.05_v2` let's filter in the same way as above

output
```bash
# minDP10
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 9231 out of a possible 9231 Sites
Run Time = 2.00 seconds

#maxmiss0.8
After filtering, kept 117 out of 117 Individuals
Outputting VCF file...
After filtering, kept 6171 out of a possible 9231 Sites
Run Time = 1.00 seconds

# remove individuals with > 40% missing
Excluding individuals in 'exclude' list
After filtering, kept 109 out of 117 Individuals
Outputting VCF file...
After filtering, kept 6171 out of a possible 6171 Sites
Run Time = 1.00 seconds
```

What individuals were removed? 
```bash
cat remove_individuals.txt

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
The same individuals were removed regardless so I will make a choice based on a better number of loci. 

run interactive mode
```bash
srun -p short -N 1 --pty /bin/bash
module load vcftools


#filter by minimum depth per genotype (minDP = 10)
vcftools --vcf populations.snps.vcf \
         --minDP 10 \
         --recode --recode-INFO-all \
         --out minDP10

#filter for sites present in >= 80% of individuals (refiltered any SNPs found in ≤80% of individuals)
vcftools --vcf minDP10.recode.vcf \
         --max-missing 0.8 \
         --recode --recode-INFO-all \
         --out minDP10_maxmiss0.8

# Now remove individuals with >40% missing data
### 1: compute missingness per individual
vcftools --vcf minDP10_maxmiss0.8.recode.vcf \
         --missing-indv \
         --out missingness

### 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' missingness.imiss > remove_individuals.txt

### now filter out individuals
vcftools --vcf minDP10_maxmiss0.8.recode.vcf \
         --remove remove_individuals.txt \
         --recode --recode-INFO-all \
         --out minDP10_maxmiss0.8_filtInd
```


|Parameters|Loci|Sites kept after filtering|indiv kept after filtering| 
|----------|----|-----|---|
|maf 0.01_v2|33203|109|65631|* this does not add up to the number of sites found in vcftools
|p_1_maf_0.01_v2|35278|109|8601|
|p_1_maf_0.05_v2|35248|109|5474|
|p_2_maf_0.01_v2|30565|109|8611|
|p_3_maf_0.01_v2|27521|109|8625|
|p_4_maf_0.01_v2|24883|109|8652|
|p_5_maf_0.01_v2|20562|109|8652|
|p_6_maf_0.01_v2|15438|109|9147|**
|p_6_maf_0.05_v2|15438|109|6171|
|p_6_maf_0.01_v2 (max missing 0.7)|15438|108|11658|***
|p_6_maf_0.05_v2 (max missing 0.7)|15438|107|7812|


Let's relax the missingness across individuals to 70% (i.e., a site needs to be found in at least 70% of individuals instead of 80% which is pretty high)
```bash
#filter for sites present in >= 70% of individuals (refiltered any SNPs found in ≤80% of individuals)
vcftools --vcf minDP10.recode.vcf \
         --max-missing 0.7 \
         --recode --recode-INFO-all \
         --out minDP10_maxmiss0.7

# Now remove individuals with >40% missing data
### 1: compute missingness per individual
vcftools --vcf minDP10_maxmiss0.7.recode.vcf \
         --missing-indv \
         --out missingness_maxmiss0.7

### 2: generate a list of individuals to remove
awk '$5 > 0.4 {print $1}' missingness_maxmiss0.7.imiss > remove_individuals_maxmiss0.7.txt

### now filter out individuals
vcftools --vcf minDP10_maxmiss0.7.recode.vcf \
         --remove remove_individuals_maxmiss0.7.txt \
         --recode --recode-INFO-all \
         --out minDP10_maxmiss0.7_filtInd

```


What individuals were removed
```bash
(miniconda3) [r.gatins@c0645 p6_maf_0.01_v2]$ cat remove_individuals_0.7.txt
INDV
MA_303_aligned_sorted
MA_304_aligned_sorted
MA_314_aligned_sorted
MA_316_aligned_sorted
MA_318_aligned_sorted
MA_324_aligned_sorted
MD_163_aligned_sorted
ME_253_aligned_sorted
ME_254_aligned_sorted
```
Even if we lose one individual, we gain ~2,000 sites, so we will go with this filtering 


# Final Script

## 1. Stacks
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
  -o $WOR_DIR/stacks/final \
  --popmap $WOR_DIR/popmap/BSB_all \
  --samples $WOR_DIR/samples \
  -X "populations: -r 0.80 --min-maf 0.01 -p 6 --write-single-snp --fstats --vcf --genepop --hwe --structure"

#--------- END Diagnostics/Logging Information---------------
echo = `date` job $JOB_NAME done
echo “using $NSLOTS CPUs”
```


## 2. Filtering with vcftools

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
INPUT_VCF="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/stacks/final/populations.snps.vcf"
OUTDIR="/projects/lotterhos/2020_NOAA_BlackSeaBass_ddRADb/Lotterhos_Project_001/stacks_ref_v2/stacks/final/"

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

### Output

gstacks
```bash
Genotyped 509504 loci:
  effective per-sample coverage: mean=31.0x, stdev=9.5x, min=3.2x, max=48.5x
  mean number of sites per locus: 170.8
  a consistent phasing was found for 731471 of out 819475 (89.3%) diploid loci needing phasing
```

VCF Filtering 
```

```


