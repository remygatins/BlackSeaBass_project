# Black Sea Bass Genomics Project 
--------------------
-------------------
## Sample selection
First round of DNA extractions (Feb 2020)
(*we are hoping to collect more fish from New England and Maine in Summer 2020 for a second round of genomcis*):

Sample selection **criteria**: total lenght of fish < than 30cm (except those from ME).

| Location | Number of individuals | Notes |
|----------|:-----------------------:|-------|
|Maryland, MD| 20 | 20 selected out of 28, pick smallest; all < 30cm |
|Maine, ME| 5 | all 5 samples available in Feb 2020; 2 < 30cm, 3 > 30cm |
|North Carolina, NC | 15 | all 15 samples available; all < 30cm |
|New Jersey, NJ | 20 | 20 selected out of 30, pick smallest; all < 30cm |
|Southern New England, SN | 4 | selected 4 out of 7 pick smallest; all < 30cm |
| **Total** | **64** | *Goal = 64 samples; multiples of 8 are easier to work with*|


## **PS: this graph is outdated. The updated graph is in 2020-02-BSB-DNA-extractions-summary.md and shows the samples were extracted successfully!!!!!**

![firstroundsamples](https://github.com/thais-neu/BSBproject.md/blob/master/first-round-selected-samples.png)

----------------------------------------------
## Sample preservation method

Samples are in the -80oC freezer, except SN (-20oc freezer).

**Samples in EtOH**: ME, SN (assuming, since liquid is not frozen)

**Samples in RNAlater**: none

**Samples in OGL buffer**: from NC, MD and NJ.

------------------------------------------------
--------------------------------------------------------------------
## DNA extraction from fish fin clips.
----------------------------------------

Detailed protocol here: (https://github.com/thais-neu/BSBproject.md/blob/master/protocols/DNA-extraction-kits.md)

Method used: sample-dependent, noted in each round.

Modifications from the above method: 
   * elution volume varied, see notes on each round of extraction.
   * elution is done in molecular grade water, not buffer from kit.
   
Initial amount of tissue: 15 to 20 mg, unless otherwise noted.

-----------------------------------------------

### Gel prep quick reference 

Detailed protocol here: (https://github.com/thais-neu/BSBproject.md/blob/master/protocols/agarose-gel.md)

Agarose 1.5%
> 1.05 g in 70 mL TBE buffer (small gel)

> 3.15 g in 210 mL TBE buffer (medium gel)

Add 7uL GelRed.

Ladder/sample: ~4uL + 2uL loading dye.

Run at 100V for 60min.

---------------------------------------------------

### 26-Feb-2020

### BSB DNA extractions, round 1.

*Note: upon thawing samples for extraction, I noticed that samples from NC and NJ have very small amount of tissue. Sara recommends not using them in the first round of extraction (practice more before dealing with such small tissues)*

**Round 1 - digestion**

Solutions used: buffer "1", protease "1" (we have two kits in the lab, randomly labelling kits "1" and "2" just to keep track of which ones I'm using for which samples).

Samples selected for 1st round of extraction, ++plate++ method:

| Location | UniqueID | Fin_clip_vial_ID|
|:---------|:--------:|:--------------:|
| MD | Cs_MD_139 | Cs_MD_064 |
| MD | Cs_MD_143 | Cs_MD_068 |
| ME | Cs_ME_165 | 103116 |
| MD | Cs_MD_150 | Cs_MD_076 |
| MD | Cs_MD_140 | Cs_MD_065 |
| MD | Cs_MD_162 | Cs_MD_088 |
| SN | Cs_SN_179 | 79 |
| NEG| NA | NA |

*NEG = negative control, extraction solutions only, no tissue.*

This took ~3h but including the time it takes to thaw out samples in OGL buffer; processing took ~1.5h.

----------------------
### 27-Feb-2020

**Round 1 - purification and gel**

Final elution volume: 400 uL (2x 200uL, combined), visible pellets in the eluate (normal, but maybe too big).

Gel to check, ran at 105V, for 60min.

From left to right: ladder, #139, #143, #165, #150, #140, #162, #179, NEG (all 4uL,except #143 and #150, 6uL).

![gelimageBSBround1](https://github.com/thais-neu/BSBproject.md/blob/master/img/2020-02-27-BSBgel1.jpg)

Most samples didn't work. Could be proteinase 1 is too old, or final DNA was too diluted in 400uL.

Sample #165 is the only one that worked - finclip looks different than the others, has 'cartilage', maybe that is where DNA come from as oppose to the 'soft' tissue of the fin?

**Troubleshooting**

Problem could be that enzyme 1 is too old/not working (expiration date: 28-Mar-2020) and/or final elution volume (400uL) was too high, diluting the extracted DNA too much.

**Round 2 - digestion**

Samples selected for 2nd round of digestion (troubleshooting round), ++plate++ method:

| Location | UniqueID | Fin_clip_vial_ID|enzyme|note|
|:---------|:--------:|:--------------:|-------|----|
| Sara's Iceland cod | NA | #24 | Lotterhos lab "2"| positive ctrl in OGL buffer, at room temp |
| Sara's Iceland cod | NA | #24 | Hughes lab | positive ctrl, in OGL buffer, at room temp |
| MD | Cs_MD_136 | Cs_MD_061 | Lotterhos lab "2"| didnt run in round 1|
| MD | Cs_MD_136 | Cs_MD_061 | Hughes lab | didn't run in round 1|
| SN | Cs_SN_179 | 79 | Lotterhos lab "1" | didn't work in round 1|
| NEG| NA | NA |

Testing the enzymes: samples #24 and #136 are run side by side with Lotterhos enzyme 2 (**not** used in round 1) and enzyme from Hughes lab (newer kit, should work). 
> If enzyme 1 (used in round 1) was bad, these 4 samples will work; if enzyme 2 is also bad, only samples digested with enzyme from Hughes lab will work.

Testing the elution volume: all samples are going to be eluted tomorrow in a final volume of 100uL (2x 50uL).
> If they all work, all enzymes are working, including enzyme 1, and we just need to elute in less water.
> If all but #79 works, enzyme 1 is bad regardless of how much elution water we use.

-------------------------------------------------
### 28-Feb-2020

**Round 2 - purification and gel**

Final elution volume = 100uL (2x 50uL, combined).

From left to right: ladder, +ctrl Lotterhos 2 (L2) enzyme, +ctrl Hughes (H) enzyme, #136 MD sample L2 enzyme, #136 MD sample H enzyme, #179 SN sample L1 enzyme, -ctrl, +ctrl Lotterhos 2 (L2) enzyme (repeat).

![gelimageBSBround2](https://github.com/thais-neu/BSBproject.md/blob/master/img/2020-02-28-BSBgel2.jpg)

**Notes from round 2 gel results:**

* +ctrl worked with both L2 and H enzymes, so both these enzymes are working;
* sample #179 worked with L1 enzyme, eluted in 100uL (not in 400uL from round 1), so L1 enzyme works and **final elution volume needs to be 100uL (as opposed to 400uL)**
* sample #136 (from MD) did not work at all with either enzyme (L2 and H) so the sample is the problem. Why:
   * sample #136 is in OGLfix but so is +ctrl, could it be that BSB DNA is harder to extract from finclips (compared to cod DNA)?
   * finclips from #136 (as well as all others from MD in both rounds) and #179 (and #165 round 1) look different from one another; finclips from MD seem to have no cartilage, just the 'soft' tissue - could that be a problem? 
   * the tissue pellet remaining in lysis plate after extraction is bigger and darker in samples that didn't extract well (all from MD in rounds 1 and 2). Could the DNA be trapped there?
      * extract for longer?
      * use more than 25uL enzyme?
      * higher temperature?

**Image below is from left over tissue pellet in lysis plate** (order matches that of gels)**:**

Solid black pellets correspond to samples that didn't extract well.

![lysisplateimage](https://github.com/thais-neu/BSBproject.md/blob/master/img/2020-02-28-BSBlysisplate.jpg) 

Plan for next week:
* troubleshoot extraction using a MD sample of fish >30cm - one that does **not** meet the selection criteria.
* get advice from Katie on how to tweak digestion.

-----------------------
### 2-March-2020

Summary results from last week: 
* enzymes work;
* positive control (cod sample) and samples from ME and SN worked (eluting in 100uL);
* samples from MD didn't work regardless of enzyme used or elution volume.

*Notes from meeting with Katie: samples from MD that didn't extract are probable the issue; we traced them back and found that those fish were collected in 2016, kept frozen (likely at -20oC), thawed out to collect fin clip in Oct 2019 - this probably led to DNA degradation of the finclips sent to us. We will quantify the DNA in the extracts using Nanodrop and depending on results, contact the sequencing facility to see if it is usable. If so, we extract the remaining MD samples. If not, we drop the MD sample set.*

Next steps:
* Measure DNA concentration of extracted samples using OGL's Nanodrop.
* Move on with extracting the other samples - for the small fin clips (samples from NJ and NC), do a test with a NJ sample that is outside of our selection criteria (fish >30cm) and a positive control (cod) to test if the tubes method of extraction work with these small fin clips.

This is how small a fin clip we have from NJ (and NC); couldn't get a weigh using our scale.

![example_NJfinclip](https://github.com/thais-neu/BSBproject.md/blob/master/img/NJ-finclip-example.jpg)

**Round 3 - digestion**

Samples selected for round 3 of extractions (testing small finclips), ++tube++ method:

| Location | UniqueID | Fin_clip_vial_ID|note|
|:---------|:--------:|:---------------:|----|
| Sara's Iceland cod | NA | #24 |used 5mg|
| NJ | Cs_NJ_107 | Cs_NJ_002 | weigh did not register, fin too small |

Solutions used: Protease "1", Buffer "2" 

-----------------------------------------

### 3-March-2020

**Round 3 - purification**

Final elution volume: 100uL (2x 50uL, combined).

Nanodrop quantification of extracts (rounds 1, 2 and 3; the dotted line is the threshold for 'good quality' DNA (= high DNA to contaminant ratio)).

![Nanodrop_results](https://github.com/thais-neu/BSBproject.md/blob/master/img/2020-03-03-resultssummary.jpg)

* Round 1 = most samples eluted in 400uL are too diluted or extraction didn't work;
* Round 2 = first positive ctrl is surprisingly low concentration, given there's a clear band on the gel;
* Round 2 = concentrations of both #136 is surprisingly high, given no bands on the gel;
* Round 3 = tubes method worked reasonably well for the tiny fin clip (#107); 
* Round 3 = I expected the positive ctrl to be at higher concentration, given that the start material was higher. 

-----------------------------------------------

### 4-March-2020

Checking Nanodrop results against Qubit quantification of samples extracted in round 2 (except NEG) and round 3.

Blue = Nanodrop, Red = Qubit

![qubit_results](https://github.com/thais-neu/BSBproject.md/blob/master/img/2020-03-04-qubitresults.jpg)

--------------------------------
### 6-March-2020

From the sequencing facility: 
> "We request DNA to be sent with volumes >20 ul and concentrations >25 ng/ul as measured by Picogreen or Qubit since Nanodrop will tend to overestimate the amount of double stranded DNA available.  However, we input 10 ul of DNA at 10 ng/ul into our GBS library preparation, so DNA >10 ng/ul as measured by Picogreen or Qubit is acceptable. Our library prep rate includes an initial Picogreen QC step, and we can still attempt library prep on sequencing on samples with concentrations <10 ng/ul. We have successfully sequenced samples with initial concentrations <5 ng/ul, however, the results are more variable and we cannot guarantee that the samples will be successful."

Plan for next week:
* use another positive control (ask Alan for oyster tissue, flash frozen/kept at -80oC; ask Sara for a fin clip that has been kept in EtOH and at -80oC). Positive controls used so far are cod fin clip kept in OGL fix **at room temp** so maybe that's why the extraction yield is low.
* run another 'practice' round with new positive controls - should be getting ~100ng/ul (ask Sara's yield).
* ran a quick test today and it seems like I need at least 30uL elution vol to coat the entire membrane in the extractin column (less volume might not be enough to collect all DNA from the membrane); so when running the samples, elute first round in 30uL and second in 70uL and store the two elutions **separately**. Quantify the first one using Qubit - if >25ng/ul, great; if not, we can concentrate the second elution and add to the first one to increase yield.


--------------------------------

### 11-March-2020

Positive control used in rounds 2 and 3 was a cod tissue sample kept in OGLfix at room temperature. Perhaps room temp wasn't enough to preserve DNA? In round 4, I'm using positive controls kept at -80oC.

Positive controls to work with from now on are 
* (1) cod fin clip tissue preserved in EtOH at -80oC (from Sara), sample ID 17_267 092617 (Box 80_41), concentration in Sara's extract was 91 ng/uL. = **POS_A**
* (2) oyster tissue preserved dry (frozen tissue, no liquid) at -80oC (from Alan), sample ID 221_RNA_Mantle 080116. = **POS_B**

**Round 4 - digestion**

Samples selected for round 4 of extraction (troubleshooting low yield extractions from previous rounds) ++tubes method++: 


| Location | UniqueID | Fin_clip_vial_ID|note|
|:---------|:--------:|:---------------:|----|
| Sara's Gulf of Maine Cod fin | POS_A | #267 |preserved in EtOH, -80oC|
| Alan's Oyster mantle | POS_B | #221 | preserved dry at -80oC | 
| NJ | Cs_NJ_111 | Cs_NJ_006 | used 8mg |

Solutions used: buffer "2", protease "1"

### 12-March-2020

**Round 4 - purification**

Elution volume varied, see Qubit results table below.

Qubit results:

| Location | UniqueID | Fin_clip_vial_ID|note|Qubit conc. (ng/uL) in vol elution 1 (uL) | Qubit conc. (ng/uL) in vol elution 2 (uL) |
|:---------|:--------:|:---------------:|----|:----------------------------------------:|:-----------------------------------------:|
| Sara's Gulf of Maine Cod fin | POS_A | #267 |preserved in EtOH, -80oC| 160 in 50 | NA |
| Alan's Oyster mantle | POS_B | #221 | preserved dry, -80oC | 108 in 50 | NA |
| NJ | Cs_NJ_111 | Cs_NJ_006 | preserved in OGLfix, -80oC, used 8mg | 57.0 in 30 | 4.12 in 70 |

**Notes on Round 4 results:**

Positive controls from samples stored at -80oC yielded much higher concentrations than those kept in OGLfix room temp. Sara's previous extraction of POS_A was 91 ng/uL in 100uL elution and Alan says he typically gets 50-100 in 100uL elution (ethanol precipitation method). So stick with those samples as positive controls when troubleshooting. 

NJ sample, although starting material was lower than what is typically considered ideal, yielded enough DNA to meet sequencing facility requirements.

### 14-April-2020 ###

**Round 5 - digestion**

Extracting 20 samples from NJ.

Solutions used: Buffer "2", Protease "1".

### 15-April-2020 ###

**Round 5 - purification**

Two elutions, 1st elution in 30uL, 2nd elution in 70uL, stored separately at -80oC.

Quantification via Qubit: tried 1uL from first sample, TOO HIGH (outside of detection limit); then all samples were quantified using 1uL of sample  diluted 1:1 with 70oC molecular grade water (dilution factor 2x); all samples quantified are from 1st elution; 2nd elution not quantified.

Round 5 sample info and DNA concentration results:

| Location | UniqueID | Fin_clip vial_ID|extraction tube ID|initial amount (fin clip in mg)|Qubit (ng/uL)|DNA concentration considering dilution factor 2x (ng/uL)| vol. left in 1st elution (uL)| enough for ddRAD? |
|:---------|:--------:|:---------------:|:----------------:|:---------------:|:-----------:|:----------------:|:--------------:|:-----:|
| NJ| Cs_NJ_106 | Cs_NJ_001 | 1  | 2  | 62.4 | 124.8 | 28 | plenty |
| NJ| Cs_NJ_108 | Cs_NJ_003 | 2  | 1  | 51.8 | 103.6 | 29 | plenty |
| NJ| Cs_NJ_109 | Cs_NJ_004 | 3  | 2  | 79.0 | 158.0 | 29 | plenty |
| NJ| Cs_NJ_112 | Cs_NJ_007 | 4  | <1 | 24.4 | 48.8  | 29 | plenty |
| NJ| Cs_NJ_113 | Cs_NJ_008 | 5  | 4  | 64.6 | 129.2 | 29 | plenty |
| NJ| Cs_NJ_114 | Cs_NJ_009 | 6  | <1 | 16.8 | 33.6  | 29 | plenty |
| NJ| Cs_NJ_115 | Cs_NJ_010 | 7  | <1 | 4.42 | 8.84  | 29 | attempt |
| NJ| Cs_NJ_118 | Cs_NJ_013 | 8  | 2  | 60.8 | 121.6 | 29 | plenty |
| NJ| Cs_NJ_119 | Cs_NJ_014 | 9  | 1  | 56.0 | 112.0 | 29 | plenty |
| NJ| Cs_NJ_121 | Cs_NJ_016 | 10 | <1 | 25.0 | 50.0  | 29 | plenty |
| NJ| Cs_NJ_122 | Cs_NJ_017 | 11 | 1  | 42.4 | 84.8  | 29 | plenty |
| NJ| Cs_NJ_124 | Cs_NJ_019 | 12 | <1 | 27.8 | 55.6  | 29 | plenty |
| NJ| Cs_NJ_128 | Cs_NJ_023 | 13 | 2  | 53.8 | 107.6 | 29 | plenty |
| NJ| Cs_NJ_129 | Cs_NJ_024 | 14 | 2  | 68.0 | 136.0 | 29 | plenty |
| NJ| Cs_NJ_130 | Cs_NJ_025 | 15 | <1 | 29.4 | 58.8  | 29 | plenty |
| NJ| Cs_NJ_131 | Cs_NJ_026 | 16 | <1 | 22.0 | 44.0  | 29 | plenty |
| NJ| Cs_NJ_132 | Cs_NJ_027 | 17 | 2  | 44.6 | 89.2  | 29 | plenty |
| NJ| Cs_NJ_133 | Cs_NJ_028 | 18 | <1 | 35.2 | 70.4  | 29 | plenty |
| NJ| Cs_NJ_134 | Cs_NJ_029 | 19 | <1 | 4.26 | 8.52  | 29 | attempt |
| NJ| Cs_NJ_135 | Cs_NJ_030 | 20 | <1 | 3.14 | 6.28  | 29 | attempt |
| NA| NA | NA | blank | NA | TOO LOW | NA | NA | NA |

> <1mg = fins were so small, no weight registered on the scale

> blank = the water used for elution and dilution of samples for Qubit quantification; not a proper NEG control

> standard 1 (low): 60.53; standard 2 (high): 20355.80

### 28-April-2020 ###

**Round 6 - digestion**

Extracting 24 samples, 2 from ME, 7 from SN, and 15 from NC.

Solutions used: Buffer "2", Protease "1".

***Surprise**: when I defrosted the tubes, I see that the samples for NC are muscle tissue, not fin. Consulted with Sara Schaal and she said it is ok to extract with the same protocol. Double-checked with Jon Grabowski to make sure it is from Black Sea Bass, he confirmed. So I moved on with the extractions, exacly as with the fin samples.*

### 29-April-2020 ###

**Round 6 - purification**

Two elutions, 1st elution in 30uL, 2nd elution in 70uL, stored separately at -80oC.

Quantification via Qubit: tried with 1uL of sample first, it worked for all except two samples - #22 (too high) and #33 (too low). I repeated these two, #22 using 1uL sample + 1uL MGwater (dil factor 2x) and #33 using 2uL sample.

Round 6 sample info and DNA concentration results:

| Location | UniqueID | Fin_clip vial_ID|extraction tube ID|initial amount (fin clip in mg)|Qubit (ng/uL)|DNA concentration considering dilution factor 2x (ng/uL)| vol. left in 1st elution (uL)| enough for ddRAD? |
|:---------|:--------:|:---------------:|:----------------:|:---------------:|:-----------:|:----------------:|:--------------:|:-----:|
| ME| Cs_ME_164 | 1027162  | 21 | 0.018 | 120   |NA |29|plenty|
| ME| Cs_ME_165 | 103116   | 22 | 0.019 | 120   |240|28|plenty|
| SN| Cs_SN_179 | 79       | 23 | 0.020 | 108   |NA |29|plenty|
| SN| Cs_SN_182 | 15       | 24 | 0.017 | 88.6  |NA |29|plenty|
| SN| Cs_SN_185 | 87       | 25 | 0.018 | 98.6  |NA |29|plenty|
| SN| Cs_SN_190 | 47       | 26 | 0.016 | 86.6  |NA |29|plenty|
| SN| Cs_SN_191 | 58       | 27 | 0.020 | 88.4  |NA |29|plenty|
| SN| Cs_SN_189 | 4        | 28 | 0.017 | 79.6  |NA |29|plenty|
| SN| Cs_SN_009 | Cs_C_009_DNA | 29 | 0.015 | 100.0 |NA|29|plenty|
| NC| Cs_NC_233 | Cs_NC_31 | 30 | 0.019 | 9.50  |NA |29|maybe|
| NC| Cs_NC_234 | Cs_NC_32 | 31 | 0.018 | 13.3  |NA |29|acceptable|
| NC| Cs_NC_235 | Cs_NC_33 | 32 | 0.018 | 10.1  |NA |29|acceptable|
| NC| Cs_NC_236 | Cs_NC_34 | 33 | 0.019 | 0.102 |NA |27| no |
| NC| Cs_NC_237 | Cs_NC_35 | 34 | 0.018 | 4.82  |NA |29| ?? |
| NC| Cs_NC_238 | Cs_NC_36 | 35 | 0.019 | 4.92  |NA |29| ?? |
| NC| Cs_NC_239 | Cs_NC_37 | 36 | 0.020 | 4.94  |NA |29| ?? |
| NC| Cs_NC_240 | Cs_NC_38 | 37 | 0.019 | 8.50  |NA |29| attempt |
| NC| Cs_NC_241 | Cs_NC_39 | 38 | 0.019 | 3.22  |NA |29| attempt |
| NC| Cs_NC_242 | Cs_NC_40 | 39 | 0.018 | 6.36  |NA |29| attempt |
| NC| Cs_NC_243 | Cs_NC_41 | 40 | 0.018 | 4.54  |NA |29| ?? |
| NC| Cs_NC_244 | Cs_NC_42 | 41 | 0.019 | 0.864 |NA |29| no |
| NC| Cs_NC_245 | Cs_NC_43 | 42 | 0.020 | 5.58  |NA |29| attemtp |
| NC| Cs_NC_246 | Cs_NC_44 | 43 | 0.019 | 4.44  |NA |29| ?? |

> samples from NC are muscle tissue, not fin tissue; sample #38 (exraction tube ID) looked a little slimy, degraded, and had less OGLfix liquid in the tube.

> standard 1 (low): 67.94 ; standard 2 (high): 22828.89

**Notes on Round 6**

Muscle tissue samples did not extract well. There is a lot of material left for a second attempt.

### 19-May-2020 ###

**Round 7 - digestion**

Extracting 20 samples, 3 from ME, 15 from NC, and 2 from MD.

Solutions used: Buffer "2", Protease "1".

> For samples from NC and MD, I used all the leftover muscle (NC) and fin (MD) tissue, given that using 15-20mg of tissue did not yield enough DNA in round 6; volumes for TL buffer and OB protease were adjusted to 1000 uL and 125 uL per tube, respectively. Typically, these volumes are 200 uL and 25 uL, respectively. 

### 20-May-2020 ###

**Round 7 - purification**

> Due to the extra volume of TL and protease needed in the digestion of some of these samples, the samples came out of the thremomixer at a volume of ~1200 uL. After the centrifugation to remove the insoluble pellet, each sample was divided into 3 aliquots of equal (~400uL) volumes. This was needed because equal volumes of BL buffer and 100% ethanol must be added to the sample, which adds up to 1200 uL and the volume of the microcentrifuge tube is 1500 uL. That means that for steps 7-12, corresponding to DNA binding, several rounds of centrifugation were needed to concentrate all 3 aliquots plus the added BL buffer and ethanol into a single spin column per sample. Remaining steps were not adjusted.

Three elutions, 1st elution in 30uL, 2nd elution in 70uL, 3rd elution in 100 uL, stored separately at -80oC.

Quantification via Qubit: tried with 1uL of sample first, it worked for all except five samples - #44, #45, #46, #54 and #59 (too high). I repeated these five using 1uL sample + 1uL MGwater (dil factor 2x). Randomly selected a sample #55 to quantify 2nd and 3rd elutions.

Round 7 sample info and DNA concentration results:

| Location | UniqueID | Fin_clip vial_ID|extraction tube ID|initial amount (fin clip in mg)|Qubit (ng/uL)|DNA concentration considering dilution factor 2x (ng/uL)| vol. left in 1st elution (uL)| enough for ddRAD? |
|:---------|:--------:|:---------------:|:----------------:|:---------------:|:-----------:|:----------------:|:--------------:|:-----:|
| ME| Cs_ME_166 | 102016   | 44 | 0.019 | 98.4  | 196.8 | 28 |plenty|
| ME| Cs_ME_167 | 1027161  | 45 | 0.021 | 106.0 | 212.0 | 29 |plenty|
| ME| Cs_ME_176 | 630      | 46 | 0.018 | 8.12  | 16.24 | 29 | acceptable |
| NC| Cs_NC_233 | Cs_NC_31 | 47 | 0.288 | 108.0 |NA     | 29 |plenty|
| NC| Cs_NC_234 | Cs_NC_32 | 48 | 0.132 | 38.4  |NA     | 29 |plenty|
| NC| Cs_NC_235 | Cs_NC_33 | 49 | 0.194 | 71.2  |NA     | 29 |plenty|
| NC| Cs_NC_236 | Cs_NC_34 | 50 | 0.182 | 0.280 |NA     | 27 |**no**|
| NC| Cs_NC_237 | Cs_NC_35 | 51 | 0.300 | 116.0 |NA     | 29 |plenty|
| NC| Cs_NC_238 | Cs_NC_36 | 52 | 0.136 | 32.0  |NA     | 29 |plenty|
| NC| Cs_NC_239 | Cs_NC_37 | 53 | 0.258 | 108.0 |NA     | 29 |plenty|
| NC| Cs_NC_240 | Cs_NC_38 | 54 | 0.193 | 61.0  |122.0  | 28 |plenty|
| NC| Cs_NC_241 | Cs_NC_39 | 55 | 0.087 | 30.6  |NA     | 29 |plenty|
| NC| Cs_NC_242 | Cs_NC_40 | 56 | 0.176 | 89.6  |NA     | 29 |plenty|
| NC| Cs_NC_243 | Cs_NC_41 | 57 | 0.309 | 106.0 |NA     | 29 |plenty|
| NC| Cs_NC_244 | Cs_NC_42 | 58 | 0.344 | 50.4  |NA     | 29 |plenty|
| NC| Cs_NC_245 | Cs_NC_43 | 59 | 0.097 | 95.6  |191.2  | 28 |plenty|
| NC| Cs_NC_246 | Cs_NC_44 | 60 | 0.140 | 56.2  |NA     | 29 |plenty|
| NC| Cs_NC_247 | Cs_NC_45 | 61 | 0.041 | 7.94  |NA     | 29 |attempt|
| MD| Cs_MD_136 | Cs_MD_61 | 62 | 0.094 | 98.4  |NA     | 29 |plenty|
| MD| Cs_MD_139 | Cs_MD_64 | 63 | 0.177 | 114.0 |NA     | 29 |plenty|

2nd and 3rd elution, extraction tube ID #55
|Elution | Unique ID | Vial ID | Extraction ID | Qubit DNA concentration ug/uL |
|:--------:|:-----------:|:---------:|:---------------:|:------:|
| 2nd | Cs_NC_241 | Cs_NC_39 | 55 | 1.82 | 
| 3rd | Cs_NC_241 | Cs_NC_39 | 55 | 0.170 |

> standard 1 (low): 62.31 ; standard 2 (high): 22262.95

**Notes on Round 7 reagents**

BL buffer from another kit was used for samples 63 (all 3 aliquots) and 62 (2 out of 3 aliquots).
DNA Wash from another kit was used for samples 62 and 63 and wash #1 and all samples in wash #2.


### 25-May-2020 ###

**Round 8 - digestion**

Extracting 19 samples, 18 from MD and 1 from ME (repeat from round 7).

Solutions used: Buffer "2", Protease "1" (extraction ID 64, 65, 66) and Protease "2" (all others).

> For all samples in this round, I used ~0.100-0.150 fin tissue, given that this amount yielded good amounts of DNA in previous rounds form MD samples. ME samples are precious, so I repeated the ME sample that did not in round 7. Volumes for TL buffer and OB protease were adjusted to 1000 uL and 125 uL per tube, respectively. Typically, these volumes are 200 uL and 25 uL, respectively.  

### 26-May-2020 ###

**Round 8 - purification**

> Due to the extra volume of TL and protease needed in the digestion of these samples, the samples came out of the thremomixer at a volume of ~1200 uL. After the centrifugation to remove the insoluble pellet, each sample was divided into 3 aliquots of equal (~400uL) volumes. This was needed because equal volumes of BL buffer and 100% ethanol must be added to the sample, which adds up to 1200 uL and the volume of the microcentrifuge tube is 1500 uL. That means that for steps 7-12, corresponding to DNA binding, several rounds of centrifugation were needed to concentrate all 3 aliquots plus the added BL buffer and ethanol into a single spin column per sample. Remaining steps were not adjusted.

Three elutions, 1st elution in 30uL, 2nd elution in 70uL, 3rd elution in 100 uL, stored separately at -80oC.

Quantification via Qubit: tried with 1uL of sample first, it worked for all except 1 samples - #82 which was diluted 2x, 3x, 4x, all yielding too high. Sample #82 was used to quantify 2nd and 3rd elutions.

Round 8 sample info and DNA concentration results:

| Location | UniqueID | Fin_clip vial_ID|extraction tube ID|initial amount (fin clip in mg)|Qubit (ng/uL)|DNA concentration considering dilution factor 2x (ng/uL)| vol. left in 1st elution (uL)| enough for ddRAD? |
|:---------:|:--------:|:---------------:|:----------------:|:---------------:|:-----------:|:----------------:|:--------------:|:-----:|
| MD| Cs_MD_137 | Cs_MD_62 | 64 | 0.128 | 74.8 |NA| 29 |plenty|
| MD| Cs_MD_138 | Cs_MD_63 | 65 | 0.135 | 96.8 |NA| 29 |plenty|
| MD| Cs_MD_140 | Cs_MD_65 | 66 | 0.114 | 24.0 |NA| 29 |acceptable|
| MD| Cs_MD_141 | Cs_MD_66 | 67 | 0.136 | 82.2 |NA| 29 |plenty|
| MD| Cs_MD_142 | Cs_MD_67 | 68 | 0.136 | 32.0 |NA| 29 |plenty|
| MD| Cs_MD_143 | Cs_MD_68 | 69 | 0.121 | 28.4 |NA| 29 |plenty|
| MD| Cs_MD_145 | Cs_MD_71 | 70 | 0.108 | 30.6 |NA| 29 |plenty|
| MD| Cs_MD_149 | Cs_MD_75 | 71 | 0.155 | 32.0 |NA| 29 |plenty|
| MD| Cs_MD_150 | Cs_MD_76 | 72 | 0.097 | 26.2 |NA| 29 |plenty|
| MD| Cs_MD_151 | Cs_MD_77 | 73 | 0.148 | 33.2 |NA| 29 |plenty|
| MD| Cs_MD_152 | Cs_MD_78 | 74 | 0.146 | 19.9 |NA| 29 |acceptable|
| MD| Cs_MD_154 | Cs_MD_80 | 75 | 0.153 | 38.8 |NA| 29 |plenty|
| MD| Cs_MD_158 | Cs_MD_84 | 76 | 0.153 | 14.8 |NA| 29 |acceptable|
| MD| Cs_MD_159 | Cs_MD_85 | 77 | 0.142 | 19.5 |NA| 29 |acceptable|
| MD| Cs_MD_160 | Cs_MD_86 | 78 | 0.132 | 45.8 |NA| 29 |plent|
| MD| Cs_MD_161 | Cs_MD_87 | 79 | 0.151 | 22.0 |NA| 29 |acceptable|
| MD| Cs_MD_162 | Cs_MD_88 | 80 | 0.108 | 46.4 |NA| 29 |plenty|
| MD| Cs_MD_163 | Cs_MD_89 | 81 | 0.125 | 13.8 |NA| 29 |acceptable|
| ME| Cs_ME_176 | 630      | 82 | 0.129 | 51.0 | 510.0 |  49*  | plenty |

2nd and 3rd elution, extraction tube ID #82
|Elution | Unique ID | Vial ID | Extraction ID | Qubit DNA concentration ug/uL |
|:--------:|:-----------:|:---------:|:---------------:|:------:|
| 2nd | Cs_ME_176 | 630 | 82 | too high | 
| 3rd | Cs_ME_176 | 630 | 82 | 14.5 |

> standard 1 (low): 66.79 ; standard 2 (high): 21501.22

> *this sample had really high concentration, and I needed 5uL of volume to finally get a value within the calibration range, bringing the volume down from 30ul to 25uL. So I added 25uL of Molecular Grade Water to the elution tube and measured one more time, for a leftover volume of 49uL.

**Notes on Round 8 reagents**

HBC buffer from another kit was used for all samples.

## Extraction of 2020-Feb-11th - Mass, Maine, Rhode Island (ADW)

I extracted the 67 samples collected from MA,ME, and RI in 2020 using standard extraction protocols. I used an omega plate DNA Tissue Kit. I eluted samples in 200uL of extraction buffer.

**Break down of samples by location**

| Location | Number of samples |
|:--------:|:-----------------:|
| ME | 15 |
| MA | 30 |
| RI | 22 |
| Total | 67 |

**Extraction Notes**

* Performed overnight digestion using plates and centrifuge in Hugher lab.
* I used silicon plates to cover plate for overnight digestion. I notice some preciptation on the top of plate in morning, which I tried to carefully wipe away as i pull off silicon mat.
* I realized that sample that a peice of ```Cs_MA_150``` fell into well ```H3``` so I decided to just extract the sample twice (wells ```G3``` and ```H3```)
* Samples were eluted into 200uL of elution buffer heated to 57C. 
* I quantified samples and then froze at -20C.
* I requantified samples that were less than <20ng/ul on March 10th 2021. 
* The estimated amount and threshold criteria for >20ng/ul are based on these requants where available.
* Samples in omega elution plate labeled ```Black seabass DNA Extractions Feb 11th 2021 ADW``` in the -20C freezer.

**Sample Information**

| UniqueID | Fin_Clip_Vial_ID | Location | Plate_row | Plate_column | Qubit (ng/uL) | Re-Quant (ng_ul) | Amount (ug, based on 190uL elution) | > 20ng/ul | 
 | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | 
 | Cs_ME_248 | Cs_ME_91 | ME | A | 1 | 20 | NA | 3.8 | GOOD | 
 | Cs_ME_249 | Cs_ME_92 | ME | B | 1 | 29.2 | NA | 5.548 | GOOD | 
 | Cs_ME_250 | Cs_ME_93 | ME | C | 1 | 21.4 | NA | 4.066 | GOOD | 
 | Cs_ME_251 | Cs_ME_94 | ME | D | 1 | 42.4 | NA | 8.056 | GOOD | 
 | Cs_ME_252 | Cs_ME_95 | ME | E | 1 | 39.6 | NA | 7.524 | GOOD | 
 | Cs_ME_253 | Cs_ME_96 | ME | F | 1 | 14.4 | 22.2 | 4.218 | GOOD | 
 | Cs_ME_254 | Cs_ME_97 | ME | G | 1 | 5.7 | 6.08 | 1.1552 | RE-EXTRACT | 
 | Cs_ME_255 | Cs_ME_98 | ME | H | 1 | 25.2 | NA | 4.788 | GOOD | 
 | Cs_ME_256 | Cs_ME_99 | ME | A | 2 | 37.2 | NA | 7.068 | GOOD | 
 | Cs_ME_257 | Cs_ME_100 | ME | B | 2 | 40.6 | NA | 7.714 | GOOD | 
 | Cs_ME_258 | Cs_ME_101 | ME | C | 2 | 0.3 | 69.04 | 13.1176 | GOOD | 
 | Cs_ME_259 | Cs_ME_102 | ME | D | 2 | 17.9 | 10.4 | 1.976 | RE-EXTRACT | 
 | Cs_ME_260 | Cs_ME_103 | ME | E | 2 | 15.2 | 14.4 | 2.736 | RE-EXTRACT | 
 | Cs_ME_261 | Cs_ME_104 | ME | F | 2 | 36.6 | NA | 6.954 | GOOD | 
 | Cs_ME_262 | Cs_ME_105 | ME | G | 2 | 46.6 | NA | 8.854 | GOOD | 
 | Cs_MA_298 | Cs_MA_141 | MA | H | 2 | 14.6 | 19.3 | 3.667 | GOOD | 
 | Cs_MA_299 | Cs_MA_142 | MA | A | 3 | 9.18 | 27.4 | 5.206 | GOOD | 
 | Cs_MA_300 | Cs_MA_143 | MA | B | 3 | 28.6 | NA | 5.434 | GOOD | 
 | Cs_MA_302 | Cs_MA_145 | MA | C | 3 | 27.2 | NA | 5.168 | GOOD | 
 | Cs_MA_303 | Cs_MA_146 | MA | D | 3 | 27.4 | NA | 5.206 | GOOD | 
 | Cs_MA_304 | Cs_MA_147 | MA | E | 3 | 40.8 | NA | 7.752 | GOOD | 
 | Cs_MA_306 | Cs_MA_149 | MA | F | 3 | 29.3 | NA | 5.567 | GOOD | 
 | Cs_MA_307 | Cs_MA_150 | MA | G | 3 | 32.6 | NA | 6.194 | GOOD | 
 | Cs_MA_307 | Cs_MA_150 | MA | H | 3 | 17.4 | NA | 3.306 | GOOD | 
 | Cs_MA_309 | Cs_MA_152 | MA | A | 4 | 22.6 | NA | 4.294 | GOOD | 
 | Cs_MA_310 | Cs_MA_153 | MA | B | 4 | 38 | NA | 7.22 | GOOD | 
 | Cs_MA_311 | Cs_MA_154 | MA | C | 4 | 19.9 | 24 | 4.56 | GOOD | 
 | Cs_MA_313 | Cs_MA_156 | MA | D | 4 | 48.8 | NA | 9.272 | GOOD | 
 | Cs_MA_314 | Cs_MA_157 | MA | E | 4 | 39 | NA | 7.41 | GOOD | 
 | Cs_MA_315 | Cs_MA_158 | MA | F | 4 | 40.4 | NA | 7.676 | GOOD | 
 | Cs_MA_316 | Cs_MA_159 | MA | G | 4 | 39.4 | NA | 7.486 | GOOD | 
 | Cs_MA_317 | Cs_MA_160 | MA | H | 4 | 3.64 | 4.22 | 0.8018 | RE-EXTRACT | 
 | Cs_MA_318 | Cs_MA_161 | MA | A | 5 | 36.3 | NA | 6.897 | GOOD | 
 | Cs_MA_319 | Cs_MA_162 | MA | B | 5 | 18.2 | 23 | 4.37 | GOOD | 
 | Cs_MA_320 | Cs_MA_163 | MA | C | 5 | 41.2 | NA | 7.828 | GOOD | 
 | Cs_MA_321 | Cs_MA_164 | MA | D | 5 | 57.8 | NA | 10.982 | GOOD | 
 | Cs_MA_323 | Cs_MA_166 | MA | E | 5 | 52.6 | NA | 9.994 | GOOD | 
 | Cs_MA_324 | Cs_MA_167 | MA | F | 5 | 38.6 | NA | 7.334 | GOOD | 
 | Cs_MA_325 | Cs_MA_168 | MA | G | 5 | 42 | NA | 7.98 | GOOD | 
 | Cs_MA_327 | Cs_MA_170 | MA | H | 5 | 37.6 | NA | 7.144 | GOOD | 
 | Cs_RI_328 | Cs_RI_171 | RI | A | 6 | 50 | NA | 9.5 | GOOD | 
 | Cs_RI_329 | Cs_RI_172 | RI | B | 6 | 45.4 | NA | 8.626 | GOOD | 
 | Cs_RI_330 | Cs_RI_173 | RI | C | 6 | 61.4 | NA | 11.666 | GOOD | 
 | Cs_RI_331 | Cs_RI_174 | RI | D | 6 | 12.4 | 12 | 2.28 | RE-EXTRACT | 
 | Cs_RI_332 | Cs_RI_175 | RI | E | 6 | 91.2 | NA | 17.328 | GOOD | 
 | Cs_RI_333 | Cs_RI_176 | RI | F | 6 | 114 | NA | 21.66 | GOOD | 
 | Cs_RI_334 | Cs_RI_177 | RI | G | 6 | 47 | NA | 8.93 | GOOD | 
 | Cs_RI_335 | Cs_RI_178 | RI | H | 6 | 17.3 | 18.4 | 3.496 | RE-EXTRACT | 
 | Cs_RI_336 | Cs_RI_179 | RI | A | 7 | 15.8 | 25.4 | 4.826 | GOOD | 
 | Cs_RI_337 | Cs_RI_180 | RI | B | 7 | 88.2 | NA | 16.758 | GOOD | 
 | Cs_RI_338 | Cs_RI_181 | RI | C | 7 | 83.6 | NA | 15.884 | GOOD | 
 | Cs_RI_339 | Cs_RI_182 | RI | D | 7 | 98.2 | NA | 18.658 | GOOD | 
 | Cs_RI_340 | Cs_RI_183 | RI | E | 7 | 90.4 | NA | 17.176 | GOOD | 
 | Cs_RI_341 | Cs_RI_184 | RI | F | 7 | 98.8 | NA | 18.772 | GOOD | 
 | Cs_RI_342 | Cs_RI_185 | RI | G | 7 | 88.2 | NA | 16.758 | GOOD | 
 | Cs_RI_343 | Cs_RI_186 | RI | H | 7 | 19.3 | 31.8 | 6.042 | GOOD | 
 | Cs_RI_344 | Cs_RI_187 | RI | A | 8 | 62.4 | NA | 11.856 | GOOD | 
 | Cs_RI_345 | Cs_RI_188 | RI | B | 8 | 15.1 | 17.2 | 3.268 | RE-EXTRACT | 
 | Cs_RI_346 | Cs_RI_189 | RI | C | 8 | 78 | NA | 14.82 | GOOD | 
 | Cs_RI_347 | Cs_RI_190 | RI | D | 8 | 61.4 | NA | 11.666 | GOOD | 
 | Cs_RI_348 | Cs_RI_191 | RI | E | 8 | 3.1 | 3.1 | 0.589 | RE-EXTRACT | 
 | Cs_RI_349 | Cs_RI_192 | RI | F | 8 | 21 | NA | 3.99 | GOOD | 

## Re-extraction of previous samples - 2020-Mar-12th - Mass, Maine, Rhode Island (ADW)

I re-extracted samples that failed during the Feb-2021 extraction. Standard protocols from the omega kit were used.

I eluted into 100uL of molecular grade water. Samples were stored in the 9th row of the Feb-2021 extraction plate.

| UniqueID | Fin_Clip_Vial_ID | Location | Plate_row | Plate_column | Qubit (ng/uL) | Re-Quant (ng_ul) | Amount (ng, based on 100uL elution) | > 20ng/ul | 
| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | 
| Cs_ME_254 | Cs_ME_97 | ME | A | 9 | 46.2 | NA | 4620 | GOOD | 
| Cs_ME_259 | Cs_ME_102 | ME | B | 9 | 1.76 | NA | 176 | RE-EXTRACT | 
| Cs_ME_260 | Cs_ME_103 | ME | C | 9 | 15.4 | NA | 1540 | RE-EXTRACT | 
| Cs_MA_317 | Cs_MA_160 | MA | D | 9 | 6.88 | NA | 688 | RE-EXTRACT | 
| Cs_RI_331 | Cs_RI_174 | RI | E | 9 | 80.6 | NA | 8060 | GOOD | 
| Cs_RI_335 | Cs_RI_178 | RI | F | 9 | 90.0 | NA | 9000 | GOOD | 
| Cs_RI_345 | Cs_RI_188 | RI | G | 9 | TOO HIGH | NA | NA | GOOD | 
| Cs_RI_348 | Cs_RI_191 | RI | H | 9 | 25.6 | NA | 2560 | GOOD | 

There were still three samples that could be re-extracted. I believe there is enough sample to do this. We may want to look into either borrowing a couple of individual extraction preps from another lab or buying a new kit if we want to extract these samples.

Alternatively, there should be enough material among the two extracts if we want to try and concentrate it with a bead wash.


## 12 April 2021 

## List of samples submitted to University of Minessota Genomics Center for ddRAD library prep and sequencing. 

The samples below are distributed in two 96-well plates ID 2021_Centropristis_striata_BlaSeaBass_ddRADSeq_Plate1 and Plate2.

28 ul of each sample was sent to library prep and sequencing; concentration of each sample (based on Qubit) and plate layout can be found in 

![UMN_samplesheet](https://github.com/thais-neu/BSBproject.md/blob/master/img/20210412_UMNsamplesheet)

| Location | Number of samples |
|:--------:|:-----------------:|
|       MA |        23         |
|       MD |        20         |
|       ME |        20         |
|       NC |        13         |
|       NJ |        17         |
|       RI |         7         |
|       SN |        22         |
|    Total |       122         |

