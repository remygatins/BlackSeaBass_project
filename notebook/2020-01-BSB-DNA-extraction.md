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

Positive controls from samples stored at -80oc yielded much higher concentrations than those kept in OGLfix room temp. Sara's previous extraction of POS_A was 91 ng/uL in 100uL elution and Alan says he typically gets 50-100 in 100uL elution (ethanol precipitation method), so stick with those samples as positive controls when troubleshooting. 

NJ sample, although starting material was lower than what is typically considered ideal, yielded enough DNA to meet sequencing facility requirements.
