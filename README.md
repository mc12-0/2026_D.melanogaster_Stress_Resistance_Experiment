# D.melanogaster_Stress_Experiment
# README - Data for IaSC EE meta-analysis, from HT et al.(2026) 

## Title of Dataset
Measures of male reproductive fitness, success in mate choice trials, fecundity induced by successful males, sperm offense; mate harm and productivity induced in a male limited selection experiment. 

## Associated Preprint
- **Title:** Male-benefit adaptation under sex-limited selection shaped by compensatory evolution in Drosophila melanogaster
- **Authors:** 
	- Harshavardhan Thyagarajan, 
	- Mindy G. Baroody, 
	- Imran Sayyed, 
	- Joshua A. Kowal, 
	- Troy Day,  
	- Adam K. Chippindale
- **DOI:** doi.org/10.64898/2026.03.03.709222 

## Dataset Overview
This dataset contains experimental data on fitness, mate choice, mate harm and sperm competition measurements for male-limited and control selection populations. The data pertains to male flies from experimentally evolved *Drosophila melanogaster* populations. The data are provided as CSV files.  

The 4 files correspond to the following datasets:  
1. `CRF.csv` — Competitive reproductive fitness 
2. `MC.csv` — Mate choice and fecundity induction
3. `SC.csv` — Sperm offense 
4. `MH_CG.csv` — Mate harm   

Common variables:
1. Selection 
	ML — Male limited
	MC — Matched control (matched by replicate number, for example ML1 and MC1 are derived from the same parental population)
2. Background - Genetic background of test males
  HC - Home court (ML like genetics), derived from clone generator dams
  WT - Wild type (MC like genetics), derived from control dams
3. Female - Female used for testing male fitness / mating success / sperm offense etc.
  DTP - Clone generator female
  CBr (or) Cp - Recessively marked control females
4. Replicate — This describes the replicate pairs of selected and control populations. It is necessary to clarify here the method of setting up replicates. Commonly, a single source population is used to derive n statistical replicate populations under selection / control conditions. Here rather, we use 3 long separated populations (themselves originally derived as statistical replicates over 20 years prior) to initiate each selection / control pair. 
5. Contest — An identifier for a distinct chamber in which an assay was conducted. 
---

## File Details

### 1. `CRF.csv`

* **Dimensions:** 456 rows × 7 columns
* **Columns:**

  * `Selection` — Described above, categorical
  * `Female` — Described above, categorical
  * `Background` — Described above, categorical
  * `Replicate` — Described above, categorical
  * `Contest` — Described above, categorical
  * `Red` — Number of red-eyed progeny, count
  * `U` — Number of mutant/recessively marked progeny, count
  
### 2. `MC.csv`

* **Dimensions:** 864 rows × 11 columns
* **Columns:**

  * `Selection` — Described above, categorical
  * `Background` — Described above, categorical
  * `Replicate` — Described above, categorical
  * `Female` — Described above, categorical
  * `Contest` — — Described above, categorical
  * `Obs.start.time` — Observation start time, numeric
  * `Mating.start.time` — Time when mating began, numeric
  * `Mating.end.time` — Time when mating ended, numeric
  * (All time variables recorded in the form of minutes since 00h (for example, 11am = 660 mins))
  * `Mated.male` — Identity of mating male (e.g. wt or competitor type), categorical
  * `Offspring_M` — Number of male offspring produced, count
  * `Offspring_F` — Number of female offspring produced, count

### 3. `MH_CG.csv`

* **Dimensions:** 120 rows × 13 columns
* **Columns:**

  * `Selection` — Described above, categorical
  * `Background` — Described above, categorical
  * `Replicate` — Described above, categorical
  * `Contest` — Described above, categorical
  * `Productivity` — Total offspring productivity count, count
  * `Day0`-`Day7` — Surviving females measurement recorded on day 0-day 7, count
  
### 4. `SC.csv`

* **Dimensions:** 831 rows × 8 columns
* **Columns:**

  * `Selection` — Described above, categorical
  * `Female` — Described above, categorical
  * `Background` — Described above, categorical
  * `Replicate` — Described above, categorical
  * `Vial` — Identifier of an offspring collection unit from a contest, categorical
  * `WT` — Number of red-eyed progeny, count
  * `U` — Number of mutant/recessively marked progeny, count
  * `Contest` — Described above, categorical

## Methods

Full details are provided in associated preprint. Across assays, test males were collected at 10–11 days post-oviposition, and experiments were conducted between days 12-14.

1. *Competitive reproductive fitness* Groups of 30 MC or ML flies were combined with 30 recessively marked competitor (Cr) males and 50 Cr virgin females in ventilated “mini-cages” with yeast for 2 days. 

2. *Mate choice trials -> mating success*: Single test males were combined with a Cr male and a virgin Cr female. When one pair initiated copulation, the rival male was removed without disturbance, and the mating event was monitored until separation. Mating latency, mating duration, and the offspring produced (sex ratio, fecundity induction, and eye-color marker) were recorded to identify the successful sire. 

3. *Sperm offense (P2)*: Groups of 12 virgin Cr females were combined with Cr males, and matings were closely observed. After mating, Cr males were removed. On day 13, 10 “target” males (MC or ML) were added to each vial for a second mating, observed for ~3 hours. Ten females per vial were then isolated for 24 hours to lay eggs. Offspring were scored by eye color to determine the proportion sired by the second male. 

4. *Mate harm*: Groups of 12 virgin Cr females were combined with 10 target male, days 12-14, and males were removed after. Female survivor counts assayed for 7 days starting with introduction of males, and productivity was estimated as pupal counts from eggs laid across the 7 day period.
---
