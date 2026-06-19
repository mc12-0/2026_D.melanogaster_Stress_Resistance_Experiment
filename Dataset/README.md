# README - Data Structure
This section contains all the data files (.csv format) used to analyze the impact of 3 different chromosomal inversions (karyotypes) on stress resistance  of Zambia originating *Drosophila melanogaster*. Multiple stress experimental setup were performed: 
1. Thermal stress (heat/cold shock) on both larvae and adult flies.
2. Non-thermal stress (desiccation/starvation) on adult flies.
3. Measurement of female fertility after a non-thermal stress (desiccation/starvation).

## Table of contents
- [Files Name](#files-name)
- [Common variables across all files](#common-variable-across-all-files)
- [Thermal Stress](#thermal-stress)
- [Non-thermal Stress](#non-thermal-stress)
- [Female Fertility post-stress](#female-fertility-post-stress)

## Files Name
The name of the different files follow the next convention:
- First letter(s): Developmental Stage (`A` = Adult, `L` = Larvae) or Fertility (`F`)
- Middle letters: Type of Stress (`CS` = Cold Shock, `HS` = Heat Shock, `DS` = Desiccation Stress, `SS`= Starvation Stress )
- Last letter (if applicable): Sex (`_F` = Females, `_M` = Males)

The 12 files correspond to the follwing dataset:
* **Thermal Stress:**
    1. `ACS_F.csv`: Resitance of Adult Females to Cold Shock  
    2. `ACS_M.csv`: Resitance of Adult Males to Cold Shock
    3. `AHS_F.csv`: Resitance of Adult Females to Heat Shock 
    4. `AHS_M.csv`: Resitance of Adult Males to Heat Shock 
    5. `LCS.csv`: Rsistance of Larvae to Cold Shock
    6. `LHS.csv`: Resistance of Larvae to Heat Shock
* **Non-thermal Stress:**
    1. `ADS_F.csv`: Resitance of Adult Females to Desiccation Stress
    2. `ADS_M.csv`: Resitance of Adult Males to Desiccation Stress
    3. `ASS_F.csv`: Resitance of Adult Females to Starvation Stress
    4. `ASS_M.csv`: Resitance of Adult Males to Starvation Stress
* **Female fertility post-stress:**
    1. `FADS.csv`: Fertility after Adult Desiccation Stress
    2. `FASS.csv`: Fertility after Adult Starvation Stress
 
## Common variables across all files
1. `Karyotype`: The chromosomal inversion group of the tested flies (e.g., *STD*, *3RP*, *2Lt*, *3RK*, and their corresponding heterozygotes like *3RP_HET_1*). *STD* serves as the standard control group.
   
2. `Vial_rep._nr.`: The replicate number for the given karyotype.
   
3. `Name_of_tube`: A unique identifier string combining the karyotype, and replicate number (e.g., 3RP_5).

## Thermal Stress
Different thermal stress were conduct: 
1. **Cold Shock on adult flies** (males and females): Vials containing generally 15 males and 15 female were placed in ice at 0°C for 18 hours.
2. **Heat Shock on adult flies** (males and females): Vials containing generally 15 males and 15 female were placed in an oven at 39°C for 75 minutes.
3. **Cold Shock on larvae :** Vials containg 60 hours-old larvae were placed in ice at 0°C for 12 hours.
4. **Heat Shock on larvae:** Vials containg 60 hours-old larvae were placed in ice at 39°C for 5.5 hours.

### File details

#### 1. `ACS_F.csv`:
* **Dimensions:** 193 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the cold exposure (15 in the majority), count
  * `ALIVE_F` — Number of females surviving after the cold exposure and the 4 hour recovery, count
  * `DEAD_F` — Number of females that died during the cold exposure, count

#### 2. `ACS_M.csv`:
* **Dimensions:** 193 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the cold exposure (15 in the majority), count
  * `ALIVE_F` — Number of males surviving after the cold exposure and the 4 hour recovery, count
  * `DEAD_F` — Number of males that died during the cold exposure, count
  
#### 3. `AHS_F.csv`:
* **Dimensions:** 197 rows × 7 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the heat shock (15 in the majority), count
  * `ALIVE_F` — Number of females surviving after the heat shock and the 4 hour recovery, count
  * `DEAD_F` — Number of females that died during the heat shock, count
  * `BORDERS` —Indicate if the vials experienced edge effects during the heat shock (e.g., "YES" or "NO"). Vials marked "YES" were excluded from the final statistical analysis, categorial
 
#### 4. `AHS_M.csv`:
* **Dimensions:** 197 rows × 7 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the heat shock (15 in the majority), count
  * `ALIVE_F` — Number of males surviving after the heat shock and the 4 hour recovery, count
  * `DEAD_F` — Number of males that died during the heat shock, count
  * `BORDERS` —Indicate if the vials experienced edge effects during the heat shock (e.g., "YES" or "NO"), categorial

#### 5. `LCS.csv`:
* **Dimensions:** 219 rows × 7 columns
* **Columns:**
  * `No.` —  Index number for the row
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Total_eggs` — Number of eggs initially placed in the vial (typically 30), count
  * `Adults` — Number of individuals that successfully developed into adult flies after the larval cold exposure, count
  * `Pupae` — Number of individuals that reached the pupal stage but do not emerge as adults, count
 
#### 6. `LHS.csv`:
* **Dimensions:** 201 rows × 7 columns
* **Columns:**
  * `No.` —  Index number for the row
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Total_eggs` — Number of eggs initially placed in the vial (typically 30), count
  * `Adults` — Number of individuals that successfully developed into adult flies after the larval cold exposure, count
  * `Pupae` — Number of individuals that reached the pupal stage but do not emerge as adults, count

## Non-thermal Stress
Two different non-thermal stress were conduct: 
1. **Desiccation stress** on adult flies (males and females): sexed under gas (CO2) flies were maintain for 18.5 hours in vials containing neither food or water.
2. **Starvation stress** on adult flies (males and females): flies (some sexed under gas (CO2), some not) were maintain for 41 hours in empty vials. Water was provided by moistening the plug.

### File details

#### 1. `ADS_F.csv` : 
* **Dimensions:** 277 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the desiccation stress exposure, count
  * `ALIVE_F` — Number of females surviving after the desiccation stress exposure, count
  * `DEAD_F` — Number of females that died during the desiccation stress exposure, count

#### 2. `ADS_M.csv` : 
* **Dimensions:** 277 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the desiccation stress exposure (15 in the majority), count
  * `ALIVE_F` — Number of males surviving after the desiccation stress exposure, count
  * `DEAD_F` — Number of males that died during the desiccation stress exposure, count
  
#### 3. `ASS_F.csv` : 
* **Dimensions:** 242 rows × 9 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the starvation stress exposure, count
  * `ALIVE_F` — Number of females surviving after the stravation stress exposure, count
  * `DEAD_F` — Number of females that died during the stravation stress exposure, count
  * `GAS` — Indicate if the flies in the vials were sexed under gas or not prior the starvation stress (e.g., "YES" or "NO"). Only vials marked "YES" were keep for the final statistical analysis, categorial
  * `BORDERS` — Indicate if the vials experienced edge effects during the heat shock (e.g., "YES" or "NO"). Vials marked "YES" were excluded from the final statistical analysis, categorial
 
#### 4. `ASS_M.csv` :
* **Dimensions:** 242 rows × 9 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the starvation stress exposure, count
  * `ALIVE_F` — Number of males surviving after the stravation stress exposure, count
  * `DEAD_F` — Number of males that died during the stravation stress exposure, count
  * `GAS` — Indicate if the flies in the vials were sexed under gas or not prior the starvation stress (e.g., "YES" or "NO"). Only vials marked "YES" were keep for the final statistical analysis, categorial
  * `BORDERS` — Indicate if the vials experienced edge effects during the heat shock (e.g., "YES" or "NO"). Vials marked "YES" were excluded from the final statistical analysis, categorial


## Female Fertility post-stress
Female fertility was measured after two non-thermal stresses: 
1. Desiccation stress on adult flies (males and females): sexed under gas (CO2) flies were maintain for 18.5 hours in vials containing neither food or water.
2. Starvation stress on adult flies (males and females): flies (some sexed under gas (CO2), some not) were maintain for 41 hours in empty vials. Water was provided by moistening the plug.

The females that survived the non-thermal stresses were placed individually in vials containing food. After isolating the females, we kept them at 25°C for 3 days to allow them to lay eggs. After 3 days, the females were discarded. 

We allowed their offspring to develop, and after 14 days (for post-desiccation fertility) and 15 days (for post-starvation fertility), we froze the offspring and counted them.  


### File details

#### 1. `FADS.csv` : 
* **Dimensions:** 342 rows × 6 columns
* **Columns:**
  * `No.`— Index number for the row
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Adults` — Number of offspring that successfully emerged as adult flies following the parents' exposure to the desiccation stress, count
  * `Pupae_(dark/dead)` — Number of offspring that reached the pupal stage following the parents' exposure to the desiccation stress but died or failed to emerge as adults, count
  
#### 2. `FASS.csv` : 
* **Dimensions:** 380 rows × 7 columns
* **Columns:**
  * `No.`— Index number for the row
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Adults` — Number of offspring that successfully emerged as adult flies following the parents' exposure to the desiccation stress, count
  * `Pupae` — Number of offspring that reached the pupal stage following the parents' exposure to the desiccation stress but died or failed to emerge as adults, count
  * `GAS` — Indicate if the flies in the vials were sexed under gas or not prior the starvation stress (e.g., "YES" or "NO"), categorical
 

