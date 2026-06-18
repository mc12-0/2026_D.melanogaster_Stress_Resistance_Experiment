# README - Thermal stress Data structure
This section contains all the data files (.csv format) used to analyze the impact of 3 different chromosomal inversions (karyotypes) on the survival rates of Zambia originating *Drosophila melanogaster* across various developmental stages and under distinct thermal stress conditions.

## Experimental design
Different thermal stress were conduct: 
1. Cold Shock on adult flies (males and females): Vials containing generally 15 males and 15 female were placed in ice at 0°C for 18 hours.
2. Heat Shock on adult flies (males and females): Vials containing generally 15 males and 15 female were placed in an oven at 39°C for 75 minutes.
3. Cold Shock on larvae: Vials containg 60 hours-old larvae were placed in ice at 0°C for 12 hours.
4. Heat Shock on larvae: Vials containg 60 hours-old larvae were placed in ice at 39°C for 5.5 hours.

## Data overview
The files in this repository follow a strict naming convention to quickly identify the specific experiment:
- First letter(s): Developmental Stage (A = Adult, L = Larvae)
- Middle letters: Type of Stress (CS = Cold Shock, HS = Heat Shock)
- Last letter (if applicable): Sex (_F = Females, _M = Males)
  
The 6 files correspond to the following datasets:
1. `ACS_F.csv`: Resitance of Adult Females to Cold Shock  
2. `ACS_M.csv`: Resitance of Adult Males to Cold Shock
3. `AHS_F.csv`: Resitance of Adult Females to Heat Shock 
4. `AHS_M.csv`: Resitance of Adult Males to Heat Shock 
5. `LCS.csv`: Rsistance of Larvae to Cold Shock
6. `LHS.csv`: Resistance of Larvae to Heat Shock

### Common Variables 
1. `Karyotype`: The chromosomal inversion group of the tested flies (e.g., *STD*, *3RP*, *2Lt*, *3RK*, and their corresponding heterozygotes like *3RP_HET_1*). *STD* serves as the standard control group.
   
2. `Vial_rep._nr.`: The replicate number for the given karyotype.
   
3. `Name_of_tube`: A unique identifier string combining the experiment type, karyotype, and replicate number (e.g., ACS 3RP_HET_1_5).

## File details

### 1. `ACS_F.csv` : Adult Cold Shock - Females
* **Dimensions:** 193 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the cold exposure (15 in the majority), count
  * `ALIVE_F` — Number of females surviving after the cold exposure and the 4 hour recovery, count
  * `DEAD_F` — Number of females that died during the cold exposure, count

### 2. `ACS_M.csv` : Adult Cold Shock - Males
* **Dimensions:** 193 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the cold exposure (15 in the majority), count
  * `ALIVE_F` — Number of males surviving after the cold exposure and the 4 hour recovery, count
  * `DEAD_F` — Number of males that died during the cold exposure, count
  
### 3. `AHS_F.csv` : Adult Heat Shock - Females
* **Dimensions:** 197 rows × 7 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the heat shock (15 in the majority), count
  * `ALIVE_F` — Number of females surviving after the heat shock and the 4 hour recovery, count
  * `DEAD_F` — Number of females that died during the heat shock, count
  * `BORDERS` —Indicate if the vials experienced edge effects during the heat shock (e.g., "YES" or "NO"). Vials marked "YES" were excluded from the final statistical analysis, categorial
 
### 4. `AHS_M.csv` : Adult Heat Shock - Males
* **Dimensions:** 197 rows × 7 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the heat shock (15 in the majority), count
  * `ALIVE_F` — Number of males surviving after the heat shock and the 4 hour recovery, count
  * `DEAD_F` — Number of males that died during the heat shock, count
  * `BORDERS` —Indicate if the vials experienced edge effects during the heat shock (e.g., "YES" or "NO"), categorial

### 5. `LCS.csv` : Larval Cold Shock
* **Dimensions:** 219 rows × 7 columns
* **Columns:**
  * `No.` —  Index number for the row
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Total_eggs` — Number of eggs initially placed in the vial (typically 30), count
  * `Adults` — Number of individuals that successfully developed into adult flies after the larval cold exposure, count
  * `Pupae` — Number of individuals that reached the pupal stage but do not emerge as adults, count
 
 
### 6. `LHS.csv` : Larval Heat Shock
* **Dimensions:** 201 rows × 7 columns
* **Columns:**
  * `No.` —  Index number for the row
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Total_eggs` — Number of eggs initially placed in the vial (typically 30), count
  * `Adults` — Number of individuals that successfully developed into adult flies after the larval cold exposure, count
  * `Pupae` — Number of individuals that reached the pupal stage but do not emerge as adults, count
