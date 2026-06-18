# README - Non-thermal stress Data structure
This section contains all the data files (.csv format) used to analyze the impact of 3 different chromosomal inversions (karyotypes) on the survival rates of Zambia originating *Drosophila melanogaster* under distinct non-thermal stress conditions.

## Experimental design
Two different non-thermal stress were conduct: 
1. Desiccation stress on adult flies (males and females): sexed under gas (CO2) flies were maintain for 18.5 hours in vials containing neither food or water.
2. Starvation stress on adult flies (males and females): flies (some sexed under gas (CO2), some not) were maintain for 41 hours in empty vials. Water was provided by moistening the plug.

## Data overview
The name of the diffeerent files follow the next convention:
- First letters: Developmental Stage (A for Adult)
- Middle letters: Type of Stress (DS = Desiccation Stress, SS= Starvation Stress)
- Last letter: Sex (_F = Females, _M = Males)
  
The 6 files correspond to the following datasets:
1. `ADS_F.csv`: Resitance of Adult Females to Desiccation Stress 
2. `ADS_M.csv`: Resitance of Adult Males to Desiccation Stress
3. `ASS_F.csv`: Resitance of Adult Females to Starvation Stress 
4. `ASS_M.csv`: Resitance of Adult Males to Starvation Stress

### Common Variables 
1. `Karyotype`: The chromosomal inversion group of the tested flies (e.g., *STD*, *3RP*, *2Lt*, *3RK*, and their corresponding heterozygotes like *3RP_HET_1*). *STD* serves as the standard control group.
   
2. `Vial_rep._nr.`: The replicate number for the given karyotype.
   
3. `Name_of_tube`: A unique identifier string combining the karyotype, and replicate number (e.g., 3RP_5).

## File details

### 1. `ADS_F.csv` : 
* **Dimensions:** 277 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of females placed in the vial before the desiccation stress exposure, count
  * `ALIVE_F` — Number of females surviving after the desiccation stress exposure, count
  * `DEAD_F` — Number of females that died during the desiccation stress exposure, count

### 2. `ADS_M.csv` : 
* **Dimensions:** 277 rows × 6 columns
* **Columns:**
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `TOTAL_F` — Initial number of males placed in the vial before the desiccation stress exposure (15 in the majority), count
  * `ALIVE_F` — Number of males surviving after the desiccation stress exposure, count
  * `DEAD_F` — Number of males that died during the desiccation stress exposure, count
  
### 3. `ASS_F.csv` : 
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
 
### 4. `ASS_M.csv` :
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
 
