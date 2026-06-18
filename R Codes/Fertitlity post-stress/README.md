# README - Fertility post-stress Data structure
This section contains all the data files (.csv format) used to analyze the impact of 3 different chromosomal inversions (karyotypes) on the female fertility of Zambia originating *Drosophila melanogaster* after distinct non-thermal stress conditions.

## Experimental design
Female fertility was measured after two non-thermal stresses: 
1. Desiccation stress on adult flies (males and females): sexed under gas (CO2) flies were maintain for 18.5 hours in vials containing neither food or water.
2. Starvation stress on adult flies (males and females): flies (some sexed under gas (CO2), some not) were maintain for 41 hours in empty vials. Water was provided by moistening the plug.

The females that survived the non-thermal stresses were placed individually in vials containing food. After isolating the females, we kept them at 25°C for 3 days to allow them to lay eggs. After 3 days, the females were discarded. 
We allowed their offspring to develop, and after 14 days (for post-desiccation fertility) and 15 days (for post-starvation fertility), we froze the offspring and counted them.  

## Data overview
The name of the diffeerent files follow the next convention:
- First letter: Fertility 
- Last letters: Type of Stress on the adults (A = Adult, DS = Desiccation Stress, SS= Starvation Stress)

The 2 files correspond to the following datasets:
1. `FADS.csv`: Fertility after Adult Desiccation Stresss 
2. `FASS.csv`: Fertility after Adult Starvation Stress

### Common Variables 
1. `No.`: Index number for the row
   
2. `Karyotype`: The chromosomal inversion group of the tested flies (e.g., *STD*, *3RP*, *2Lt*, *3RK*, and their corresponding heterozygotes like *3RP_HET_1*). *STD* serves as the standard control group.
   
3. `Vial_rep._nr.`: The replicate number for the given karyotype.
  
4. `Name_of_tube`: A unique identifier string combining the karyotype, and replicate number (e.g., 3RP_5).

## File details

### 1. `FADS.csv` : 
* **Dimensions:** 342 rows × 6 columns
* **Columns:**
  * `No.`— Described above
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Adults` — Number of offspring that successfully emerged as adult flies following the parents' exposure to the desiccation stress, count
  * `Pupae_(dark/dead)` — Number of offspring that reached the pupal stage following the parents' exposure to the desiccation stress but died or failed to emerge as adults, count
  
### 2. `FASS.csv` : 
* **Dimensions:** 380 rows × 7 columns
* **Columns:**
  * `No.`— Described above
  * `Karyotype` — Described above, categorical
  * `Vial_rep._nr.` — Described above, categorical
  * `Name_of_tube` — Described above, categorical
  * `Adults` — Number of offspring that successfully emerged as adult flies following the parents' exposure to the desiccation stress, count
  * `Pupae` — Number of offspring that reached the pupal stage following the parents' exposure to the desiccation stress but died or failed to emerge as adults, count
  * `GAS` — Indicate if the flies in the vials were sexed under gas or not prior the starvation stress (e.g., "YES" or "NO"), categorial
 
