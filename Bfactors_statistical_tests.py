
##SciPy stats — t-tests, Mann-Whitney, KDE, distributions
# Test if B-factors differ between high and low resolution structures
# B-factors from RCSB — fetch 200 structures


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
import json
import seaborn as sns
from scipy import stats
from Bio.PDB import PDBParser, PDBList
import warnings
warnings.filterwarnings('ignore')

##First fetch 200 structures of certain protein from RCSB, using the same query based approach used in "RCSB_metadata_analyse.py" in Week 01###

def fetch_pdb_structure_ids_resolution(res_min, res_max, rows=100):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query={
        "query":{
            "type":"terminal",
            "service": "text",
            "parameters": {
                "attribute":"rcsb_entry_info.resolution_combined",
                "operator": "range",
                "value" : {
                    "from" : res_min,
                    "to" : res_max,
                    "include_lower": True,
                    "include_upper" : True

                }
            }   
        },

        "return_type" :"entry",
        "request_options" : {
            "paginate" :{"start":0, "rows" :rows},
            "sort" : [{"sort_by" : "score", "direction" : "desc"}]
        }
    }

    response = requests.post(url, json=query)
    results = response.json()
    # with open("res_ids.json", "w") as f:
    #     json.dump(results, f, indent= 4) (Just to debug to understand the resultant file for IDs )
    #return results
    ids = [r["identifier"] for r in results.get('result_set', [])]
    print(f"{len(ids)} found within the Resolution {res_min} - {res_max} Ang")
    return ids

high_res_ids = fetch_pdb_structure_ids_resolution(0.5,1.5,rows=100)
low_res_ids = fetch_pdb_structure_ids_resolution(3.0, 4.0, rows = 100)
# with open("high_res_ids.txt", "w") as f:
#     json.dump(high_res_ids, f, indent= 4)

# with open("low_res_ids.txt", "w") as f:
#     json.dump(low_res_ids,  f, indent= 4)

###Fetch the pdb structures and their bfactors for each pdb file###

###For this we can create a function which downloads the pdb structure, get the C-alpha atom and it's bfactor"""##

def get_b_factor(pdb_id):
    try:
        pdb_test = PDBList()
        pdb_test.retrieve_pdb_file(pdb_id, file_format = "pdb", pdir= ".")
            
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, file= f"pdb{pdb_id.lower()}.ent")

        bfactor = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if 'CA' in residue:
                        bfactor.append(residue['CA'].get_bfactor())
            break ####considers only the first model.
        if len(bfactor)>10:
            #print(f"Length of bfactor is: { len(bfactor)}, hence computing the mean")
            return np.mean(bfactor)
        else:
            return None
    except Exception as e:
        return None
####The try and except command of function will help us to make sure that even if certain PDBs are not working properly, for instance they couldn't be downloaded are not available in the requested format###
##Then, that can be ignored and moved to next PDBID##

###Collect bfactors for high_res and low_res IDs#####

high_res_bfactors = []

for id in high_res_ids[:50]:
    bf = get_b_factor(id)
    if bf is not None:
        high_res_bfactors.append(bf)
    
# with open("high_res_bfactors.txt", "w") as f:
#     json.dump(high_res_bfactors, f, indent= 4)
    

low_res_bfactors = []
for id in low_res_ids[:50]:
    bf = get_b_factor(id)
    if bf is not None:
        low_res_bfactors.append(bf)

# with open("low_res_bfactors.txt", "w") as f:
#     json.dump(low_res_bfactors, f, indent= 4)

high_res_bfactors = np.array(high_res_bfactors)
low_res_bfactors = np.array(low_res_bfactors) ###This will yield a 10 elements in a single row

print(f"\nHigh-res group: {len(high_res_bfactors)} structures")
print(f"Low-res group:  {len(low_res_bfactors)} structures")


######Statistics#####
###Here we will use Descriptive Statistics, Statistical tests which include t-tests,Mann-Whitney U tests, Shapiro Wilk tests###
###To know more about these tests check the readme file section statistical tests###


###Descriptive statistics###
"""Descriptive Statistics of High Res B-factors"""
print("Descriptive Statistics of High Res B-factors")
print("Mean:", np.mean(high_res_bfactors))
print("Std:", np.std(high_res_bfactors))
print("Min:", np.min(high_res_bfactors))
print("Max:", np.max(high_res_bfactors))

print("Descriptive Statistics of Low-Res B-factors")
print("Mean:", np.mean(low_res_bfactors))
print("Std:", np.std(low_res_bfactors))
print("Min:", np.min(low_res_bfactors))
print("Max:", np.max(low_res_bfactors))

######Statistical Tests####
###T-test#####
t_stat, t_pval = stats.ttest_ind(high_res_bfactors,low_res_bfactors)
print("T-test")
print(f"t-statistic = {t_stat}")
print(f"p-value = {t_pval}")
print(f"Results: {'SIGNIFICANT' if t_pval < 0.05 else 'NOT SIGNIFICANT'} at p < 0.05")

###Mann-Whitney Test####
u_stat, u_pval = stats.mannwhitneyu(high_res_bfactors,low_res_bfactors, alternative='two-sided')
print(f"Mann-Whitney U test")
print(f"U-statistic = {u_stat}")
print(f"p-value     = {u_pval}")
print(f"Result: {'SIGNIFICANT' if u_pval < 0.05 else 'NOT significant'} at p < 0.05")

# Test 3: Check normality (Shapiro-Wilk)
_, p_norm_high = stats.shapiro(high_res_bfactors[:50])
_, p_norm_low  = stats.shapiro(low_res_bfactors[:50])
print("Normality test (Shapiro-Wilk)")
print(f"High-res: p = {p_norm_high} "
      f"({'normal' if p_norm_high > 0.05 else 'NOT normal'})")
print(f"Low-res:  p = {p_norm_low}  "
      f"({'normal' if p_norm_low  > 0.05 else 'NOT normal'})")

######Plots######
fig, axes = plt.subplots(1,3, figsize = (15,5))
sns.set_style("whitegrid")

###Kernel Density Estimate#### - smooth distribution curve
axes[0].set_title('B-factor Distribution (KDE)', fontsize=12)
sns.kdeplot(high_res_bfactors, ax=axes[0], label='High-res (<1.5Å)',
            color='steelblue', fill=True, alpha=0.4)
sns.kdeplot(low_res_bfactors,  ax=axes[0], label='Low-res (>3.0Å)',
            color='coral', fill=True, alpha=0.4)
axes[0].set_xlabel('Mean Cα B-factor (Å²)')
axes[0].set_ylabel('Density')
axes[0].legend()

###Plot 2:Boxplot####
axes[1].set_title('B-factor Boxplot', fontsize=12)
data = [high_res_bfactors,low_res_bfactors]
bp = axes[1].boxplot(data,  labels=['High-res\n(<1.5Å)', 'Low-res\n(>3.0Å)'],
                      patch_artist=True)
bp['boxes'][0].set_facecolor('steelblue')
bp['boxes'][1].set_facecolor('coral')
axes[1].set_ylabel('Mean Cα B-factor (Å²)')

# Add p-value annotation
axes[1].text(1.5, max(high_res_bfactors.max(), low_res_bfactors.max()),
             f'p = {u_pval}', ha='center', fontsize=11,
             color='red' if u_pval < 0.05 else 'black')

# Plot 3: Histogram overlay
axes[2].set_title('B-factor Histogram', fontsize=12)
axes[2].hist(high_res_bfactors, bins=20, alpha=0.6,
             color='steelblue', label='High-res (<1.5Å)')
axes[2].hist(low_res_bfactors,  bins=20, alpha=0.6,
             color='coral',     label='Low-res (>3.0Å)')
axes[2].set_xlabel('Mean Cα B-factor (Å²)')
axes[2].set_ylabel('Count')
axes[2].legend()

plt.suptitle('Do B-factors differ between high and low resolution structures?',
             fontsize=13)
plt.tight_layout()
plt.savefig('./results_today/w02d1_bfactor_stats.png', dpi=150)
plt.show()

####Save data to CSV###
data1_df = pd.DataFrame({
    'pdb_id' : high_res_ids[:len(high_res_bfactors)],
    'group' : 'high-res',
    'mean b-factor' : high_res_bfactors
    
})

data2_df = pd.DataFrame({
    'pdb_id' : low_res_ids[:len(low_res_bfactors)],
    'group' : 'low-res',
    'mean b-factor' : low_res_bfactors
    
})

final_data = pd.concat([data1_df,data2_df])
final_data.to_csv("./results_today/BFactor_results.csv", index = False)
print("Results Saved")
