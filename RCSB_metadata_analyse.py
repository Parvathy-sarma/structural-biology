###Loading metadata of PDBs to look into their resolution, deposition year and other parameters###
##After loading the metadata of around 100 PDBs, we would load those PDB files locally. Once acheived, next step will be to analyse those metadata##

import pandas as pd
import matplotlib.pyplot as plt
import requests
import json

url = "https://search.rcsb.org/rcsbsearch/v2/query"

query = {
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "P09874",
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "operator": "exact_match",
          "value": "UniProt",
          "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
        }
      }
    ]
  },
  "request_options": {
    "paginate": {
      "start": 0,
      "rows": 115
    }
  },
  "return_type": "polymer_entity"
}


response = requests.post(url, json=query)
results = response.json()

###Read the metadata, that is fetch the PDB files listed in the json file from the RCSB site and save it/process it as required###

##Extract PDB_IDs##

pdb_ids =[r['identifier'] for r in results['result_set']]
print(f"{len(pdb_ids)} PDBs retrieved")

###The above code is the compact form of the one given below###
# pdb_ids = []
# for r in results['result_set']:
#     pdb_ids.append(r['identifier'])
# print(pdb_ids)

###Fetch details for PDB. Usually what we do is paste the PDB id in search box and download the data that is required
###We are doing the same here, we take the PDBid, paste it in the URL###

record = []
for pdb_id in pdb_ids[:len(pdb_ids)]:
    entry = pdb_id.split('_')[0]
      # Create the API URL for this specific PDB entry
    pdb_url = f"https://data.rcsb.org/rest/v1/core/entry/{entry}"

    # Send a request to the server to fetch data (like opening a webpage programmatically)
    r = requests.get(pdb_url)

    #### Request module acts like wget ##
    
    # Check if the request worked successfully (200 = OK)
    if r.status_code == 200:

        # Convert the response (JSON text) into a Python dictionary
        data = r.json()

        # Add extracted information to the 'record' list
        record.append({
            
            'pdb_id' : entry,
            
            # Get resolution (takes first value from a list, or None if missing)
            'resolution': data.get('rcsb_entry_info', {}).get('resolution_combined', [None])[0],

            # Get experimental method (e.g., X-ray, Cryo-EM)
            'method': data.get('rcsb_entry_info', {}).get('experimental_method', None),

            # Get deposit year (extract first 4 characters from date string)
            'year': data.get('rcsb_accession_info', {}).get('deposit_date', 'unknown')[:4],

            # Get organism taxonomy ID (may be None if not available)
            'organism': data.get('rcsb_entry_info', {}).get('organism_taxid', None)
        })
with open ("./records.json", "w") as f:
    json.dump(record, f, indent=4)

print("Data saved")


###Analysing this data, is first making a data frame or make it into a table format##

df = pd.DataFrame(record)

####Make it into a CSV file####
df.to_csv("./metada.csv", index="False")

###Plot histogram of Resolution data###
df['resolution'].dropna().plot(kind = "hist", bins = 20, title = 'Resolution distribuition', xlabel= "Resolution (Ang)", edgecolor='black', linewidth=1.5)
plt.savefig("./resolution.png", format = 'png')
plt.show()

####Plot distribuition year

df.groupby('year')['pdb_id'].count().plot(
    kind = 'bar', title = 'Structure per year', xlabel = "Year", ylabel = "Structure", edgecolor='black', linewidth=1.5)
plt.savefig("./Structure_deposited.png", format = 'png')
plt.show()
