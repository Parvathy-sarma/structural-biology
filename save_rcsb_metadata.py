###First goal is to load the metadata from RCSB. So RCSB have scripts in their query website to download the metadata, which we can process to get the metadata###
##For instance lets take the metadata for protein which resembles the UNIPROT ID P09874 (PARP1-protein). I chose this coz, I majorly work with PARP1 protein###

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
  "return_type": "polymer_entity"
}


response = requests.post(url, json=query)
results = response.json()


with open("./PDB.json", "w") as f:
    json.dump(results, f, indent=4)

print("Data Saved")
