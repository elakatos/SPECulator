"""Extracting mutation information from ensembl with ensembl rest api"""

# TODO: Change the name of the file

# Import modules
import requests
import json


### Example mutations ###

#ENST00000645992:c.37C>T
#ENST00000477459:c.716G>A
#ENST00000379568:c.1149G>A
#ENST00000643996:c.407G>A
#ENST00000431539:c.124G>A
#ENST00000510305:c.349C>T
#ENST00000635950:c.1690G>A
#ENST00000711487:c.5212C>T
#ENST00000361725:c.705C>T


# TODO: Fix so it includes introns
# TODO: Write 2 functions separately based on argument iincluding/excluding introns
# TODO: Update so it only retrieves HGVS Genomic


### HGVS notations ###

# WORKS:

#hgvs_notation = "ENST00000560442:c.822C>T"
#hgvs_notation = "ENST00000169551:c.609G>A"
#hgvs_notation = "ENST00000692002:c.79G>A"
#hgvs_notation = "ENST00000519026:c.1396G>A"
#hgvs_notation = "ENST00000560442:c.4G>A"
#hgvs_notation = "ENST00000635950:c.388C>T"
#hgvs_notation = "ENST00000169551:c.37G>A"
#hgvs_notation = "ENST00000635950:c.1004C>T"


# DOESNT WORK:

#hgvs_notation = "ENST00000558353:c.323G>A"
#hgvs_notation = "ENST00000431539:c.757C>T"


# Variables
hgvs_notation = "ENST00000169551:c.609G>A"
url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
headers = {"Content-Type": "application/json", "Accept": "application/json"}
data = json.dumps({"ids": [hgvs_notation]})

# Send a POST request
response = requests.post(url, headers=headers, data=data)

def get_hgvsg():
    """Uses HGVS-coding to get HGVS-genomic"""
    
    # Check if the request was successful
    if response.status_code == 200:
        
        # Parse the response JSON
        result = response.json()
        print(json.dumps(result, indent=2))
        
    else:
        print(f"Failed to retrieve data: {response.status_code}")
        print(response.text)


# Test program
if __name__ == "__main__":
    get_hgvsg()