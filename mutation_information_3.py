"""Extracting mutation information from ensembl with ensembl rest api"""

# TODO: Change the name of the file (once its finished)

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


# TODO: Write 2 functions separately based on argument including/excluding introns

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

# TODO: Change so user inputs ID
# TODO: Double check list of IDs
# TODO: Do I need a control step for HGVS coding sequence?
# TODO: Is it possible to exclude all unecessary information right from the beginning to save time?
# TODO: Install local DB?
# TODO: Integrate with main program
# TODO: Transform HGVS genomic to relevant data.

def get_hgvsg():
    """Extracts and prints HGVS genomic (HGVSG) corresponding to the input HGVS coding (HGVSC)."""
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the response JSON
        result = response.json()
        
        # Iterate over each allele result
        for allele_result in result:
            for _, variant_data in allele_result.items():
                
                # Check for both 'hgvsc' and 'hgvsg' in the variant details
                if 'hgvsc' in variant_data and 'hgvsg' in variant_data:
                    
                    # Modify input notation for comparison without version number
                    base_hgvs_coding = hgvs_notation.split(':')[0].split('.')[0]
                    
                    # Check if the input HGVS coding matches any returned HGVSC values (ignoring version)
                    if any(base_hgvs_coding in s.split(':')[0].split('.')[0] for s in variant_data['hgvsc']):
                        print(f"HGVS genomic (HGVSG): {variant_data['hgvsg'][0]}")  # Assuming only one HGVSG is returned
                        
                        return
                    
        print("No matching HGVS coding found.")
        
    else:
        print(f"Failed to retrieve data: {response.status_code}")
        print(response.text)


# Test program
if __name__ == "__main__":
    get_hgvsg()