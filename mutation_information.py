"""Extracting mutation information from ensembl with ensembl rest api"""

# TODO: Change the name of the file

# Import modules
import requests
import sys


### Example mutations ###

#ENST00000558353:c.323G>A
#ENST00000560442:c.822C>T
#ENST00000431539:c.757C>T
#ENST00000169551:c.609G>A
#ENST00000692002:c.79G>A
#ENST00000519026:c.1396G>A
#ENST00000560442:c.4G>A
#ENST00000635950:c.388C>T
#ENST00000169551:c.37G>A
#ENST00000635950:c.1004C>T
#ENST00000645992:c.37C>T
#ENST00000477459:c.716G>A
#ENST00000379568:c.1149G>A
#ENST00000643996:c.407G>A
#ENST00000431539:c.124G>A
#ENST00000510305:c.349C>T
#ENST00000635950:c.1690G>A
#ENST00000711487:c.5212C>T
#ENST00000361725:c.705C>T

# Access variables
server = "https://rest.ensembl.org"
ext = "/variant_recoder/homo_sapiens"                                                # For now, homo sapiens as default
headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
#hgvs_notation = "ENST00000558353:c.323G>A"

### Defining functions to retrieve data from ensembl rest api ###
# TODO: Wrap in function
# TODO: change so it takes my program as input

# Make the POST request to the server
#r = requests.post(server+ext, headers=headers, data='{ "ids" : ["rs56116432", "rs1042779" ] }')
r = requests.post(server+ext, headers=headers, data='{ "ids" : ["ENST00000558353:c.323G>A"] }')

def get_mut_info():
    """Retrieves data from ensembl"""
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    print(repr(decoded))



### test program ###

if __name__ == "__main__":
    get_mut_info()