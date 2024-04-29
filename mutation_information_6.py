"""Extracting mutation information from ensembl with ensembl rest api"""

# TODO: Change the name of the file (once its finished)

# Import modules
import requests
import json
import re

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


# TODO: Write 2 functions separately based on argument including/excluding introns?

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

# X and Y chromosomes
#ENST00000603986:c.1050C>T
#ENST00000603986:c.1223G>A
#ENST00000603986:c.2005G>A
#ENST00000603986:c.904C>T
#ENST00000603986:c.10G>A
#ENST00000603986:c.684G>A
#ENST00000383070:c.55C>T
#ENST00000603986:c.1621C>T
#ENST00000603986:c.987T>C
#ENST00000383070:c.217G>A

# 2 first X, 2 last Y
#hgvs_notation = ["ENST00000603986:c.1050C>T", "ENST00000603986:c.1223G>A", "ENST00000383070:c.55C>T", "ENST00000383070:c.217G>A"]

# Variables
#hgvs_notation = ["ENST00000558353:c.323G>A"]
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000560442:c.822C>T"]        #  2 working
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A"]        #  1 working and 2 not 
#hgvs_notation = ["ENST00000558353:c.323G>A", "ENST00000169551:c.609G>A"]         #  1 not and 2 working
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A"]     # 1 working 2 not 3 working
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A", "ENST00000431539:c.757C>T"]   # 1 working 2 not 3 working 4 not

#ENSG00000184831     x-transcript    APOO-gene
#ENST00000603986.6   x-transcript    CCDC120-gene
#ENST00000383070.2   y-transcript    SRY-gene

# TODO: Format these correct
#ENST00000603986:c.1050C>T, : NC_000023.11:g.49067337G>A    X
#ENST00000603986:c.1223G>A, : NC_000023.11:g.49067337G>A    x
#ENST00000383070:c.55C>T, : NC_000024.10:g.2787387C>T       Y
#ENST00000383070:c.217G>A, : NC_000024.10:g.2787387C>T      Y


# TODO: Do I need a control step for HGVS coding sequence?    - Last

# TODO: Install local DB?                                     - Last
# TODO: Add so that if reference base not matching,           - Last
#           change the reference to input base?


# TODO: Integrate with main program

# TODO: Cleanup function especially prints

# Solving failed inputs:
# TODO: Re-request the failed inputs
# TODO: extract positions from failed -> change reference base -> re-request
#               - Then either keep reference change or change back.
# TODO: later merge so failed that are rerun are appended to hgvs_genomic


### Main coding to genomic - request ###
# TODO: Fails to remove from list for the X and Y transcripts - Double for the X and Y -probably why
# TODO: Change 23 to X and 24 to Y
def get_hgvs_genomic(hgvs_input, url, headers):
    """
    Extracts and HGVS genomic (HGVSG) corresponding to the input HGVS coding (HGVSC).
    """
    
    ### Ensembl REST API (variant_recorder)
    # JSON data for the POST request
    data = json.dumps({"ids": hgvs_input})   
    
    # Send a POST request
    response = requests.post(url, headers=headers, data=data)
    
    # Storage
    hgvs_genomic = {}
    hgvs_failed = hgvs_input[:]                   # A copy of the input in which successfull inputs will be removed.
    
    # Check if the request was successful
    if response.status_code == 200:                                   
        
        # Parse the response JSON
        result = response.json()
        found_match = False                                              # Flag for matching HGVS coding
        
        #### Iterate over each HGVS coding  ####
        for i, allele_result in enumerate(result):                       # Enumerate adds index "i" to the object
            for _, variant_data in allele_result.items():                # Iterate over each dictionary in the list
                if type(variant_data) == dict:                           # Errors saved as lists, hits saved as dictionaries.
                    
                    # Iterate through each HGVS input to match with responses
                    for hgvs in hgvs_input:                                           # List of HGVS coding
                        
                        # Remove version number
                        base_hgvs_coding = hgvs.split(':')[0].split('.')[0]              # Remove version number
                        
                        # Check for both 'hgvsc' and 'hgvsg' in the variant details
                        if 'hgvsc' in variant_data and 'hgvsg' in variant_data:                                          # Check if both keys are present
                            
                            # Check if the input HGVS coding matches any returned HGVSC values (ignoring version)
                            if any(base_hgvs_coding in s.split(':')[0].split('.')[0] for s in variant_data['hgvsc']):    # Check if any HGVSC matches
                                
                                hgvs_genomic[hgvs] = variant_data['hgvsg'][0]     # Save successfull genomic HGVS
                                # TODO: The remove step fails for the X and Y transcripts
                                # TODO: ALL X and Y are printed double times, probably the reason for the failed removal
                                # TODO: Only fails when same transcript appears 2 times - Probably due to split of string -> Removes both with same transcript?
                                print(f"hgvs to be removed: {hgvs}\n")
                                hgvs_failed.remove(hgvs)                          # Remove successfull input from failed list
                                
                                found_match = True           
                                
                            #TODO: Here an else if it does not match
                else:
                    print("Input contains errors.")
                    
        if not found_match:
            print("No matching HGVS coding found.")
        
    else:
        print(f"Failed to retrieve data: {response.status_code}\n")
        print(f"{response.text}\n")
        
    return hgvs_genomic, hgvs_failed


### failed input list to genomic ###

def failed_to_genomic(hgvs_failed, url, headers):
    """Re-request the failed inputs and append to the genomic dictionary"""


# TODO: Double check numbers for XY
#           - X: 23  y: 24
#           - saved as numebrs or XY?
# TODO: Adding as a nested dictionary instead

### Genomic to mutation information ###
genomic_test = {"ENST00000169551:c.609G>A": "NC_000018.10:g.74158243G>A"}

def hgvs_converter(hgvs_genomic):
    """ Takes the HGVS dictionary as input and returns mutation information for VCF creation"""
    chr_info = {}
    
    pattern = re.compile(r'(\d{2})\..*\.(\d+)[ATCGU]')      # Patterns for chromosome and locus
    
    for key, genomic in hgvs_genomic.items():
        match = pattern.search(genomic)
        if match:
            chromosome = match.group(1).lstrip("0")         # lstrip() removes the potential "0" from the first position.
            locus = match.group(2)
            chr_info[key] = (chromosome, locus)
        else:
            print(f"No match found in: {genomic}")
            
    return chr_info




#### Test program ####
if __name__ == "__main__":
    
    
    # REST API URL and headers
    url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    #hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A", "ENST00000431539:c.757C>T"]   # 1 working 2 not 3 working 4 not
    
    ###   These are X an Y - First does not work, second work    ###
    
    #hgvs_notation = ["ENST00000603986:c.1050C>T", "ENST00000603986:c.1223G>A", "ENST00000383070:c.55C>T", "ENST00000383070:c.217G>A"] # X and Y
    hgvs_notation = ["ENST00000603986:c.1050C>T", "ENST00000383070:c.55C>T"] # X and Y
    
    
    # Result objects (genomic in dict, failed in list)
    hgvs_genomic_test, hgvs_failed_test = get_hgvs_genomic(hgvs_notation, url, headers)
    
    # Test print of the program
    print("\nHere starts the coding-genomic:\n")
    
    # Print the successful HGVS coding
    for key, value in hgvs_genomic_test.items():
        print(f"coding: {key}, genomic: {value}")
        
        
    print("\nHere starts the failed:\n")
    
    # Print the failed HGVS coding
    for failed in hgvs_failed_test:
        print(f"Failed: {failed}")
    
    ### genomic hgvs -> chromosome info
    genomic_test = {
        "ENST00000169551:c.609G>A": "NC_000018.10:g.74158243G>A",
        "ENST00000519026:c.1396G>A": "NC_000008.11:g.20145849C>T",
        "ENST00000603986:c.1223G>A": "NC_000023.11:g.49067337G>A",  # X - chromosome
        "ENST00000383070:c.55C>T": "NC_000024.10:g.2787387C>T",     # Y - chromosome
        }
    
    chr_info = hgvs_converter(genomic_test)
    
    for key, value in chr_info.items():
        print(f"key: {key}\nvalue: {value}")


