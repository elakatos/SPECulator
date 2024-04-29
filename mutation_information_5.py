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


# Variables
#hgvs_notation = ["ENST00000558353:c.323G>A"]
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000560442:c.822C>T"]        #  2 working
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A"]        #  1 working and 2 not 
#hgvs_notation = ["ENST00000558353:c.323G>A", "ENST00000169551:c.609G>A"]         #  1 not and 2 working
#hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A"]     # 1 working 2 not 3 working
hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A", "ENST00000431539:c.757C>T"]   # 1 working 2 not 3 working 4 not



# TODO: Do I need a control step for HGVS coding sequence?    - Last
# TODO: Is it possible to exclude all unecessary information  - Last
#           right from the beginning to save time?
# TODO: Install local DB?                                     - Last
# TODO: Add so that if reference base not matching,           - Last
#           change the reference to input base?

# TODO: Doesnt work if only 1 input and that one fails.
# TODO: Double check list of coding IDs
#               - Working successfully saved in dictionary
#                        - ONLY 1 "NOT WORKING" registered
#                           - Fix and save in dictionary
#                  
# TODO: Double check numbers for XY
# TODO: Transform HGVS genomic to relevant data.
# TODO: Make it save errors in a dictionary and return
# TODO: Integrate with main program

# TODO: Cleanup function especially prints
# TODO: Bryt ut i flera funktioner

def get_hgvs_genomic(hgvs_input):
    """
    Extracts and HGVS genomic (HGVSG) corresponding to the input HGVS coding (HGVSC).
    """
    
    # TODO: shouold these be in the function? main program file?
    # Ensembl REST API (variant_recorder)
    url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = json.dumps({"ids": hgvs_input})
    
    # Send a POST request
    response = requests.post(url, headers=headers, data=data)
    
    # Dictionaries to save in 
    hgvs_genomic = {}
    hgvs_failed = {}
    
    # Check if the request was successful
    if response.status_code == 200:
        
        # Parse the response JSON
        result = response.json()
        
        #### Iterate over each HGVS coding  ####
        found_match = False                                                          # Flag for matching HGVS coding
        
        for i, allele_result in enumerate(result):         # Enumerate adds index "i" to the object
            for _, variant_data in allele_result.items():  # Iterate over each dictionary in the list
                
                if type(variant_data) == list:  # Errors saved as lists, hits saved as dictionaries.
                    
                    # Add failed to dictionary
                    hgvs_failed[hgvs_input[i]] = variant_data
                    
                    # Test for failed HGVS coding
                    
                    
                # Continues with successfull HGVS coding
                else:
                    # Iterate through each HGVS input to match with responses
                    for hgvs in hgvs_input:                                           # List of HGVS coding
                        # Modify input notation for comparison without version number
                        base_hgvs_coding = hgvs.split(':')[0].split('.')[0]              # Remove version number
                        
                        # Check for both 'hgvsc' and 'hgvsg' in the variant details
                        if 'hgvsc' in variant_data and 'hgvsg' in variant_data:                                          # Check if both keys are present
                            # Check if the input HGVS coding matches any returned HGVSC values (ignoring version)
                            if any(base_hgvs_coding in s.split(':')[0].split('.')[0] for s in variant_data['hgvsc']):    # Check if any HGVSC matches
                                print(f"Input: {hgvs} -> HGVS genomic (HGVSG): {variant_data['hgvsg'][0]}")
                                
                                hgvs_genomic[hgvs] = variant_data['hgvsg'][0]
                                
                                found_match = True           
                                
                                # Flag for matching HGVS coding
                            #TODO: Here an else if it does not match
        
        if not found_match:
            print("No matching HGVS coding found.")
        
    # TODO: This should also be saved in the dictionary
    else:
        print(f"Failed to retrieve data: {response.status_code}")
        print(response.text)
        
    return hgvs_genomic, hgvs_failed


# Test dictionary
genomic_to_mut_info = {"ENST00000169551:c.609G>A": "NC_000018.10:g.74158243G>A"}

def hgvs_converter(genomic_to_mut_info):
    """
    Takes the HGVS dictionary as input and returns 
    mutation information for VCF creation
    """
    



#### Test program ####
if __name__ == "__main__":
    hgvs_genomic_test, hgvs_failed_test = get_hgvs_genomic(hgvs_notation)
    
    # Test print of the program
    print("\nHere starts the coding-genomic:\n")
    
    for key, value in hgvs_genomic_test.items():
        print(f"coding: {key}, genomic: {value}")
        
        
    print("\nHere starts the failed:\n")
    
    for key, values in hgvs_failed_test.items():
        print(f"Key: {key}\n")
        for value in values:
            print(f"Value: {value}\n")


