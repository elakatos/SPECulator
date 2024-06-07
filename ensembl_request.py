"""Extracting mutation information from ensembl with ensembl rest api"""

# Import modules
import requests
import json
import re


# Test hgvs coding tests and extra print statements are included in the "ensembl_request_with_prints.py" file.

# Takes the list of hgvs_coding if it contains more than 200 items and splits it into smaller lists.
# Integrated into get_hgvs_genomic function
def split_list(list, size):
    """
    Splits the list into smaller lists.
    """
    
    # Return a list of sublists with the specified size
    return [list[i:i + size] for i in range(0, len(list), size)] 

# TODO: This function needs optimization. It loops every input multiple times.
# Takes a list of HGVS coding, url to ensembl,and Ensembl REST API (variant_recorder) 
# headers as input and returns a dictionary with the corresponding HGVS genomic (HGVSG) notation.
def get_hgvs_genomic(hgvs_input, url, headers, batch_size):
    """
    Extracts and HGVS genomic (HGVSG) corresponding to the input HGVS coding (HGVSC).
    """
    
    hgvs_input = list(set(hgvs_input))                                 # Remove duplicates
    
    # Check if input needs to be split into smaller lists
    if len(hgvs_input) > batch_size:                  
        sublists = split_list(hgvs_input, batch_size)              # Split the list into smaller lists   
    else:
        sublists = [hgvs_input]                                         # If the list is smaller than the maximum batch size, save it as a sublist
        
    # Storage
    hgvs_genomic = {}
    hgvs_failed = [(i, hgvs) for i, hgvs in enumerate(hgvs_input)]       # A list of tuples (index, value)
    index_map = {hgvs: i for i, hgvs in enumerate(hgvs_input)}           # A dictionary mapping hgvs to its original index
    
    #Start request loop for sublists
    for sublist in sublists:
        
        # TODO: Remove testprint when fixed
        # Testprint Before request
        print("\nSublist entering request")
        #for coding in sublist:
        #    print(coding)
            
        print("Accessing Ensembl REST API...")
        # Ensembl REST API (variant_recorder)
        data = json.dumps({"ids": sublist}) 
        
        # Execute request
        try:
            response = requests.post(url, headers=headers, data=data) 
        except requests.exceptions.RequestException as e:
            print(f"An error occurred while making the POST request: {e}")
            continue
        # Proceed with successful requests
        if response.status_code == 200:                                      # Status code 200 means successfull request                         
            
            result = response.json()                                         # Parse the response JSON
            # Testprint
            #print("\nOutput JSON from ensembl:")
            #print(json.dumps(result, indent=4))
            processed = set()                                                # A set to keep track of processed items
            found_match = False                                              # Flag for matching HGVS coding
            
            print("\nTesting matching coding/genomic:")
            # Iterate over each HGVS coding
            for allele in result:                       
                for _, variant_data in allele.items():                       # Iterate over each dictionary in the list
                    if type(variant_data) == dict:                           # Errors saved as lists, hits saved as dictionaries.
                        
                        # Iterate through each HGVS input to match with responses
                        for j, hgvs in enumerate(sublist):
                            base_hgvs_coding = hgvs.split(':')[0].split('.')[0]  # Extract transcript id
                            mutation_detail = hgvs.split(':')[1]  # Extract mutation detail
                            
                            for variant in variant_data['hgvsc']:
                                variant_transcript_id = variant.split(':')[0].split('.')[0]
                                variant_mutation_detail = variant.split(':')[1]
                                
                                if base_hgvs_coding == variant_transcript_id and mutation_detail == variant_mutation_detail:
                                    #print(f"Exact match found for: {hgvs}, adding: {variant_data['hgvsg'][0]}")
                                    hgvs_genomic[hgvs] = variant_data['hgvsg'][0]
                                    hgvs_failed.remove((index_map[hgvs], hgvs)) # Remove the matching HGVS coding
                                    processed.add(hgvs)
                                    # TODO: found_match only used for basic print. Elaborate or remove.
                                    found_match = True
                                    break                                       # Flag for matching HGVS coding
                                    
                    else:
                        print("Input contains errors.")
                        
            if not found_match:
                print("No matching HGVS coding found.")
                
        # Print error status code and message if the request failed
        else:
            print(f"Failed to retrieve data: {response.status_code}\n")
            print(f"{response.text}\n")
            
    return hgvs_genomic, hgvs_failed


# Converts the genomic HGVS notation to chromosome and locus information. Also returns the reference and alternative alleles.
def hgvs_converter(hgvs_genomic):
    """ 
    Takes the HGVS dictionary as input and returns mutation information in for VCF creation.
    Uses the regular expression module.
    """
    
    chr_info = {}
    pattern_genomic = re.compile(r'(\d{2})\..*\.(\d+)[ATCGU]')      # Patterns for chromosome number and locus
    pattern_key = re.compile(r'([ATCGU])>([ATCGU])')                # Pattern for letters before and after ">"
    
    # Save chromosome and locus information
    for coding, genomic in hgvs_genomic.items():
        match_genomic = pattern_genomic.search(genomic)
        match_coding = pattern_key.search(coding)
        
        if match_genomic and match_coding:
            # Check if the chromosome is somatic, X or Y
            if match_genomic.group(1) == "23":
                chromosome = "X"
            elif match_genomic.group(1) == "24":
                chromosome = "Y"
            else:
                chromosome = match_genomic.group(1).lstrip("0")    # lstrip() removes the potential "0" from the first position.
                
            locus = match_genomic.group(2)
            reference = match_coding.group(1)
            alternative = match_coding.group(2)
            chr_info[coding] = (chromosome, locus, reference, alternative)
        else:
            print(f"No match found in: {genomic}")
            
    return chr_info


#### Test program ####
if __name__ == "__main__":
    
    # REST API URL and headers
    url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # test hgvs coding inputs
    #hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A", "ENST00000431539:c.757C>T"]   # 1 working 2 not 3 working 4 not
    #hgvs_notation = ["ENST00000603986:c.1050C>T", "ENST00000603986:c.1223G>A", "ENST00000383070:c.55C>T", "ENST00000383070:c.217G>A"] # X and Y
    #hgvs_notation = ["ENST00000603986:c.1050C>T", "ENST00000383070:c.55C>T"] # X and Y
    hgvs_notation = ["ENST00000169551:c.609G>A", "ENST00000558353:c.323G>A", "ENST00000519026:c.1396G>A", "ENST00000431539:c.757C>T", "ENST00000603986:c.1050C>T", "ENST00000603986:c.1223G>A", "ENST00000383070:c.55C>T", "ENST00000383070:c.217G>A", "ENST00000383070:c.217G>A"]
    
    
    # Test importing hgvsc list
    path_to_file = "tests/example_transcript_list_test_signature_1_of_1.txt"
    
    #with open(path_to_file, "r") as file:
        #hgvs_notation = file.readlines()
        
    #hgvs_notation = [line.strip() for line in hgvs_notation]
    
    # Check input size:
    print(f"\nInput size: {len(hgvs_notation)}\n")
    
    
    # Call HGVS genomic function (genomic in dict, failed in list)
    hgvs_genomic_test, hgvs_failed_test = get_hgvs_genomic(hgvs_notation, url, headers)
    
    # Test print of the program
    print("\nMatching coding-genomic:\n")
    
    # Print the successful HGVS coding
    for key, value in hgvs_genomic_test.items():
        print(f"coding: {key}, genomic: {value}")
        
    print("\nHere starts the failed:\n")
    
    # Print the failed HGVS coding
    for failed in hgvs_failed_test:
        print(f"Failed: {failed}")
    
    
    # Call the HGVS converter function
    chr_info = hgvs_converter(hgvs_genomic_test)
    
    # Test prints for the chromosome information
    print("\nHere starts the chromosome information:\n")
    
    for key, value in chr_info.items():
        print(f"key: {key}\nvalue: {value}")