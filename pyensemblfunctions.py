"""
Importing metadata and data from Ensembl database.
"""


import requests                     # requests is a Python library for making HTTP requests


#TODO add more fix argumentparse for transcript ID
#TODO Connect the files
#TODO double check REF base
#TODO remove end position
#TODO add any other info?

#### Accessing ENSEMBL online

# Specify the transcript ID
transcript_id = "ENST00000367800"

# Function to retrieve the chromosome and start position of a transcript from the ENSEMBL REST API
def get_chr_pos(transcript_id):
    """
    Retrieve the chromosome and start position of a transcript from the ENSEMBL REST API.
    """
    chr_info = {}
    
    # Define the ENSEMBL API endpoint
    url = "https://rest.ensembl.org"
    
    #TODO add argumentparse for transcript ID
    # Construct the API request URL
    request_url = f"{url}/lookup/id/{transcript_id}?content-type=application/json"  # f-string is used to insert the transcript ID into the URL
    
    # Make the GET request to the ENSEMBL REST API
    response = requests.get(request_url) # function from the requests library to retrieve data from the specified URL. Stored as raw data to make JSON.
    
    # Check if the request was successful
    if response.status_code == 200:              # 200 is the HTTP status code for a successful request
        
        # Parse the JSON response
        data = response.json()      # function from the requests library to parse the JSON data
        chr_info["chromosome"] = int(data['seq_region_name'])
        chr_info["start_pos"] = data['start']
        chr_info["end_pos"] = data['end']
    else:
        print("Failed to retrieve data from Ensembl:", response.status_code)
        
    return chr_info

transcript_info = get_chr_pos(transcript_id)  # Call the function with the specified transcript ID
for key, value in transcript_info.items():
    print(key, value)  # Print the chromosome and start position of the transcript

print(transcript_info)



#TODO better descriptions
#TODO double check positions

# Function to get mutated position
def set_mutated_pos(start_pos, transcript_pos):
    """
    Setting position on chromosome
    """
    
    mutated_pos = start_pos + transcript_pos - 1
    return mutated_pos

#Testing mutated mosition
transcript_pos = 2
mutated_pos = set_mutated_pos(transcript_info["start_pos"], transcript_pos)
print("Mutated position:", mutated_pos)