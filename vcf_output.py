"""
Importing metadata and data from Ensembl database.
"""

import requests                     # requests is a Python library for making HTTP requests

# Todos to fix
#TODO add more fix argumentparse for transcript ID
#TODO Connect the files
#TODO double check REF base
#TODO remove end position
#TODO add any other info?
#TODO fix newline

#### Accessing ENSEMBL online

# Specify the transcript ID
transcript_id = "ENST00000367800"
#transcript_id_2 = "ENST00000327430"
transcript_pos = 2

# Function to retrieve the chromosome and start position of a transcript from the ENSEMBL REST API
def get_chr_pos(transcript_id, transcript_pos):
    """
    Retrieve the chromosome and start position of a transcript from the ENSEMBL REST API.
    """
    
    #TODO loop through transcript id to find the newest version.
    
    chr_info = {}
    
    # Make the GET request to the ENSEMBL REST API
    request_url = f"{url}/lookup/id/{transcript_id}?content-type=application/json" # Construct the API request URL - f-string is used to insert the transcript ID into the URL
    url = "https://rest.ensembl.org"                                               # Define the ENSEMBL API endpoint
    response = requests.get(request_url)                                           # function from the requests library to retrieve data from the specified URL. Stored as raw data to make JSON.
    
    # Check if the request was successful
    if response.status_code == 200:              # 200 is the HTTP status code for a successful request
        
        # Parse the JSON response
        data = response.json()      # function from the requests library to parse the JSON data
        chr_info["chromosome"] = data['seq_region_name']
        chr_info["start_pos"] = data['start']
        chr_info["end_pos"] = data['end']
        
        # Calculate chromosome position
        chr_info["chr_pos"] = chr_info["start_pos"] + transcript_pos - 1
    else:
        print("Failed to retrieve data from Ensembl:", response.status_code, "\n")
        
    return chr_info


# Test 
if __name__ == "__main__":
    transcript_info = get_chr_pos(transcript_id, transcript_pos)  # Call the function with the specified transcript ID
    
    for key, value in transcript_info.items():
        print(key, value)  # Print the chromosome and start position of the transcript
        
    print(transcript_info)
