"""
Importing metadata and data from Ensembl database.
"""

import requests                     # requests is a Python library for making HTTP requests

# Todos to fix
#TODO Connect the files
#TODO double check REF base
#TODO remove end position
#TODO add any other info?
#TODO fix newline
#TODO fix printing faults when searching actual version

#### Accessing ENSEMBL online

# Transcript ID tests
#transcript_id = "ENST00000367800"      # Example existing  
#transcript_id = "ENST00000327430"      # Example not existing
#transcript_id = "ENST00000327430.3"    # Example existing with version

transcript_pos = 2

# Function to retrieve the chromosome and start position of a transcript from the ENSEMBL REST API
def get_transcript_info(transcript_id):
    """
    Retrieve the chromosome and start position of a transcript from the ENSEMBL REST API.
    """
    
    # Make the GET request to the ENSEMBL REST API
    url = "https://rest.ensembl.org"                                               # Define the ENSEMBL API endpoint    
    request_url = f"{url}/lookup/id/{transcript_id}?content-type=application/json" # Construct the API request URL - f-string is used to insert the transcript ID into the URL
    response = requests.get(request_url)                                           # function from the requests library to retrieve data from the specified URL. Stored as raw data to make JSON.
    
    # Check if the request was successful
    if response.status_code == 200:              # 200 is the HTTP status code for a successful request
        # Parse the JSON response
        chr_info = {}
        data = response.json()                   # function from the requests library to parse the JSON data
        chr_info["chromosome"] = data['seq_region_name']
        chr_info["start_pos"] = data['start']
        chr_info["end_pos"] = data['end']
    else:
        print("Failed to retrieve data from Ensembl:", response.status_code, "\n")
        
        return False
        
    return chr_info



# Function to find the newest version of the transcript ID if the input is not existing.
def find_transcript_version(transcript_id):
    """
    Find the newest version of the transcript ID and retrieve data.
    """
    
    try:
        for version in range(1, 11):
            transcript_id_with_version = f"{transcript_id}.{version}"
            chr_info = get_transcript_info(transcript_id_with_version)
                
            if chr_info:
                break
        print(transcript_id_with_version)
    
        return chr_info
    except:
        print("Transcript ID not found.")
        return False



# Calculates the chromosome position from the chromosome starting position.
def calc_chr_pos(chr_info, transcript_pos):
    """
    Calculate the chromosome position.
    """
    
    chr_info["chr_pos"] = chr_info["start_pos"] + transcript_pos - 1
    return chr_info


# Test program
if __name__ == "__main__":
    
    transcript_info = get_transcript_info(transcript_id)  # Call the function with the specified transcript ID
    
    if not transcript_info:
        transcript_info = find_transcript_version(transcript_id)
    
    if transcript_info:
        transcript_info = calc_chr_pos(transcript_info, transcript_pos)
        
        for key, value in transcript_info.items():
            print(key, value)  # Print the chromosome and start position of the transcript
        
    print(transcript_info)