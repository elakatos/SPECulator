"""Processing FASTA files and mutational profiles to simulate mutational signatures."""


# Import built-in Python libraries
import argparse                          # Library for parsing command-line arguments
import random                            # Library for generating random numbers
from datetime import date                # Library for handling dates
import os                                # Library for interacting with the operating system
import re                                # Library for regular expressions
# TODO: use or remove the following library
import sys                               # Library for system-specific parameters and functions
import pickle                            # This library is particularly useful when you need to save complex Python data structures like lists, dictionaries, or class instances to a file that can be later retrieved

# Import external libraries
from Bio import SeqIO                    # Biopython library for handling FASTA files
import numpy as np                       # Library for numerical computing

# Import functions from other files
# TODO: fix typo! (now commented out as I don't have vcf)
# from vfc_output import get_chr_pos       # Import the function accessing ensembl - from the file vcf_output.py
from ensembl_request import get_hgvs_genomic, hgvs_converter 

# ====================
# INPUT OPERATIONS 
# ====================

#region Input


# Function to parse command-line arguments
def parse_arguments():
    """
    Parses command-line arguments and returns them as a dictionary.
    """
    
    parser = argparse.ArgumentParser(description="Script for processing FASTA files and mutational profiles.")   # Creates a new ArgumentParser object with a description
    
    parser.add_argument('-i', required=False, help="FASTA transcript files")                     # Adds -i FASTA file
    parser.add_argument('-t', required=False, help="Transcript list (file)")                     # Adds -t transcript list file
    parser.add_argument('-d', required=True, help="Database folder with pre-processed transcript pkl files")      # Adds -d database folder
    parser.add_argument('-c', required=False, help="Pre-processed transcript counts file")                        # Adds -c transcript counts file
    parser.add_argument('-f', required=True, help="Mutational profile")                         # Adds -f signature file
    parser.add_argument('-n', type=int, required=True, help="Number of simulated mutations")    # Adds -n number of mutations
    parser.add_argument('-r', type=int, required=True, help="Number of runs")                   # Adds -r number of runs
    
    args = parser.parse_args()                # Reads the command line arguments 
    return vars(args)                         # Returns the arguments as a dictionary
            # output example: {'i': 'input.fasta', 'f': 'profile.txt', 'n': 100, 'r': 10}

# Function to read a FASTA file and return sequences
def read_fasta_list(fasta_file):
    """
    Reads a FASTA file and returns a list of tuples with (ID, sequence) format.
    """
    
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(( ((record.id.split('|')[0]).split(' ')[0]).split('.')[0],str(record.seq)))
    return sequences                                   # Example output: [('ENSTxxx': 'ATGCGACTGATCGATCGTACG')]

# Function to read a FASTA file and return sequences
def read_fasta(fasta_file):
    """
    Reads a FASTA file and returns a dictionary of sequences with their IDs.
    """
    
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):    # SeqIO.parse() reads FASTA
        sequences[record.id] = str(record.seq)         # Adds the ID and sequence to the dictionary
        
    return sequences                                   # Example output: {'sequenceID': 'ATGCGACTGATCGATCGTACG',}

def read_transcript_list(transcript_list_file):
    """
    Reads in list of transcript IDs to consider
    """
    seqIDs = []
    with open(transcript_list_file, "r") as list_in:
        lines = list_in.readlines()
    seqIDs = [(l.strip().split('|')[0]).split('.')[0] for l in lines] # processing to ignore .version and other info
    return seqIDs                                                     # example output: ['ENST000001', 'ENST000002']

#endregion



# ====================
# FREQUENCY OPERATIONS
# ====================


#region Frequency

# Function to calculate frequencies from the mutational profile
def get_freq(profile_file):
    """
    Reads the mutational profile file and calculates frequencies.
    Returns a dictionary with triplet information and their frequencies.
    """
    freq = {}
    total_count = 0
    
    with open(profile_file, 'r') as file:                          # Opens the mutational count file
        for line in file:
            
            # triplet, change, count = The splited lines is assigned to these variables in order.
            triplet, change, count = line.strip().split('\t')      # Separates the three columns in the count file at "tab" and saves it to its respective variable 
            
            count = int(count)                                     # Converts the count string to an integer
            if triplet not in freq:                                # Checking if triplet is not in the dictionary - Might seem unecessary since each triplet should only occur once can prevent errors if the count files are wrong
                freq[triplet] = {}                                 # If not, creates a new dictionary for the triplet as value
                                                        # output example: {'CGA': {}} - The empty dictionary will be "change" 
            freq[triplet][change] = {'count': count}               # Adds new key-value pair to the dictionary
            total_count += count                        # output example: {'CGA': {'G/A': {'count': 10}}} - Triple nested dictionary
            
    # Calculating frequencies
    for triplet in freq:                                           # Iterates over the triplets in the dictionary               
        for change in freq[triplet]:     #TODO This is unecessary since there should only be one change per triplet
            freq[triplet][change]['freq'] = freq[triplet][change]['count'] / total_count   # Calculates the frequency of the change and adds it to the dictionary
            
    return freq                           # output example: {'CGA': {'G/A': {'count': 10, 'freq': 0.1}}}

# TODO: Redo so It calculates probabilities and not saves freq.

# Function to calculate triplet counts and probabilities
def calculate_triplet_counts(sequences):
    """
    Calculates triplet counts and probabilities for each sequence.
    Returns dictionaries containing triplet counts and probabilities.
    """
    
    gene_length = {}
    triplet_count = {}  
    pos_in_gene = {}    
    counting = {}       
    
    # Iterate sequences dictionary extracted from FASTA transcripts        Example:     {'sequenceID': 'ATGCGACTGATCGATCGTACG',}
    for transcript, sequence in sequences.items():                         # Iterates over the sequences dicrionary input
        gene_length[transcript] = len(sequence)                            # Adds the length of the sequence to the dictionary as value, ex: {'sequenceID': 21}
        triplets = [sequence[i:i+3] for i in range(len(sequence) - 2)]     # Creates a list of triplets from the sequence. [i:i+3] is the slice of the sequence.
                                                                        # -2 since it is triplets so lenght is two less than the full sequence. [i:i+3] gives 3 positions. 
        # Count triplets (per transcript and total), triplet frequencies (per transcript), store positions of triplets
        count = {}
        for i, triplet in enumerate(triplets):                                         # Iterates over the triplets. enumerate() function is used to get both the index (i) and the value (triplet) of each element.
            count[triplet] = count.get(triplet, 0) + 1                                 # counting the occurrences of each triplet. The get() method of the dictionary is used to retrieve the current count of the triplet (or 0 if the triplet is not yet in the dictionary), and then 1 is added to this count. The result is then stored back in the count dictionary with the triplet as the key.
            pos_in_gene.setdefault(transcript, {}).setdefault(triplet, []).append(i)   # Adds the position of the triplet in the gene to the dictionary. This line is storing the positions of each triplet in the pos_in_gene dictionary. The setdefault() method is used to ensure that there is a dictionary for each transcript and a list for each triplet. The index i is then appended to the list of positions for the triplet. The index i is then appended to the list of positions for the triplet. This results in a nested dictionary structure where the keys are the transcript values, each transcript maps to another dictionary where the keys are the triplet values, and each triplet maps to a list of positions.
                        # example: {'transcript1': {'ATG': [0, 5, 10], 'TGC': [1, 6]}, 'transcript2': {'ATG': [0, 3, 7], 'TGC': [1]}}
        triplet_count[transcript] = {triplet: {'count': c, 'freq': c/(len(sequence)-2)} for triplet, c in count.items()}    # Adds the triplet count and frequency to the dictionary. This line of code is creating a dictionary comprehension that maps each triplet to another dictionary. This inner dictionary contains two keys: 'count' and 'freq'.
                        # example: {'transcript1': {'ATG': {'count': 3, 'freq': 0.3}, 'TGC': {'count': 2, 'freq': 0.2}}}
        for triplet in triplet_count[transcript]:                                                          # Iterates over the triplets in the dictionary
            counting[triplet] = counting.get(triplet, 0) + triplet_count[transcript][triplet]['count']     # Adds the triplet to the counting dictionary and counts it. If in the dictionary, adds the count of the triplet. If not, adds the triplet to the dictionary count 1. The overall effect of this loop is to calculate the total count of each triplet across all transcripts and store these counts in the counting dictionary.
                        # example: {'ATG': 3, 'TGC': 2}
    return triplet_count, gene_length, pos_in_gene, counting               


def process_triplet_positions(sequences, db_folder, triplet_counts={}):
    """
    Calculates triplet counts and probabilities for each sequence.
    Saves triplet positions within each transcript as it processes them.
    Returns a dictionary containig triplet counts.
    Triplets counts can be supplied (if one wants to append to a previously computed dictionary) or started from empty
    """
    
    for (transcript, sequence) in sequences:                               # Iterates over the sequences
        triplets = [sequence[i:i+3] for i in range(len(sequence) - 2)]     # Creates a list of triplets from the sequence
        count = {}
        pos_in_transcript = {}
        
        for i, triplet in enumerate(triplets):                                         # Iterates over the triplets
            count[triplet] = count.get(triplet, 0) + 1                                 # Adds the triplet to the dictionary and counts it. If in the dictionary, adds 1 to the count. If not, adds the triplet to the dictionary count 1.
            pos_in_transcript.setdefault(triplet, []).append(i)                        # Adds the position of the triplet in the gene to the (transcript-specific) dictionary
        del triplets      # Delete triplets to help reduce memory usage
        
        dbname = db_folder+'/'+transcript+'.pkl'                                      # Used to construct a file path for a pickle file that is specific to each transcript.
        with open(dbname, 'wb') as saved_dict:
            pickle.dump(pos_in_transcript, saved_dict, pickle.HIGHEST_PROTOCOL)       # Pickle the processed triplet positions dictionary using the highest protocol available. pickle.HIGHEST_PROTOCOL is a constant that represents the highest protocol number available in your Python interpreter. Using the highest protocol allows you to pickle more kinds of Python objects
        [triplet_counts.setdefault(triplet, [[],[]])[1].append(c) for triplet, c in count.items()]  # Add triplet count info to total counts. If info on triplet does not exist yet, then initialises it as two lists
        [triplet_counts[triplet][0].append(transcript) for triplet in count.keys()]
        
    return triplet_counts             # Returns one count dictionary

def compute_triplet_counts(seqIDs, db_folder):
    """
    Calculates triplet counts for given list of transcripts, when transcripts are already in database
    """
    triplet_counts={}    
    
    for transcript in seqIDs:                                              # Iterates over the sequences
        with open(db_folder+'/'+transcript+'.pkl', 'rb') as read_db:       # Pathing together the file path for the pickle file that is specific to each transcript. Opens the pickle file in read mode. The 'r' stands for read mode, and the 'b' stands for binary mode. So 'rb' means the file is opened for reading in binary mode. This is typically used for non-text files like images or executable files, or in this case, pickle files.
            posingene = pickle.load(read_db)                               # Read pickled transcript dictionary
        [triplet_counts.setdefault(triplet, [[],[]])[1].append(len(pos)) for triplet, pos in posingene.items()]  # Add triplet count info to total counts. If info on triplet does not exist yet, then initialises it as two lists
        [triplet_counts[triplet][0].append(transcript) for triplet in posingene.keys()]
                    
    return triplet_counts 

# Function to calculate probabilities for each triplet in each transcript
def calculate_probabilities(triplet_count, counting):   # Takes the previous dictionaries as arguments
    """
    Calculates the probability for each triplet in each transcript.
    Returns a dictionary containing these probabilities.
    """
    
    probabilities = {}
    
    for transcript in triplet_count:
        for triplet in triplet_count[transcript]:                                   # iterates the transcript in the first level and the triplet in the second level
            prob = triplet_count[transcript][triplet]['count'] / counting[triplet]  # Calculates the probability of the triplet - Takes the transcript triplet count and divides it by the total count
            triplet_count[transcript][triplet]['prob'] = prob                       # Adds the probability to the dictionary - This gives the probability of the triplet appearing in this specific transcript given its overall frequency in all transcripts.
                        # Example: {'transcript1': {'triplet1': {'count': 10, 'freq': 0.5, 'prob': 0.2}}, 'transcript2': {'triplet1': {'count': 20, 'freq': 0.4, 'prob': 0.4}}}
            if triplet not in probabilities:                                        # If the triplet is not in the dictionary
                probabilities[triplet] = {'name': [], 'prob': []}                   # Creates a new dictionary for the triplet with two empty lists as values
                
            probabilities[triplet]['name'].append(transcript)                       # Adds the transcript 
            probabilities[triplet]['prob'].append(prob)                             # Adds the probability
                        # Example: {'CGA': {'name': ['transcript1', 'transcript2'], 'prob': [0.2, 0.4]}}
    return probabilities

#endregion



# ====================
# RANDOMIZED OPERATIONS
# ====================

#region Random

# Function to perform random sampling based on frequencies
def random_sampling(frequencies, n_sim):
    """
    Takes a dictionary of frequencies (freq) and performs random sampling (args=n times) based on these frequencies.
    Returns a list of sampled elements.
    """
    
    elements = []
    probs = []
    
    for triplet in frequencies:                                 # From get_freq() - # Example: {'CGA': {'G/A': {'count': 10, 'freq': 0.1}}}
        for change in frequencies[triplet]:                     # Iterates over the substitutions (second level) in the dictionary
            element = f"{triplet}_{change}"                     # Creates a string with the triplet and the substitution separated by an underscore
            elements.append(element)                            # Adds the string to the elements list
            probs.append(frequencies[triplet][change]['freq'])  # Adds the frequency to the probs list. This frequency will be used as the probability of selecting the corresponding element during the random sampling.
                                # Examples: ['CGA_G/A', 'CGA_C/A'], [0.1, 0.2]
    return np.random.choice(elements, size=n_sim, p=probs)      # This is an array of n_sim elements, each of which is a string from the elements list, selected randomly based on the probabilities in the probs list. The exact contents of the array will vary each time the function is called
                                # Example: array(['CGA_G/A', 'CGA_G/A', 'CGA_C/A', ..., 'CGA_G/A', 'CGA_C/A', 'CGA_G/A'])


# Function to get a random position in a gene
def get_random_position_in_gene(posingene, transcript, triplet):
    """
    Gets a random position in a gene for a given transcript and triplet. 
    Returns a random position that will be mutated
    """
    # This function returns a random position in a gene where a given triplet occurs in a given transcript.
    # Input: pos_in_gene{} - {'transcript1': {'ATG': [0, 5, 10], 'TGC': [1, 6]}, 'transcript2': {'ATG': [0, 3, 7], 'TGC': [1]}}
    positions = posingene[transcript][triplet]   # Retrieves the list of positions where the given triplet occurs in the given transcript from the dictionary.
    return random.choice(positions) + 1          # A random position is selected from this list using the random.choice() function. +1 is added since the triplet start at the first base and it is the second that will be mutated. 

def get_transcript_position(transcript, triplet, db_folder):
    """
    Loads a transcript and randomly chooses a position with a given triplet
    """
    with open(db_folder+'/'+transcript+'.pkl', 'rb') as read_db:
        posingene = pickle.load(read_db)
    
    positions = posingene[triplet]   
    return random.choice(positions) + 1 

#endregion


# ====================
# OUTPUT OPERATIONS
# ====================

#region Output


# Function to write output directories
def get_output_folders(fasta_name):
    """
    Returns output directory and subdirectory for today's date
    """
    
    # Metadata for the output directory name
    today = date.today().strftime('%Y%m%d')                                
    signature_name = '.'.join(args['f'].split('/')[-1].split('.')[:-1])   # may have to change if adding standard signatures.
    n_value = str(args['n'])                                              # makes a string of the number of mutations argument.
    
    #### output index #####
    
    # Get a list of all directories that start with the date string
    existing_dirs = [d for d in os.listdir('./output') if d.startswith('output_' + today)]  # Creates a new list of names of files and directories in ./output that start with 'output_' followed by today.
    
    # Find the highest run index among the existing directories
    highest_run_index = 0
    
    for d in existing_dirs:
        match = re.search(r'\((\d+)\)', d) # From regular expressions, first argument = pattern, second argument = string to be searched. r'\((\d+)\)' is the pattern. \d+ is a regular expression pattern that matches one or more digits. The parentheses are used to capture the digits as a group. earches the string "d" for a pair of parentheses enclosing one or more digits
        if match:
            highest_run_index = max(highest_run_index, int(match.group(1))) # The group() method returns the part of the string where there was a match. The max() function is used to find the highest run index among the existing directories.
    
    # output directories paths
    output_directories = f"output/output_{today}_({highest_run_index + 1})_{fasta_name}_{signature_name}_n{n_value}/"
    
    return output_directories



# TODO: change so it takes fasta_name and signature_name as arguments

# Function to write output path
def get_output_path(directories_paths, run_index, fasta_name):
    """
    Returns the output path and sets file name with run index.
    """
    
    # Metadata for file name                                       
    signature_name = '.'.join(args['f'].split('/')[-1].split('.')[:-1])  
    r_value = str(args['r']) 
    
    # Output path
    output_path = (                                                                                             
        f"{directories_paths}"                                            # output folder for run
        f"{fasta_name}_{signature_name}_{run_index}_of_{r_value}.txt"     # Run file.
    )
    
    return output_path



# Function to write output to a file

# TODO: rewrite to VCF 
# TODO: remove intermediate

def write_output(output, output_directories, fasta_name):
    """
    Writes an intermediate file and then appends it to an output file.
    """
    
    # Output paths
    filepath_inter = "output/intermediate.txt"
    filepath_out = get_output_path(output_directories, i+1, fasta_name)     # i+1 is the run index
    
    # Creating or overwriting intermediate file
    with open(filepath_inter, "w") as f:
        f.write(output + '\n')
        
    # Appending intermediate file to output file
    with open(filepath_inter, "r") as rf:       # Read mode intermediate file
        with open(filepath_out, "a") as wf:     # Append mode for output file
            for line in rf:
                wf.write(line)

#endregion


# ====================
# MAIN PROGRAM
# ====================


if __name__ == "__main__":                                                  # Checks if the script is executed as the main program or imported as module. If the script is executed as the main program, the code block is executed, not if it is imported as a module.
    
    # Parse command-line arguments
    args = parse_arguments()                                                # Calls the function to parse command-line arguments and return the arguments as a dictionary
    db_folder = args['d']                                                   # Assign database folder to a variable    
    # TODO: Check if db_folder exist, otherwise, create it.
    
    # Read input files and calculate frequencies
    freq = get_freq(args['f'])                                              # Calculate frequencies from mutational profile file and store in freq dictionary
    if args['i'] is not None:     # if a fasta file is provided
        sequences = read_fasta_list(args['i'])                              # Read FASTA file and return list of tuples of IDs and sequences
        triplet_counts = process_triplet_positions(sequences, db_folder)    # In this case new database is processed based on fasta
        run_name = '.'.join(args['i'].split('/')[-1].split('.')[:-1])       # split('/')[-1] - Picks the string after the last "/". split('.')[:-1] and '.'.join() removes the file type and joins the string.
    elif args['t'] is not None:   # if a list of transcripts is provided
        transcripts = read_transcript_list(args['t'])
        triplet_counts = compute_triplet_counts(transcripts, db_folder)
        run_name = '.'.join(args['t'].split('/')[-1].split('.')[:-1])       # split('/')[-1] - Picks the string after the last "/". split('.')[:-1] and '.'.join() removes the file type and joins the string.
    else:                         # use all entries in pre-processed count file
        with open(args['c'], 'rb') as read_db:     # Read in saved count file
            triplet_counts = pickle.load(read_db)
    # TODO: add checks for ensuring database and triplet_counts contain the same transcripts
    
    # Calculate triplet counts, gene lengths, and probabilities
    # triplet_count, gene_length, pos_in_gene, counting = calculate_triplet_counts(sequences)   # The function returns 4 dictionaries
    # probabilities = calculate_probabilities(triplet_count, counting)                          # Calculate probabilities for each triplet in each transcript and store in probabilities dictionary 
    
    # Perform random sampling based on triplet frequencies
    sampled_triplets = random_sampling(freq, args['n'])      # takes -n mutations as input. Returns a list of sampled elements.
    
    #Creating output directory
    directories = get_output_folders(run_name)             # Returns the intended output directory path as a string.
    os.makedirs(directories, exist_ok=True)                # Creates the output directories if they do not exist already.
    
    
    #### Main simulation part ####: 
    
    # Iterating over sampled triplets and looping runs
    for i in range(args['r']):                            # Run the simulation -r times. Numbers series created with range(). The loop will run the number of times specified by the -r argument.
        
        hgvsc_list = []                                   # Create an empty list to store the output strings
        
        for run in sampled_triplets:                      # Each element in this list is a string representing a triplet and a substitution separated by an underscore "_". The loop iterates over these strings. 
            trip, sub = run.split("_")                    # split() splits the string  at the underscore. The parts are unpacked into the variables trip (triplet) and sub (substitution). If run was a string like "CGA_G/A", then trip would become "CGA" and sub would become "G/A".
            
            # probabilities is a dictionary containing information about each triplet, including names and probabilities.
            # names = probabilities[trip]['name']         # names is a list of transcripts
            # probs = probabilities[trip]['prob']         # probs is a list of probabilities
            names = triplet_counts[trip][0]
            counts = triplet_counts[trip][1]              # counts are used as weights to not have to compute probabilities (faster)
            
            # TODO: isnt this already performed in a function?
            # Randomly select a transcript based on probabilities
            # selected_transcript = np.random.choice(names, p=probs)       # np.random.choice() returns a random transcript based on the probabilities. p=probs is the list of probabilities
            selected_transcript = random.choices(names, weights=counts)[0] # using random.choices (since faster), returns a randomly selected one weighted by the number of times each transcript has the triplet
            
            # Get a random position in the selected gene
            # position = get_random_position_in_gene(pos_in_gene, selected_transcript, trip)   # Inside the function, it first retrieves the list of positions where the given triplet occurs in the given transcript from the posingene dictionary. Then it uses the random.choice() function to select a random position from this list.
            position = get_transcript_position(selected_transcript, trip, db_folder)           # Reads in the corresponding transcript dictionary file, and chooses a random position where that triplet can be found
            
            # Output the result to screen
            ch = sub.split("/")                          # split() splits the string at the slash and saves into a list. The parts are unpacked into the variables ch[0] and ch[1]. If sub was a string like "A/B", then ch would become the list ['A', 'B'].
            output_string = f"{selected_transcript}:c.{position + 1}{ch[0]}>{ch[1]}"
            hgvsc_list.append(output_string)             # Append the output string to the list hgvsc_list
            
            #print(f"\n{output_string}")                  # The printed string will have the format "selected_transcript:c.positionch[0]>ch[1]". If selected_transcript was "transcript1", position was 123, and ch was ['A', 'B'], the printed string would be "transcript1:c.123A>B"
            
            # Output to file
            write_output(output_string,directories, run_name)
            
            
        # TODO: Clean up all prints
        # test for hgvsc_list
        
        #print("\nHGVSC test:\n")
        #for hgvsc in hgvsc_list:
            #print(hgvsc)
            
        #### Retrieve chromosome and chromosome position from ENSEMBL REST API ####
        
        url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        
        hgvs_genomic, hgvs_failed = get_hgvs_genomic(hgvsc_list, url, headers)        # Calls the function to get the HGVS genomic notation from the REST API. The function returns 1 dictionary, 1 list: hgvs_genomic and hgvs_failed.
        
        #print("\nMatching coding-genomic:\n")
        
        # Print the successful HGVS coding
        #for key, value in hgvs_genomic.items():
            #print(f"coding: {key}, genomic: {value}")
            
        print("\nHere starts the failed:\n")
        
        # Print the failed HGVS coding
        for failed in hgvs_failed:
            print(f"Failed: {failed}")
            
        # Call the HGVS converter function
        chr_info = hgvs_converter(hgvs_genomic)
        
        # Test prints for the chromosome information
        #print("\nHere starts the chromosome information:\n")
        
        #for key, value in chr_info.items():
            #print(f"key: {key}\nvalue: {value}")