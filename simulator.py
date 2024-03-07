"""Processing FASTA files and mutational profiles to simulate mutational signatures."""


# Import Python libraries
import argparse                          # Library for parsing command-line arguments
import random                            # Library for generating random numbers
from Bio import SeqIO                    # Biopython library for handling FASTA files
import numpy as np                       # Library for numerical computing
from datetime import date                # Library for handling dates
import os                                # Library for interacting with the operating system
import re                                # Library for regular expressions

#TODO use or remove the following library

import sys                               # Library for system-specific parameters and functions



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
    
    parser.add_argument('-i', required=True, help="FASTA transcript files")                     # Adds -i FASTA file
    parser.add_argument('-f', required=True, help="Mutational profile")                         # Adds -f signature file
    parser.add_argument('-n', type=int, required=True, help="Number of simulated mutations")    # Adds -n number of mutations
    parser.add_argument('-r', type=int, required=True, help="Number of runs")                   # Adds -r number of runs
    
    args = parser.parse_args()                # Reads the command line arguments 
    return vars(args)                         # Returns the arguments as a dictionary



# Function to read a FASTA file and return sequences
def read_fasta(fasta_file):
    """
    Reads a FASTA file and returns a dictionary of sequences with their IDs.
    """
    
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):    # SeqIO.parse() reads FASTA
        sequences[record.id] = str(record.seq)         # Adds the ID and sequence to the dictionary
    return sequences

#endregion


# ====================
# FREQUENCY OPERATIONS
# ====================

#TODO Split these ones up into several blocks

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
            triplet, change, count = line.strip().split('\t')      # split('\t') splits the line where it finds a tab into a list of strings
            
            count = int(count)                                     # Converts the count string to an integer
            if triplet not in freq:                                # Checking if triplet is not in the dictionary
                freq[triplet] = {}                                 # If not, creates a new dictionary for the triplet as value
                
            freq[triplet][change] = {'count': count}               # Adds new key-value pair to the dictionary
            total_count += count
            
    # Calculating frequencies
    for triplet in freq:                                           # triplet is the key of the dictionary
        for change in freq[triplet]:                               # change is the key of the dictionary (value = triplet)
            freq[triplet][change]['freq'] = freq[triplet][change]['count'] / total_count   # Calculates the frequency of the change and adds it to the dictionary
            
    return freq



# Function to calculate triplet counts and probabilities
def calculate_triplet_counts(sequences):
    """
    Calculates triplet counts and probabilities for each sequence.
    Returns dictionaries containing triplet counts and probabilities.
    """
    
    triplet_count = {}  # Store triplet counts
    gene_length = {}    # Store gene lengths
    pos_in_gene = {}    # Store positions in genes
    counting = {}       # Store counting of triplets
    
    for transcript, sequence in sequences.items():                         # Iterates over the sequences
        gene_length[transcript] = len(sequence)                            # Adds the length of the sequence to the dictionary
        triplets = [sequence[i:i+3] for i in range(len(sequence) - 2)]     # Creates a list of triplets from the sequence
        count = {}
        for i, triplet in enumerate(triplets):                                         # Iterates over the triplets
            count[triplet] = count.get(triplet, 0) + 1                                 # Adds the triplet to the dictionary and counts it. If in the dictionary, adds 1 to the count. If not, adds the triplet to the dictionary count 1.
            pos_in_gene.setdefault(transcript, {}).setdefault(triplet, []).append(i)   # Adds the position of the triplet in the gene to the dictionary
            
        triplet_count[transcript] = {triplet: {'count': c, 'freq': c/(len(sequence)-2)} for triplet, c in count.items()}    # Adds the triplet count and frequency to the dictionary
        
        for triplet in triplet_count[transcript]:                                                          # Iterates over the triplets in the dictionary
            counting[triplet] = counting.get(triplet, 0) + triplet_count[transcript][triplet]['count']     # Adds the triplet to the counting dictionary and counts it. If in the dictionary, adds the count of the triplet. If not, adds the triplet to the dictionary count 1.
            
    return triplet_count, gene_length, pos_in_gene, counting               # Returns the dictionaries



# Function to calculate probabilities for each triplet in each transcript
def calculate_probabilities(triplet_count, counting):   # Takes the previous dictionaries as arguments
    """
    Calculates the probability for each triplet in each transcript.
    Returns a dictionary containing these probabilities.
    """
    
    probabilities = {}
    for transcript in triplet_count:
        for triplet in triplet_count[transcript]:
            prob = triplet_count[transcript][triplet]['count'] / counting[triplet]  # Calculates the probability of the triplet
            triplet_count[transcript][triplet]['prob'] = prob                       # Adds the probability to the dictionary
            
            if triplet not in probabilities:                                        # If the triplet is not in the dictionary
                probabilities[triplet] = {'name': [], 'prob': []}                   # Creates a new dictionary for the triplet
                
            probabilities[triplet]['name'].append(transcript)                       # Adds the transcript 
            probabilities[triplet]['prob'].append(prob)                             # Adds the probability
            
    return probabilities

#endregion

#TODO change name of Randomized...

# ====================
# RANDOMIZED OPERATIONS
# ====================

#region Random

# Function to perform random sampling based on frequencies
def random_sampling(frequencies, n_sim):
    
    elements = []
    probs = []
    
    for triplet in frequencies:
        for change in frequencies[triplet]:
            element = f"{triplet}_{change}"
            elements.append(element)
            probs.append(frequencies[triplet][change]['freq'])
            
    return np.random.choice(elements, size=n_sim, p=probs)



# Function to get a random position in a gene
def get_random_position_in_gene(posingene, transcript, triplet):
    """
    Gets a random position in a gene for a given transcript and triplet.
    """
    
    positions = posingene[transcript][triplet]   
    return random.choice(positions) + 1          

#endregion


# ====================
# OUTPUT OPERATIONS
# ====================

#region Output

#TODO cleanup when working

# Function to write output directories
def get_output_folders():
    """
    Returns output directory and subdirectory for today's date
    """
    
    # Metadata for the output directory name
    today = date.today().strftime('%Y%m%d')                                
    fasta_name = args['i'].split('/')[-1].split('.')[0]       # split('/')[-1] if  the file is given as a path. split('.')[0] returns the first part of the file name before the dot.
    signature_name = args['f'].split('/')[-1].split('.')[0]   # may have to change if adding standard signatures
    n_value = str(args['n'])
    
    print(f"Current working directory: {os.getcwd()}")
    print(f"Today's date (formatted): {today}")
    print(f"All directories: {os.listdir('.')}")
    
    #### output index #####
    
    # Get a list of all directories that start with the date string
    existing_dirs = [d for d in os.listdir('./output') if d.startswith('output_' + today)]
    print(f"Existing directories: {existing_dirs}")
    
    # Find the highest run index among the existing directories
    highest_run_index = 0
    for d in existing_dirs:
        match = re.search(r'\((\d+)\)', d) # Search for a number in parentheses
        print(f"Directory: {d}")  # Print the current directory name
        if match:
            print(f"Match: {match.group(1)}")  # 4. Print the current match
            highest_run_index = max(highest_run_index, int(match.group(1)))
    
    # output directories paths
    output_directories = f"output/output_{today}_({highest_run_index + 1})_{fasta_name}_{signature_name}_n{n_value}/"
    print(f"Output directories: {output_directories}")
    
    return output_directories



# Function to write output path
def get_output_path(directories_paths, run_index):
    """
    Returns the output path and sets file name with run index.
    """
    
    # Metadata for file name                                 # may have to change: .split('.')[0]  -  if file names usually contains more than one dot
    fasta_name = args['i'].split('/')[-1].split('.')[0]      # split('/')[-1] if  the file is given as a path. split('.')[0] returns the first part of the file name before the dot.
    signature_name = args['f'].split('/')[-1].split('.')[0]  # may have to change if adding standard signatures
    r_value = str(args['r'])
    
    # Output path
    #TODO add index to output folder
    output_path = (                                                                                             
        f"{directories_paths}"                                            # output folder for run
        f"{fasta_name}_{signature_name}_{run_index}_of_{r_value}.txt"     # Run file.
    )
    
    return output_path



# Function to write output to a file

#TODO rewrite to VCF 
#TODO remove intermediate

def write_output(output, output_directories):
    """
    Writes an intermediate file and then appends it to an output file.
    """
    
    # Output paths
    filepath_inter = "output/intermediate.txt"
    filepath_out = get_output_path(output_directories, i+1)     # i+1 is the run index
    
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


if __name__ == "__main__":                             # Checks if the script is executed as the main program or imported as module
    
    # Parse command-line arguments
    args = parse_arguments()                           # Calls the function to parse command-line arguments
    
    # Read input files and calculate frequencies
    freq = get_freq(args['f'])                         # Calculate frequencies from mutational profile
    sequences = read_fasta(args['i'])                  # Read FASTA file and return sequences
    
    # Calculate triplet counts, gene lengths, and probabilities
    triplet_count, gene_length, pos_in_gene, counting = calculate_triplet_counts(sequences)   # The function returns 4 dictionaries
    probabilities = calculate_probabilities(triplet_count, counting)
    
    # Perform random sampling based on triplet frequencies
    sampled_triplets = random_sampling(freq, args['n'])      # takes -n mutations as input
    
    #Creating output directory
    directories = get_output_folders()             # Returns the intended output directory path
    os.makedirs(directories, exist_ok=True)        # Creates the output directories if they does not exist
    
    # Main simulation part: 
    # iterating over sampled triplets and looping runs
    for i in range(args['r']):                            # Run the simulation -r times. Numbers series created with range().
        for run in sampled_triplets:                      # Each element in this list is a string representing a triplet and a substitution separated by an underscore "_"
            trip, sub = run.split("_")                    # split() splits the string  at the underscore. The parts are unpacked into the variables trip (triplet) and sub (substitution) 
            
            # probabilities is a dictionary containing information 
            # about each triplet, including names and probabilities.
            names = probabilities[trip]['name']                         # names is a list of transcripts
            probs = probabilities[trip]['prob']                         # probs is a list of probabilities
            
            # Randomly select a transcript based on probabilities
            selected_transcript = np.random.choice(names, p=probs)      # np.random.choice() returns a random transcript based on the probabilities. p=probs is the list of probabilities
            
            # Get a random position in the selected gene
            position = get_random_position_in_gene(pos_in_gene, selected_transcript, trip)   # Inside the function, it first retrieves the list of positions where the given triplet occurs in the given transcript from the posingene dictionary. Then it uses the random.choice() function to select a random position from this list.
            
            #TODO add output file (VCF format)
            
            # Output the result to screen
            ch = sub.split("/")              # split() splits the string at the slash and saves into a list. The parts are unpacked into the variables ch[0] and ch[1]. If sub was a string like "A/B", then ch would become the list ['A', 'B'].
            output_string = f"{selected_transcript}:c.{position}{ch[0]}>{ch[1]}"
            print(output_string)             # The printed string will have the format "selected_transcript:c.positionch[0]>ch[1]". If selected_transcript was "transcript1", position was 123, and ch was ['A', 'B'], the printed string would be "transcript1:c.123A>B"
            
            # Output to file
            write_output(output_string,directories)