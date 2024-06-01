"""Processing FASTA files and mutational profiles to simulate mutational signatures."""

# TODO: Remove modules and functions that are not used in the main program
# TODO: is "get_random_position_in_gene()" not used at all? Eszter has commented it out :)
# Import built-in Python libraries
import argparse                          # Library for parsing command-line arguments
import random                            # Library for generating random numbers
from datetime import date                # Library for handling dates
import os                                # Library for interacting with the operating system
import re                                # Library for regular expressions
import pickle                            # This library is particularly useful when you need to save complex Python data structures like lists, dictionaries, or class instances to a file that can be later retrieved

# Import external libraries
from Bio import SeqIO                    # Biopython library for handling FASTA files
import numpy as np                       # Library for numerical computing

# Import functions from other files
# TODO: fix typo! (now commented out as I don't have vcf)
from frequency import get_freq, calculate_triplet_counts, process_triplet_positions, compute_triplet_counts, calculate_probabilities 
from output_paths import get_output_folders, get_output_path, write_output
from randomized_operations import random_sampling, get_random_position_in_gene, get_transcript_position 
from ensembl_request import get_hgvs_genomic, hgvs_converter 
from vcf_output import vcf_writer

# TODO: Call in main argument for rerun failed requests - YES
# TODO: test input/output full step
# TODO: Test if request function cancels the program if internet connection is lost 
# TODO: Change naming convention of files/directories?

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
    parser.add_argument('-o', required=False, help="HGVS coding list. Overrides simulation and generate VCF from HGVSC list")          # Adds -o HGVS coding file for VCF generation
    parser.add_argument('-b', type=int, required=False, default=150, help="Batch size for API requests")     
    
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

# 
def read_transcript_list(transcript_list_file):
    """
    Reads in list of transcript IDs to consider
    """
    seqIDs = []
    with open(transcript_list_file, "r") as list_in:
        lines = list_in.readlines()
    seqIDs = [(l.strip().split('|')[0]).split('.')[0] for l in lines] # processing to ignore .version and other info
    return seqIDs                                                     # example output: ['ENST000001', 'ENST000002']


# Opens a txt file containing HGVS coding sequences and returns a list
def open_hgvsc(args):
    """
    Reads a file of HGVS coding sequences and returns a list
    """
    with open(args["o"], "r") as file:
        lines = [line.strip() for line in file.readlines()]
    
    return lines

#endregion


# ====================
# MAIN PROGRAM
# ====================

#region MainProgram

if __name__ == "__main__":                                                  # Checks if the script is executed as the main program or imported as module. If the script is executed as the main program, the code block is executed, not if it is imported as a module.
    
    print("\nParse arguments")
    
    # Parse command-line arguments
    args = parse_arguments()                                                # Calls the function to parse command-line arguments and return the arguments as a dictionary
    db_folder = args['d']                                                   # Assign database folder to a variable    
    # TODO: Check if db_folder exist, otherwise, create it.
    
    #### Override block ####
    
    if args['o'] is not None:     # if a list of HGVS coding is provided
        print("Simulation override. HGVSC list provided")
        
        hgvsc_list = open_hgvsc(args)  
        # TODO: Move these inside function instead?              
        url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        
        print("Access ensemble API")
        hgvs_genomic, hgvs_failed = get_hgvs_genomic(hgvsc_list, url, headers)        # Calls the function to get the HGVS genomic notation from the REST API. The function returns 1 dictionary, 1 list: hgvs_genomic and hgvs_failed.
        
        print("Extract data for VCF")
        chr_info = hgvs_converter(hgvs_genomic)
        
        print("\nWrite VCF:\n")
        ##### Write VCF #####
        # TODO: fix the VCF output path
        # TODO: Make output naming based on input file name?
        
        vcf_output_path = args['o'].rsplit('.', 1)[0] + '.vcf'
        print(f"This is '-o' input: {args['o']}")
        print(f"This is the changed string: {vcf_output_path}")
        vcf_writer(chr_info, vcf_output_path)
    
    
    ########## If no HGVSC list is provided, the program will continue with the simulation ##########
    else:                                                              
        print("No HGVSC list provided")
        print("Read input files and calculate frequencies")
        
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
        
        print("Perform random sampling")
        
        # Perform random sampling based on triplet frequencies
        sampled_triplets = random_sampling(freq, args['n'])      # takes -n mutations as input. Returns a list of sampled elements.
        
        print("Create output directories")
        
        # Creating output directory
        directories = get_output_folders(run_name, args)             # Returns the intended output directory path as a string.
        os.makedirs(directories, exist_ok=True)                # Creates the output directories if they do not exist already.
        
        
        #### Main simulation part ####: 
        
        print("Start simulation")
        
        # Iterating over sampled triplets and looping runs
        for i in range(args['r']):                            # Run the simulation -r times. Numbers series created with range(). The loop will run the number of times specified by the -r argument.
            
            print(f"\nRun {i+1}\n")                           # Print the run number
            print("Simulate mutations")
            
            hgvsc_list = []                                   # Create an empty list to store the output strings
            
            for run in sampled_triplets:                      # Each element in this list is a string representing a triplet and a substitution separated by an underscore "_". The loop iterates over these strings. 
                trip, sub = run.split("_")                    # split() splits the string  at the underscore. The parts are unpacked into the variables trip (triplet) and sub (substitution). If run was a string like "CGA_G/A", then trip would become "CGA" and sub would become "G/A".
                
                # probabilities is a dictionary containing information about each triplet, including names and probabilities.
                # names = probabilities[trip]['name']         # names is a list of transcripts
                # probs = probabilities[trip]['prob']         # probs is a list of probabilities
                names = triplet_counts[trip][0]
                counts = triplet_counts[trip][1]              # counts are used as weights to not have to compute probabilities (faster)
                
                # Randomly select a transcript based on probabilities
                # selected_transcript = np.random.choice(names, p=probs)       # np.random.choice() returns a random transcript based on the probabilities. p=probs is the list of probabilities
                selected_transcript = random.choices(names, weights=counts)[0] # using random.choices (since faster), returns a randomly selected one weighted by the number of times each transcript has the triplet
                
                # Get a random position in the selected gene
                # position = get_random_position_in_gene(pos_in_gene, selected_transcript, trip)   # Inside the function, it first retrieves the list of positions where the given triplet occurs in the given transcript from the posingene dictionary. Then it uses the random.choice() function to select a random position from this list.
                position = get_transcript_position(selected_transcript, trip, db_folder)           # Reads in the corresponding transcript dictionary file, and chooses a random position where that triplet can be found
                
                # Create HGVS coding string
                ch = sub.split("/")                          # split() splits the string at the slash and saves into a list. The parts are unpacked into the variables ch[0] and ch[1]. If sub was a string like "A/B", then ch would become the list ['A', 'B'].
                output_string = f"{selected_transcript}:c.{position + 1}{ch[0]}>{ch[1]}"
                hgvsc_list.append(output_string)             # Append the output string to the list hgvsc_list
                
                # Output HGVS coding to file
                write_output(output_string,directories, run_name, args, i)
            
            
            #### Retrieve chromosome and chromosome position from ENSEMBL REST API ####
            # TODO: Clean up all prints
            
            print("Access ensemble API")
            
            url = "https://rest.ensembl.org/variant_recoder/homo_sapiens"
            headers = {"Content-Type": "application/json", "Accept": "application/json"}
            
            hgvs_genomic, hgvs_failed = get_hgvs_genomic(hgvsc_list, url, headers, args['b'])        # Calls the function to get the HGVS genomic notation from the REST API. The function returns 1 dictionary, 1 list: hgvs_genomic and hgvs_failed.
            
            # Print matching/failed HGVS coding
            print(f"Number of matches: {len(hgvs_genomic)}")
            print(f"Number of failed matches: {len(hgvs_failed)}")
            
            
            #### Create VCF ####
            # TODO: Change name of simulator in vcf_output once decided.
            
            print("Create VCF output")
            chr_info = hgvs_converter(hgvs_genomic)
            
            vcf_output_path = get_output_path(directories, i+1, run_name, args)  # directories - Path to output folder, run_name - Processed name string based on input.
            vcf_output_path = vcf_output_path.rsplit('.', 1)[0] + '.vcf'        # 1 specifies the number of times split() will occur.
            
            vcf_writer(chr_info, vcf_output_path)
            
            print("\nEnd of program")

#endregion


# Fix:
# TODO: Move URL and HEADER inside function instead?
# TODO: Change required to optional?
# TODO: Simplify the naming of files?
# TODO: Retry loop for failed request?
# TODO: break out frequency functions

# Testing:
# TODO: batch size 150, 200, 250 x10-20 runs
# TODO: Run multiple 1000 sim and measure time (function for this)


# argumenr batchsize

# TODO: Fix VCF problem!
# If several mutations are in the same gene, the VCF file will only contain the first repeated several times.
# Fix by removing after creating if multiple are in the same list?