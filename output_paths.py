"""Functions to create output paths for the main program."""
# TODO: Change name of the file and explanation

# Import modules
from datetime import date                # Library for handling dates
import os                                # Library for interacting with the operating system
import re                                # Library for regular expressions


# ====================
# OUTPUT OPERATIONS
# ====================

#region Output


# Function to write output directories
def get_output_folders(fasta_name, args):
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
def get_output_path(directories_paths, run_index, fasta_name, args):
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

# TODO: remove intermediate

def write_output(output, output_directories, fasta_name, args, i):
    """
    Writes an intermediate file and then appends it to an output file.
    """
    
    # Output paths
    filepath_inter = "output/intermediate.txt"
    filepath_out = get_output_path(output_directories, i+1, fasta_name, args)     # i+1 is the run index
    
    # Creating or overwriting intermediate file
    with open(filepath_inter, "w") as f:
        f.write(output + '\n')
        
    # Appending intermediate file to output file
    with open(filepath_inter, "r") as rf:       # Read mode intermediate file
        with open(filepath_out, "a") as wf:     # Append mode for output file
            for line in rf:
                wf.write(line)

#endregion

