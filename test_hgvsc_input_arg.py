"""
test file to add argument '-o' 
to the command line and skip parts of the program
"""

import argparse

def parse_arguments():
    """
    Parses command-line arguments and returns them as a dictionary.
    """
    
    parser = argparse.ArgumentParser(description="Script for processing FASTA files and mutational profiles.")   # Creates a new ArgumentParser object with a description
    
    parser.add_argument('-n', type=int, required=True, help="Number of simulated mutations")    # Adds -n number of mutations
    parser.add_argument('-r', type=int, required=True, help="Number of runs")                   # Adds -r number of runs
    parser.add_argument('-o', required=False, help="HGVS coding list. Cancels simulation and generate VCF")          # Adds -o HGVS coding file for VCF generation
    
    args = parser.parse_args()                # Reads the command line arguments 
    return vars(args)                         # Returns the arguments as a dictionary


def open_hgvsc(args):
    """
    Reads a file of HGVS coding sequences and returns a list
    """
    with open(args["o"], "r") as file:
        lines = [line.strip() for line in file.readlines()]
    
    return lines


# Test program
if __name__ == "__main__":
    
    print("Wazzaaaa")
    
    args = parse_arguments()
    
    if args["o"] is not None:
        print(args["o"])
        
        hgvsc_list = open_hgvsc(args)
        for line in hgvsc_list:
            print(line)
        
        print("Skip simulation and generate VCF")
    
    else:
        print("Perform simulation")
        print(args["n"])
        print(args["r"])
    
    print("This shall always be performed")
    
    