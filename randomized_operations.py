"""Functions to randomize triplets, transcripts, substitutions"""

# Import modules
import numpy as np
import random
import pickle


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

