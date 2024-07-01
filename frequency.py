"""Frequencies and position functions"""

import pickle

# Function to calculate frequencies
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
            
            # TODO: Add separation here
            # Separate count and predefined frequency
            
            count = float(count)                                     # Converts the count string to an integer
            if triplet not in freq:                                # Checking if triplet is not in the dictionary - Might seem unecessary since each triplet should only occur once can prevent errors if the count files are wrong
                freq[triplet] = {}                                 # If not, creates a new dictionary for the triplet as value
                                                        # output example: {'CGA': {}} - The empty dictionary will be "change" 
            freq[triplet][change] = {'count': count}               # Adds new key-value pair to the dictionary
            total_count += count                        # output example: {'CGA': {'G/A': {'count': 10}}} - Triple nested dictionary
            
    # Calculating frequencies
    total_count = int(round(total_count))                        # Checks if the total count is not a whole number
    
    for triplet in freq:                                           # Iterates over the triplets in the dictionary               
        for change in freq[triplet]:          #TODO This is unecessary since there should only be one change per triplet
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
