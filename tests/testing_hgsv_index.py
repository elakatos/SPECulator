"""Testing output indexing"""

# Set the details for the search
fasta = 'sample.fasta'
transcript_id = 'ENST00000222122'  
positions = [683, 684, 685, 686, 687]                 

# Open the FASTA file
with open(fasta, 'r') as fasta_file:
    # Variables to hold the current sequence label and sequence data
    sequence_label = None
    sequence_data = []

    # Flag to indicate when the correct transcript is found
    found_transcript = False

    # Read through the file line by line
    for line in fasta_file:
        line = line.strip()  # Remove whitespace
        if not line:
            continue  # Skip empty lines
        if line.startswith(">"):  # Check if it's a header line
            if sequence_label == transcript_id:  # Check if the previous sequence was the right one
                break  # Exit the loop since we've found and processed the desired transcript
            sequence_label = line[1:]  # Remove the '>' and keep the label
            sequence_data = []  # Reset the sequence data for the new header
            found_transcript = (sequence_label == transcript_id)
        elif found_transcript:
            sequence_data.append(line)  # Add the line of sequence only if we are in the right transcript

    # Combine the collected sequence lines and extract the desired nucleotides
    if found_transcript and sequence_data:
        full_sequence = ''.join(sequence_data)
        extracted_nucleotides = {pos: full_sequence[pos - 1] for pos in positions if pos <= len(full_sequence)}
        print(f"Nucleotides at positions {positions} in transcript {transcript_id} are {extracted_nucleotides}")
    else:
        print(f"Transcript {transcript_id} not found in the file.")
