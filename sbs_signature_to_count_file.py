"""Converts SBS signatures from COSMIC into a count file."""

def convert_sbs_to_count_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        
    # Assume 'SBS1_GRCh38' is always the third column
    grch38_index = 2  
    
    with open(output_file, 'w') as outfile:
        for line in lines[1:]:                 # Skip the header line
            parts = line.strip().split()
            mutation = parts[0]
            count_value = parts[grch38_index]  
            
            # Extraction and formatting of the mutation context and change
            context = mutation[0] + mutation[2] + mutation[6]  # Extract first base, mutated base, and third base
            change = mutation[2] + '/' + mutation[4]           # Extract REF and ALT bases
            
            new_line = f"{context}\t{change}\t{count_value}\n"
            outfile.write(new_line)

# Convert SBS1 signature to a count file. 
# Put output into 'sbs_signatures' folder
input_file = 'SBS31_reference.txt'
output_file = 'sbs_signatures/sbs31.count'
convert_sbs_to_count_file(input_file, output_file)
