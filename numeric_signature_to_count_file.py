"""Converts SBS signatures from COSMIC into a count file"""

def convert_sbs_to_count_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    header = lines[0].split()
    grch38_index = header.index("SBS1_GRCh38")  # find the index of the desired column

    with open(output_file, 'w') as outfile:
        for line in lines[1:]:  # skip the header line
            parts = line.split()
            mutation = parts[0]
            count_value = float(parts[grch38_index])  # read the floating-point number from SBS1_GRCh38
            count = round(count_value * 1000000)  # Example of converting the value to an integer count (scaling factor may vary)

            # Process the mutation type from something like A[C>A]A to "ACA C/A"
            context = mutation[0] + mutation[6] + mutation[7]
            change = mutation[2:5].replace(">", "/")
            new_line = f"{context}\t{change}\t{count}\n"

            outfile.write(new_line)

# Call the function with the path to your input text file and the output file you want to create
convert_sbs_to_count_file('input.txt', 'output.count')
