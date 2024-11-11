"""
Create VCF output from the processed HGVS.
"""

import vcfpy
import re

# TODO in main program:
# TODO: Change output path to run folder

# TODO: Cleanup when done and create copy with examples

# Example input from converter (chrom, pos, ref, alt). Original
#example = {
#    "ENST00000169551:c.609G>A": ("18", "74158243", "G", "A"),
#    "ENST00000519026:c.1396G>A": ('8', '20145849', 'G', 'A'),
#    "ENST00000383070:c.217G>A": ('Y', '2787387', 'G', 'A')
#}

#output = 'output.vcf'



def chromosome_sort_key(chrom):
    return [int(text) if text.isdigit() else text for text in re.split('([0-9]+)', chrom)]


# VCF writer function
def vcf_writer(data, output):
    """
    Takes a dictionary of mutation and genomic information and creates a vcf
    """
    # Dictionary of chromosome lengths
    chromosome_lengths = {
        "1": "248956422",
        "2": "242193529",
        "3": "198295559",
        "4": "190214555",
        "5": "181538259",
        "6": "170805979",
        "7": "159345973",
        "8": "145138636",
        "9": "138394717",
        "10": "133797422",
        "11": "135086622",
        "12": "133275309",
        "13": "114364328",
        "14": "107043718",
        "15": "101991189",
        "16": "90338345",
        "17": "83257441",
        "18": "80373285",
        "19": "58617616",
        "20": "64444167",
        "21": "46709983",
        "22": "50818468",
        "X": "156040895",
        "Y": "57227415",
    }
    
    # Filter data to include only canonical chromosomes
    valid_chromosomes = set([str(i) for i in range(1, 23)] + ['X', 'Y'])
    filtered_data = {key: val for key, val in data.items() if val[0] in valid_chromosomes}
    
    # Prepare header lines for a new VCF writer
    header_lines = [
        vcfpy.HeaderLine(key='fileformat', value='VCFv4.2'),                         
        vcfpy.HeaderLine(key='source', value='SPECulator'),      
        # Automatically add contig lines for every unique chromosome sorted by numerical order and then lexicographical order for 'X', 'Y'.
        *[
            vcfpy.ContigHeaderLine.from_mapping({'ID': chrom, 'length': chromosome_lengths[chrom]}) 
            for chrom in sorted(chromosome_lengths.keys(), key=chromosome_sort_key)
        ],
        
        # 
        vcfpy.InfoHeaderLine.from_mapping({'ID': 'AC', 'Number': 'A', 'Type': 'Integer', 'Description': 'Allele count in genotypes, for autosomal chromosomes assume 1'}),
        vcfpy.FormatHeaderLine.from_mapping({'ID': 'GT', 'Number': '1', 'Type': 'String', 'Description': 'Genotype call'})
    ]
    
    # Create a new VCF writer
    writer = vcfpy.Writer.from_path(
        output,  
        vcfpy.Header(
            lines=header_lines,
            samples=vcfpy.SamplesInfos(['SimulatedSample'])
        )
    )
    
    # Add records to the VCF file
    for key, (chrom, pos, ref, alt) in filtered_data.items():
        # Assuming the type of substitution is SNV for single nucleotide variants
        alt_type = 'SNV' if len(ref) == 1 and len(alt) == 1 else 'MNV'  # Example types, adjust as necessary
        record = vcfpy.Record(
            CHROM=chrom,
            POS=int(pos),
            ID=[key],
            REF=ref,
            ALT=[vcfpy.Substitution(type_=alt_type, value=alt)],
            QUAL='.',
            FILTER=['PASS'],  # Use 'PASS' to indicate the variant passed all filters, or '.' if no filtering was applied
            INFO={},
            FORMAT=['GT'],
            calls=[vcfpy.Call(sample='SimulatedSample', data={'GT': '0/1'})]  # Using heterozygous as placeholder
        )
        writer.write_record(record)
        
    # Close the VCF writer
    writer.close()

# Test program for vcf writer.
if __name__ == "__main__":
    
    mut_info = {
    "ENST00000169551:c.609G>A": ("18", "74158243", "T", "A"),
    "ENST00000519026:c.1396G>A": ('8', '20145849', 'T', 'A'),
    "ENST00000383070:c.217G>A": ('Y', '2787387', 'T', 'A')
    }
    
    filepath = 'output.vcf'
    
    vcf_writer(mut_info, filepath)