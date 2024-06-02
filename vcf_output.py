"""
Create VCF output from the processed HGVS.
"""

import vcfpy

# TODO in main program:
# TODO: Change output path to run folder

# TODO in this file
# TODO: Change name of simulator in vcf_output once decided.

# TODO: Cleanup when done and create copy with examples
# TODO: Readup on VCF
# TODO: Understand code
# TODO: Sort chrom?

# Example input from converter (chrom, pos, ref, alt). Original
#example = {
#    "ENST00000169551:c.609G>A": ("18", "74158243", "G", "A"),
#    "ENST00000519026:c.1396G>A": ('8', '20145849', 'G', 'A'),
#    "ENST00000383070:c.217G>A": ('Y', '2787387', 'G', 'A')
#}

#output = 'output.vcf'

# VCF writer function
def vcf_writer(data, output):
    """
    Takes a dictionary of mutation and genomic information and creates a vcf
    """
    
    # Prepare header lines for a new VCF writer
    header_lines = [
        vcfpy.HeaderLine(key='fileformat', value='VCFv4.2'),                         
        vcfpy.HeaderLine(key='source', value='SpectrumSimulator'),  # Change "MutationSimulator" to your program's name
        # TODO: Change simulator to final name
        
        # Automatically add contig lines for every unique chromosome sorted by numerical order and then lexicographical order for 'X', 'Y'.
        # TODO: Not working properly. Sorts after 1 numbers (1 then 11 etc).
        *[
            vcfpy.ContigHeaderLine.from_mapping({'ID': chrom, 'length': 1})  # Replace '1' with actual length if known
            for chrom in sorted({data[0] for data in data.values()}, key=lambda x: (x.isdigit(), x))
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
            samples=vcfpy.SamplesInfos(['SimulatedSample'])  # Correctly use SamplesInfos class
        )
    )
    
    # Add records to the VCF file
    for key, (chrom, pos, ref, alt) in data.items():
        # Assuming the type of substitution is SNV for single nucleotide variants
        alt_type = 'SNV' if len(ref) == 1 and len(alt) == 1 else 'MNV'  # Example types, adjust as necessary
        record = vcfpy.Record(
            CHROM=chrom,
            POS=int(pos),
            ID='.',
            REF=ref,
            ALT=[vcfpy.Substitution(type_=alt_type, value=alt)],
            QUAL='.',
            FILTER=['PASS'],  # Use 'PASS' to indicate the variant passed all filters, or '.' if no filtering was applied
            INFO={},
            FORMAT=['GT'],
            calls=[vcfpy.Call(sample='SimulatedSample', data={'GT': '0/1'})]  # Assuming heterozygous for the example
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