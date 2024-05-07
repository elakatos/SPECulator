"""
Create VCF output from the processed HGVS.
"""

import vcfpy

# Example input from converter (chrom, pos, ref, alt).
example_input = {
    "ENST00000169551:c.609G>A": ("18", "74158243", "G", "A"),
    "ENST00000519026:c.1396G>A": ('8', '20145849', 'G', 'A'),
    "ENST00000383070:c.217G>A": ('Y', '2787387', 'G', 'A')
}

# Prepare header lines for a new VCF writer
header_lines = [
    vcfpy.HeaderLine(key='fileformat', value='VCFv4.2'),                         
    vcfpy.HeaderLine(key='source', value='MutationSimulator'),  # Change "MutationSimulator" to your program's name
    
    # Automatically add contig lines for every unique chromosome sorted by numerical order and then lexicographical order for 'X', 'Y'.
    *[
        vcfpy.ContigHeaderLine.from_mapping({'ID': chrom, 'length': 1})  # Replace '1' with actual length if known
        for chrom in sorted({data[0] for data in example_input.values()}, key=lambda x: (x.isdigit(), x))
    ],
    
    # Allele count information may not be relevant for simulated data unless simulating population data
    vcfpy.InfoHeaderLine.from_mapping({'ID': 'AC', 'Number': 'A', 'Type': 'Integer', 'Description': 'Allele count in genotypes, for autosomal chromosomes assume 1'}),
    vcfpy.FormatHeaderLine.from_mapping({'ID': 'GT', 'Number': '1', 'Type': 'String', 'Description': 'Genotype call'})
]

# Create a new VCF writer
writer = vcfpy.Writer.from_path(
    'output.vcf',  # Consider changing the file name to something more descriptive
    vcfpy.Header(
        lines=header_lines,
        samples=vcfpy.SamplesInfos(['SimulatedSample'])  # Correctly use SamplesInfos class
    )
)

# Add records to the VCF file
# Add records to the VCF file
for key, (chrom, pos, ref, alt) in example_input.items():
    record = vcfpy.Record(
        CHROM=chrom,
        POS=int(pos),
        ID='.',
        REF=ref,
        ALT=[vcfpy.Substitution(value=alt)],  # Ensure this is correctly forming an ALT entry
        QUAL='.',
        FILTER=['PASS'],  # Use 'PASS' to indicate the variant passed all filters, or '.' if no filtering was applied
        INFO={},
        FORMAT=['GT'],
        calls=[vcfpy.Call(sample='Sample1', data={'GT': '0/1'})]  # Assuming heterozygous for the example
    )
    writer.write_record(record)


# Close the VCF writer
writer.close()
