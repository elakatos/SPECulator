"""
Create VCF output from the processed HGVS.
"""

import vcfpy

# Test data
chromosome = "1"
position = 1000
ref_base = "A"
alt_base = "G"

# Example imput from converter (chrom, pos, ref, alt).
example_input = {
    "ENST00000169551:c.609G>A": ("18", "74158243", "G", "A"),
    "ENST00000519026:c.1396G>A": ('8', '20145849', 'G', 'A'),
    "ENST00000383070:c.217G>A": ('Y', '2787387', 'G', 'A')
    }

# TODO: Add so it reads the input from the converter
# TODO: Fix so it includes all "unique" chromosomes (sorted) into ContigHeaderLine
# TODO: What to do with the allele counts when the data is simulated? Is it even relevant?

# Create a new VCF writer
writer = vcfpy.Writer.from_path(
    # TODO: Change the output file name?
    'output.vcf', 
    vcfpy.Header(
        lines=[
            # TODO: Change "MutationSimulator" to whatever the program will be called
            vcfpy.HeaderLine(key='fileformat', value='VCFv4.2'),                         
            vcfpy.HeaderLine(key='source', value='MutationSimulator'),
            
            # TODO: Add contig lines for every chromosome.
            # TODO: Automate! And sorted.
            vcfpy.ContigHeaderLine.from_mapping({'ID': '18', 'length': 80373285}),
            vcfpy.ContigHeaderLine.from_mapping({'ID': '8', 'length': 145138636}),
            vcfpy.ContigHeaderLine.from_mapping({'ID': 'Y', 'length': 57227415}),
            
            # TODO: Do we need allele count info?
            # AC = Allele count, A = this INFO field can have as many values as there are alternate allelesused to provide counts of each alternate allele found in the genotypes.
            # GT = Describes the genotype of each sample at each site. The genotype is a critical part of variant representation. This part should not be changed.
            vcfpy.InfoHeaderLine.from_mapping({'ID': 'AC', 'Number': 'A', 'Type': 'Integer', 'Description': 'Allele count in genotypes, for autosomal chromosomes assume 1'}),
            vcfpy.FormatHeaderLine(id='GT', number='1', type='String', description='Genotype')
        ],
        # TODO: Change? Add so it includes all samples? Or what is it supposed to contain?
        samples=['Sample1']
    )
)

# Add records to the VCF file
for key, (chrom, pos, ref, alt) in example_input.items():
    record = vcfpy.Record(
        CHROM=chrom,
        POS=int(pos),
        ID='.',
        REF=ref,
        ALT=[vcfpy.Substitution(alt)],
        QUAL='.',
        FILTER=['.'],
        INFO={},
        FORMAT=['GT'],
        # TODO: Change into something neutral?
        # This is where the genotype data for each sample is specified in each record. {'GT': '0/1'} indicates that the sample is heterozygous at this position.
        calls=[vcfpy.Call(sample='Sample1', data={'GT': '0/1'})]  # Assuming heterozygous for the example
    )
    writer.write_record(record)

# Close the VCF writer
writer.close()