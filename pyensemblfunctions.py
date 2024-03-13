"""Importing metadata and data from Ensembl database."""

from pyensembl import EnsemblRelease

#TODO argument for ensembl version?
# Specify Ensembl release version
ensembl = EnsemblRelease(104)

# Use ENSEMBL transcript ID to get genomic location
transcript_id = "ENST00000367800"
transcript = ensembl.transcript_by_id(transcript_id)

#print the stuff
print("Chromosome:", transcript.contig)     # .contig is the chromosome
print("Start position:", transcript.start)  # .start is the start position
print("End position:", transcript.end)      # .end is the end position