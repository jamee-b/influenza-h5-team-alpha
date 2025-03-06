from Bio import SeqIO

subtypes = [
    "H5N1",
    "H5N2",
    "H5N3",
    "H5N4",
    "H5N5",
    "H5N6",
    "H5N7",
    "H5N8"
]

# Read fasta sequences into list. Each element is a biopython seqrecord.
def read_fasta(path):
    seq_records = []
    for record in SeqIO.parse(path, "fasta"):
        seq_records.append(record)
    return seq_records

# Write fasta file by passing list containing seqrecords.
def write_fasta(seq_records, path):
    with open(path, "w") as output_handle:
        for i in seq_records:
            SeqIO.write(i, output_handle, "fasta")