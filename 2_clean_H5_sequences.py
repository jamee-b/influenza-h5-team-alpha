import re
from Bio import SeqIO

def main():
    fasta_file = input("Enter fasta file: ")
    output_file = input("Enter output file: ")
    subtype = input("Enter virus subtype (ex. H5N1): ")
    seq_records = read_fasta_to_seq_records(fasta_file)
    keep_records = check_sequence_start(seq_records)
    keep_records = check_sequence_end(seq_records, keep_records)
    keep_records = check_sequence_triplet(seq_records, keep_records)
    keep_records = filter_for_term(seq_records, keep_records, term="Influenza A virus")
    keep_records = filter_for_term(seq_records, keep_records, term=subtype)
    keep_records = remove_term(seq_records, keep_records, term="partial cds")
    kept_seq_records = filter_seq_records(seq_records, keep_records)
    print("There are " + str(len(kept_seq_records)) + " " + subtype + " sequences.")
    write_fasta(kept_seq_records, output_file)

# Parse fasta file and read in as sequence record. Append to list.
def read_fasta_to_seq_records(fasta_file):
    seq_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_records.append(record)
    return seq_records

# Write out sequence records into single fasta file.
def write_fasta(seq_records, output_file):
    with open(output_file, "w") as output_handle:
        for i in seq_records:
            SeqIO.write(i, output_handle, "fasta")

# Generate new list of data that is based on keep_seqs. If true then add to new list.
def filter_seq_records(seq_records, keep_records):
    kept_seq_records = []
    for i in range(0, len(seq_records)):
        if keep_records[i] == True:
            kept_seq_records.append(seq_records[i])
        else:
            pass
    return kept_seq_records

# Filters for term in sequence record description.
def filter_for_term(seq_records, keep_records, term):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if keep_records[i] == True:
            if re.search(term, record.description):
                pass
            else:
                keep_records[i] = False
    return keep_records

# Filters for term in sequence record description.
def remove_term(seq_records, keep_records, term):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if keep_records[i] == True:
            if re.search(term, record.description):
                keep_records[i] = False
            else:
                pass
    return keep_records

# If first three nucleotides are ATG add true to keep_seq list. Else add False.
def check_sequence_start(seq_records):
    keep_records = []
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if record.seq[0:3] == "ATG":
            keep_records.append(True)
        else:
            keep_records.append(False)
    return keep_records

# If last three nucleotides are TAA, TGA, and TAG keep sequences. Else do not.
def check_sequence_end(seq_rcords, keep_records):
    for i in range(0, len(seq_rcords)):
        record = seq_rcords[i]
        if keep_records[i] == True:
            if record.seq[-3:] == "TAA":
                pass
            elif record.seq[-3:] == "TGA":
                pass
            elif record.seq[-3:] == "TAG":
                pass
            else:
                keep_records[i] = False
        else:
            pass
    return keep_records

# Keep sequence if length is divisible by 3.
def check_sequence_triplet(seq_records, keep_records):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if keep_records[i] == True:
            if len(record.seq) % 3 == 0:
                pass
            else:
                keep_records[i] = False
        else:
            pass
    return keep_records

if __name__ == "__main__":
    main()