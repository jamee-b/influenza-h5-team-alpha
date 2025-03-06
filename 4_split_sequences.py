from Bio import SeqIO
from host_keywords import H5N1_avian_hosts, H5N1_bovine_hosts, H5N1_human_hosts
from host_keywords import H5N2_avian_hosts, H5N2_swine_hosts
from host_keywords import H5N3_avian_hosts
from host_keywords import H5N4_avian_hosts
from host_keywords import H5N5_avian_hosts
from host_keywords import H5N6_avian_hosts, H5N6_swine_hosts
from host_keywords import H5N7_avian_hosts
from host_keywords import H5N8_avian_hosts
import re

def main():
    fasta_file = input("Enter FASTA file: ")
    subtype = input("Enter virus subtype: ")
    seq_records = read_fasta_to_seq_records(fasta_file)
    if subtype == "H5N1":
        H5N1_avian_seq_records = split_cds(seq_records, H5N1_avian_hosts)
        H5N1_bovine_seq_records = split_cds(seq_records, H5N1_bovine_hosts)
        H5N1_human_seq_records = split_cds(seq_records, H5N1_human_hosts)
        write_fasta(H5N1_avian_seq_records, output_file="H5N1_chicken_cds.fasta")
        write_fasta(H5N1_bovine_seq_records, output_file="H5N1_bovine_cds.fasta")
        write_fasta(H5N1_human_seq_records, output_file="H5N1_human_cds.fasta")
    elif subtype == "H5N2":
        H5N2_avian_seq_records = split_cds(seq_records, H5N2_avian_hosts)
        H5N2_swine_seq_records = split_cds(seq_records, H5N2_swine_hosts)
        write_fasta(H5N2_avian_seq_records, output_file="H5N2_chicken_cds.fasta")
        write_fasta(H5N2_swine_seq_records, output_file="H5N2_swine_cds.fasta")
    elif subtype == "H5N3":
        H5N3_avian_seq_records = split_cds(seq_records, H5N3_avian_hosts)
        write_fasta(H5N3_avian_seq_records, output_file="H5N3_duck_cds.fasta")
    elif subtype == "H5N4":
        H5N4_avian_seq_records = split_cds(seq_records, H5N4_avian_hosts)
        write_fasta(H5N4_avian_seq_records, output_file="H5N4_duck_cds.fasta")
    elif subtype == "H5N5":
        H5N5_avian_seq_records = split_cds(seq_records, H5N5_avian_hosts)
        write_fasta(H5N5_avian_seq_records, output_file="H5N5_duck_cds.fasta")
    elif subtype == "H5N6":
        H5N6_avian_seq_records = split_cds(seq_records, H5N6_avian_hosts)
        H5N6_swine_seq_records = split_cds(seq_records, H5N6_swine_hosts)
        write_fasta(H5N6_avian_seq_records, output_file="H5N6_duck_cds.fasta")
        write_fasta(H5N6_swine_seq_records, output_file="H5N6_swine_cds.fasta")
    elif subtype == "H5N7":
        H5N7_avian_seq_records = split_cds(seq_records, H5N7_avian_hosts)
        write_fasta(H5N7_avian_seq_records, output_file="H5N7_duck_cds.fasta")
    elif subtype == "H5N8":
        H5N8_avian_seq_records = split_cds(seq_records, H5N8_avian_hosts)
        write_fasta(H5N8_avian_seq_records, output_file="H5N8_chicken_cds.fasta")

# Parse fasta file and read in as sequence record. Append to list.
def read_fasta_to_seq_records(fasta_file):
    seq_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_records.append(record)
    return seq_records

# Puts sequences whose FASTA header description contains host keywords (Such as Chicken, CHICKEN, ... etc)
# into list and returns list of sequences.
def split_cds(seq_records, host_list):
    host_seq_records = []
    for host in host_list:
        for i in range(0, len(seq_records)):
            record = seq_records[i]
            if re.search(host, record.description):
                host_seq_records.append(seq_records[i])
            else:
                pass
    return host_seq_records

# Write out sequence records into single fasta file.
def write_fasta(seq_records, output_file):
    with open(output_file, "w") as output_handle:
        for i in seq_records:
            SeqIO.write(i, output_handle, "fasta")


if __name__ == "__main__":
    main()