from Bio import SeqIO
import pandas as pd
import re

def main():
    fasta_file = input("Enter fasta file: ")
    subtype = input("Enter subtype: ")
    seq_records = read_fasta_to_seq_records(fasta_file)
    write_host_organisms_csv(seq_records, subtype)

# Parse fasta file and read in as sequence record. Append to list.
def read_fasta_to_seq_records(fasta_file):
    seq_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_records.append(record)
    return seq_records

# Produces dataframe containing unique host organisms with number of sequences per organism.
def write_host_organisms_csv(seq_records, subtype):
    temp_list = []
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        try:
            group_1 = re.search(r'A/(.*?)/', record.description).group(1)
            temp_list.append(group_1)
        except:
            print("No host extracted using regex search.")
    df = pd.DataFrame({'host_organism':list(set(temp_list))})

    for i in range(0, df['host_organism'].size):
        df.loc[i, 'seq_count'] = int(temp_list.count( df.loc[i, 'host_organism']))
    df.to_csv(subtype + "_host_organisms.csv", index=False)

if __name__ == "__main__":
    main()