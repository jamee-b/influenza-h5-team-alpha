from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

def main():
    gb_input_file = input("Enter genbank file: ")
    csv_input_file = input("Enter csv file: ")
    host = input("Enter host: ")
    gb_records = read_genbank(gb_input_file)
    seq_info = read_csv(csv_input_file)
    seq_records = build_fasta(gb_records, seq_info)
    write_fasta(seq_records, host)

    
def read_genbank(filename):
    gb_records = []
    path = 'Data/Genbank/{0}'.format(filename)
    for record in SeqIO.parse(path, "genbank"):
        gb_records.append(record)
    return gb_records

def read_csv(filename):
    path = 'Data/Genbank/{0}'.format(filename)
    df = pd.read_csv(path)
    return df

def build_fasta(gb_records, seq_info):
    seq_records = []

    for i in range(0, len(gb_records)):
        record = gb_records[i]
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    try:
                        protein_id = feature.qualifiers['protein_id'][0]
                    except:
                        print("{0} not added to sequence records. Does not have protein id.".format(record.name))
                    if protein_id in seq_info['protein_id'].values:
                        try:
                            seq = feature.location.extract(record).seq
                            id = record.id
                            description = record.description
                            name = record.name
                            seq_record = SeqRecord(seq=seq, id=id, description=description, name=name)
                            seq_records.append(seq_record)
                        except:
                            print("{0} was not added to sequence records".format(protein_id))
    
    return seq_records

# Write out sequence records into single fasta file.
def write_fasta(seq_records, host):
    path = 'Data/Hosts/{0}_reference_genes.fasta'.format(host)
    with open(path, "w") as output_handle:
        for i in seq_records:
            SeqIO.write(i, output_handle, "fasta")


if __name__ == '__main__':
    main()