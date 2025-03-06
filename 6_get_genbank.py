from Bio import SeqIO
from Bio import Entrez
import pandas as pd
from config import email

def main():
    filename = input("Enter filename: ")
    host = input("Enter host: ")
    seq_info = read_csv(filename) # Dataframe containing extracted accession ids and gene products.
    handle = get_genbank_files(accessions=seq_info['accession'].to_list())
    write_genbank_file(handle, host)

def read_csv(filename):
    path = 'Data/Genbank/{0}'.format(filename)
    df = pd.read_csv(path)
    return df

def get_genbank_files(accessions):
    Entrez.email = email
    handle = Entrez.efetch(db='nucleotide', id=accessions, rettype='gb', retmode='text')
    return handle

def write_genbank_file(handle, host):
    path = 'Data/Genbank/{0}_genbank.gb'.format(host)
    with open(path, 'w') as f:
        f.write(handle.read())

if __name__ == '__main__':
    main()