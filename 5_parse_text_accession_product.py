from Bio import SeqIO
import re
import pandas as pd

def main():
    filename = input("Enter file: ")
    host = input('Enter host: ')
    seq_info = extract_sequence_info(filename)
    to_csv(seq_info, host)

def extract_sequence_info(file):
    accessions = []
    products = []
    protein_ids = []
    with open('Data/Codon Usage Database/{0}'.format(file)) as f:
        for line in f:
            if line[0] == '>':
                accession_pattern = r'>(\w*)'
                product_pattern = r'product="(.*?)"'
                protein_id_pattern = r'protein_id="(.*?)"'
                try:
                    accession = re.search(accession_pattern, line).group(1)
                    product = re.search(product_pattern, line).group(1)
                    protein_id = re.search(protein_id_pattern, line).group(1)
                    accessions.append(accession)
                    products.append(product)
                    protein_ids.append(protein_id)
                except:
                    print('Did not extract product for {0}'.format(accession))

    seq_info = pd.DataFrame({'accession':accessions, 'product':products, 'protein_id':protein_ids})
    return seq_info

def to_csv(dataframe, host):
    path = 'Data/Genbank/{0}_accession_product.csv'.format(host)
    dataframe.to_csv(path, index=False)
                
if __name__ == '__main__':
    main()
