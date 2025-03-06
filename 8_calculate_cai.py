from Bio import SeqIO
import codonbias as cb
import pandas as pd
import numpy as np
import re

def main():
    subtypes = ['H5N1', 
                'H5N2', 
                'H5N3', 
                'H5N4', 
                'H5N5', 
                'H5N6', 
                'H5N7', 
                'H5N8']
    
    h5_chicken_seq_records = read_H5_sequences("chicken", subtypes)
    h5_duck_seq_records = read_H5_sequences("duck", subtypes)
    h5_bovine_seq_records = read_H5_sequences("bovine", subtypes)
    h5_swine_seq_records = read_H5_sequences("swine", subtypes)
    h5_human_seq_records = read_H5_sequences("human", subtypes)
    chicken_seq_records = read_host_sequences('chicken')
    duck_seq_records = read_host_sequences('duck')
    swine_seq_records = read_host_sequences('swine')
    bovine_seq_records = read_host_sequences('bovine')
    human_seq_records = read_host_sequences('human')

    chicken_cai_model = build_cai_model(chicken_seq_records)
    duck_cai_model = build_cai_model(duck_seq_records)
    swine_cai_model = build_cai_model(swine_seq_records)
    bovine_cai_model = build_cai_model(bovine_seq_records)
    human_cai_model = build_cai_model(human_seq_records)

    chicken_results = calculate_cai(chicken_cai_model, h5_chicken_seq_records, 'chicken')
    duck_results = calculate_cai(duck_cai_model, h5_duck_seq_records, 'duck')
    swine_results = calculate_cai(swine_cai_model, h5_swine_seq_records, 'swine')
    bovine_results = calculate_cai(bovine_cai_model, h5_bovine_seq_records, 'bovine')
    human_results = calculate_cai(human_cai_model, h5_human_seq_records, 'human')

    all_results = concat_results([chicken_results, 
                                  duck_results, 
                                  swine_results, 
                                  bovine_results, 
                                  human_results])
    
    results_to_csv(all_results)

# Parse fasta file and read in as sequence record. Append to list.
def read_H5_sequences(host, subtypes):
    seq_records = []
    for subtype in subtypes:
        path = "Data/{0}/{0}_{1}_cds.fasta".format(subtype, host)
        try:
            for record in SeqIO.parse(path, "fasta"):
                seq_records.append(record)
        except:
            print('No {0} sequences for {1}'.format(host, subtype))
    return seq_records

# Parse fasta file and read in as sequence record. Append to list.
def read_host_sequences(host):
    seq_records = []
    path = "Data/Hosts/{0}_reference_genes.fasta".format(host)
    for record in SeqIO.parse(path, "fasta"):
        seq_records.append(record)
    return seq_records

def build_cai_model(host_seq_records):
    host_seqs = [str(record.seq) for record in host_seq_records]
    cai_model = cb.scores.CodonAdaptationIndex(ref_seq=host_seqs)
    return cai_model

def calculate_cai(cai_model, seq_records, host):
    subtypes = []
    seqs = [str(record.seq) for record in seq_records]

    for i in range(0, len(seq_records)):
        record = seq_records[i]
        subtypes.append(re.search(r'(H\dN\d)', record.description).group(1))
    
    data = {'description': [str(record.description) for record in seq_records],
            'host': [host] * len(seq_records),
            'subtype':subtypes,
            'cai':cai_model.get_score(seqs)}
    
    results = pd.DataFrame(data)
    print('CAI calculated for H5 sequences in {0}'.format(host))
    return results

def concat_results(results):
    all_results = pd.concat(results, ignore_index=True)
    return all_results

def results_to_csv(results):
    path = 'Results/cai_values.csv'
    results.to_csv(path, index=False)

if __name__ == "__main__":
    main()









