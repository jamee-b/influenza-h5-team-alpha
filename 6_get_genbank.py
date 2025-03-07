import utils

def main():
    for host in utils.host_names:
        read_path = "data/genbank/host_info/{0}.csv".format(host) # Input file location
        write_path = "data/genbank/{0}.gb".format(host) # Write file destination
        seq_info = utils.csv_to_df(read_path) # Dataframe with accession and protein ids
        handle = utils.get_genbank_handle(accessions=seq_info['accession'].to_list()) # Retrieves genbank handle using list of accession ids
        utils.write_genbank(handle, write_path) # Write genbank file

if __name__ == '__main__':
    main()