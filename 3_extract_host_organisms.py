import utils
import pandas as pd
import re

def main():
    for H5 in utils.subtypes:
        read_path = "data/subtype_sequences/{0}/01_{0}.fasta".format(H5)
        write_path = "data/subtype_sequences/{0}/{0}_hosts.csv".format(H5)
        seq_records = utils.read_fasta(read_path)
        host_df = get_hosts(seq_records)
        utils.df_to_csv(host_df, write_path)

def get_hosts(seq_records):
    temp_list = []
    for i in range(0, len(seq_records)):
        record = seq_records[i]

        # Extract host name from sequence description using regular expression.
        try:
            host_name = re.search(r'A/(.*?)/', record.description).group(1)
            temp_list.append(host_name)
        except:
            print("No host extracted for {0}.".format(record.name))

    host_df = pd.DataFrame({'host_organism':list(set(temp_list))})


    # Count occurence of host names extract to get number of sequences for host.
    for i in range(0, host_df['host_organism'].size):
        host_df.loc[i, 'seq_count'] = int(temp_list.count(host_df.loc[i, 'host_organism']))
    
    return host_df


if __name__ == "__main__":
    main()