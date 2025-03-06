from Bio import Entrez
from config import email

def main():
    search_term = input("Enter NCBI nucleotide search phrase: ")
    output_file = input("Enter output FASTA file name: ")
    Entrez.email = email
    results = get_accessions(search_term)
    acc_list = results[0]
    webenv = results[1]
    query_key = results[2]
    get_sequences(acc_list, webenv, query_key, output_file)

def get_accessions(search_term):
    stream = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y", idtype="acc", retmax=1000)
    search_results = Entrez.read(stream)
    stream.close()
    return [search_results["IdList"], search_results["WebEnv"], search_results["QueryKey"]]

def get_sequences(acc_list, webenv, query_key, output_file):
    batch_size = 100
    count = len(acc_list)
    output = open(output_file, "w")
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Going to download record %i to %i" % (start + 1, end))

        stream = Entrez.efetch(
            db="nucleotide",
            rettype="fasta",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc",
        )   
        
        data = stream.read()
        stream.close()
        output.write(data)
    output.close()

if __name__ == "__main__":
    main()
