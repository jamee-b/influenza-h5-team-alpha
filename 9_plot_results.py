import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def main():
    csv_input_file = input("Enter CSV results file: ")
    results = read_csv(csv_input_file)
    avian_results = results.loc[(results['host'] == 'bovine')]
    plot_data(avian_results)

def plot_data(dataframe):
    sns.boxplot(dataframe, x='subtype', y='cai', hue='host')
    plt.xlabel('Influenza A H5 subtypes')
    plt.ylabel('Codon adaptation index (CAI)')
    plt.show()

def read_csv(csv_input_file):
    path = 'Results/{0}'.format(csv_input_file)
    df = pd.read_csv(path)
    return df

if __name__ == "__main__":
    main()