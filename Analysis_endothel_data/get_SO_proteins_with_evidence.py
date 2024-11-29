import csv
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3


###########################################################################################################################################################
# Start of analysis
###########################################################################################################################################################
def main():

    proteinGroups = pd.read_csv('../proteinGroups.txt', sep='\t', header=0)
    print(proteinGroups.head())
    print(proteinGroups.columns)
    proteinGroups = proteinGroups.loc[:, ['Protein IDs',
                                          'Peptide counts (all)',
                                          'Peptide counts (razor+unique)',
                                          'Peptide counts (unique)']]
    print(proteinGroups.head(20))
    proteinGroups = proteinGroups[proteinGroups['Protein IDs'].str.contains(
        'ENSG')]
    print(proteinGroups.head())
    proteinGroups['Peptide counts (unique)'] = proteinGroups['Peptide counts (unique)'].apply(
        lambda x: str(x).split(';') if type(x) != 'int' else x)

    # filter for SO proteins having at least 2 peptides found
    proteinGroups['Peptide counts (all)'] = proteinGroups['Peptide counts (all)'].apply(
        lambda x: str(x).split(';') if type(x) != 'int' else x)
    proteinGroups['Protein IDs'] = proteinGroups['Protein IDs'].apply(
        lambda x: str(x).split(';') if type(x) != 'float' else x)
    print(proteinGroups[((proteinGroups['Protein IDs'].str.len() == 1) & (
        proteinGroups['Peptide counts (all)'].apply(lambda x: int(x[0])) > 1))])

    '''
    peptides_found = pd.read_excel('peptides_processed_endothel.xlsx')
    peptides_found = peptides_found[['Proteins',
                                    'Leading razor protein',
                                     'Gene names',
                                     'Protein names',
                                     'Sequence',
                                     'Start position',
                                     'End position',
                                     'Length',
                                     'Missed cleavages',
                                     'Mass',
                                     'Unique (Groups)',
                                     'Unique (Proteins)',
                                     'Charges',
                                     'PEP',
                                     'Score']]
    # Filter for posterior error probability smaller 10%, is this too lenient?
    peptides_found = peptides_found[peptides_found['PEP'] < 0.1]
    # Proteins are listed with ";", separate to obtain lists
    peptides_found['Proteins'] = peptides_found['Proteins'].apply(
        lambda x: str(x).split(';') if type(x) != 'float' else x)
    # select only groups of interest, having at least SO protein (starting with ENSG, rest has P plut number)
    peptides_found = peptides_found.loc[peptides_found['Proteins'].apply(lambda x: any(
        [y.startswith('ENSG') for y in x]) if isinstance(x, list) else False), :]
    # print(peptides_found.shape)
    # Make a copy to be able to use frame later
    pep_found_groups = peptides_found.copy()
    # Filter for Groups that have a unique protein
    peptides_found = peptides_found.loc[peptides_found['Unique (Groups)']
                                        == 'yes', :]

    # check how often there is an ORF as its own group, so cannot be a different protein
    pep_of_interest = peptides_found.loc[peptides_found['Proteins'].apply(
        lambda x: len(x) == 1), ['Proteins', 'Start position', 'End position']]
    pep_of_interest['Proteins'] = pep_of_interest['Proteins'].apply(
        lambda x: x[0])
    pep_of_interest['Proteins'] = pep_of_interest['Proteins'].str.split('_')
    pep_of_interest['Proteins'] = pep_of_interest['Proteins'].apply(
        lambda x: x[:4])
    pep_of_interest['Proteins'] = pep_of_interest['Proteins'].apply(
        lambda x: ':'.join(x))

    # save peptides that map somewhere (uniquely) to the SplitORF proteins, but that do not necessarily
    # fall within unique regions
    pep_of_interest['Proteins'].to_csv('unique_proteins_of_interest.csv')
    # print(pep_of_interest.head())'''


# run
if __name__ == "__main__":
    main()
