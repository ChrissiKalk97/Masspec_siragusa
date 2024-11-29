import csv
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3


###########################################################################################################################################################
# Start of analysis
###########################################################################################################################################################
def main():
    """proteinGroups = pd.read_csv('../proteinGroups.txt', sep='\t', header=0)
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
        proteinGroups['Peptide counts (all)'].apply(lambda x: int(x[0])) > 1))])"""

    peptides_found = pd.read_excel('../peptides_processed_endothel.xlsx')
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
    print('PEP filtered', peptides_found.shape)
    peptides_found = peptides_found[peptides_found['Score'] > 20]
    print('Score filtered', peptides_found.shape)

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
    pep_found_groups = pep_found_groups.loc[peptides_found['Unique (Groups)']
                                            == 'yes', :]

    # keep all entries where the first protein in the group is an SO (ENSG)
    # 1. make all Proteins to a list
    pep_found_groups = pep_found_groups.loc[:, ['Proteins',
                                                'Leading razor protein',
                                                'Start position',
                                                'End position',
                                                'Length',
                                                'Unique (Groups)',
                                                'Unique (Proteins)',
                                                'Score'
                                                ]]
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: list(x) if type(x) != list else x)
    pep_found_groups = pep_found_groups.loc[pep_found_groups['Proteins'].apply(
        lambda x: x[0].startswith('ENSG'))]
    print(pep_found_groups.head())
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot.split('_') for prot in x if prot.startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot[:4] for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [':'.join(prot) for prot in x if prot[0].startswith('ENSG')])
    print(pep_found_groups.head())


# run
if __name__ == "__main__":
    main()
