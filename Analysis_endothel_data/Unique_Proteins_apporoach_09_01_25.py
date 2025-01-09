import os.path
import pandas as pd
from Bio import SeqIO
from pyteomics import parser

from matplotlib import pyplot as plt
import seaborn as sbn


Cleavage_path = '/Users/christina/Documents/Riboseq/Massspec_Siragusa/Analysis_endothel_data/cleavage_approach/'


###########################################################################################################################################################
# Start of analysis
###########################################################################################################################################################
def main():

    # read in peptide evidence table and select cols of interest
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
    peptides_found = peptides_found[peptides_found['Length'] > 6]
    print('Length filtered', peptides_found.shape)

    ####### Filtering for SO ###################################################
    # Proteins are listed with ";", separate to obtain lists
    peptides_found['Proteins'] = peptides_found['Proteins'].apply(
        lambda x: str(x).split(';') if type(x) != 'float' else x)

    # select only groups of interest, having at least SO protein (starting with ENSG, rest has P plut number)
    peptides_found = peptides_found.loc[peptides_found['Proteins'].apply(lambda x: any(
        [y.startswith('ENSG') for y in x]) if isinstance(x, list) else False), :]

    # Make a copy to be able to use frame later
    pep_found_groups = peptides_found.copy()

    # Filter for Groups that have a unique protein
    pep_found_groups = pep_found_groups.loc[peptides_found['Unique (Groups)']
                                            == 'yes', :]

    ####### Preprocess protein listing and naming for SO #######################
    # keep all entries where the first protein in the group is an SO(ENSG)
    # 1. make all Proteins to a list
    pep_found_groups = pep_found_groups.loc[:, ['Proteins',
                                                'Sequence',
                                                'Leading razor protein',
                                                'Start position',
                                                'End position',
                                                'Length',
                                                'Unique (Groups)',
                                                'Unique (Proteins)',
                                                'Score'
                                                ]]

    # filter for Protein Groups with only one Protein: peptides are unique to that protein
    pep_found_groups['Unique_peptide'] = pep_found_groups['Proteins'].apply(
        lambda x: 1 if len(x) == 1 else 0)
    # Change the SO name from _ to, this throws away the information about all non SO proteins:
    # but this does not falsify anything because we defined the uniqueness above
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot.split('_') for prot in x if prot.startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot[:4] for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [':'.join(prot) for prot in x if prot[0].startswith('ENSG')])

    # divide the df in 2 parts: unique and non-unique
    # count and length of the unique and non-unique peptides
    # require at least two peptides, each of which needs to be at least 7 AA long

    unique_peptides = pep_found_groups[pep_found_groups['Unique_peptide'] == 1]
    # unlist the single protein groups
    unique_peptides['Proteins'] = unique_peptides['Proteins'].apply(
        lambda x: x[0])
    non_unique_peptides = pep_found_groups[pep_found_groups['Unique_peptide'] == 0]

    # assert that all unique peptides are actually labelled as 'yes' for 'Unique (Proteins)'
    assert sum(pep_found_groups['Unique (Proteins)']
               [pep_found_groups['Unique_peptide'] == 1] == 'yes') == unique_peptides.shape[0]

    # get number of unique peptides
    unique_peptides_summarized = unique_peptides.groupby(
        'Proteins')['Unique_peptide'].sum().to_frame()
    # get total length of unique peptides
    unique_peptides_summarized['total_length'] = unique_peptides.groupby('Proteins')[
        'Length'].sum()

    # explode protein group information to obtain information per protein
    non_unique_peptides = non_unique_peptides.explode(['Proteins'])
    non_unique_peptides = non_unique_peptides['Proteins'].value_counts()
    non_unique_peptides = non_unique_peptides[non_unique_peptides.index.isin(
        unique_peptides_summarized.index.to_list())]

    # assign the non-unique peptides to the SO proteins with the unique ones
    unique_peptides_summarized['non_unique_peptides'] = 0
    unique_peptides_summarized.loc[non_unique_peptides.index,
                                   'non_unique_peptides'] = non_unique_peptides

    # count the total nr of peptides
    unique_peptides_summarized['nr_total_peptides'] = unique_peptides_summarized['non_unique_peptides'] + \
        unique_peptides_summarized['Unique_peptide']
    unique_peptides_summarized = unique_peptides_summarized[
        unique_peptides_summarized['nr_total_peptides'] > 1]

    # write the peptides of interest to csv
    unique_peptides_summarized.to_csv(
        'Unique_(Proteins)_approach/unique_peptides_summarized.csv')

    ############################################################################
    # PLOTTING #################################################################
    ############################################################################


if __name__ == "__main__":
    main()
