import os.path
import pandas as pd
from Bio import SeqIO
from pyteomics import parser

from matplotlib import pyplot as plt
import seaborn as sbn


Cleavage_path = '/Users/christina/Documents/Riboseq/Massspec_Siragusa/Analysis_endothel_data/approach_16_12_24'


###########################################################################################################################################################
# Start of analysis
###########################################################################################################################################################
def main():

    # read in peptide evidence table and select cols of interest
    peptides_found = pd.read_excel(
        '/Users/christina/Documents/Riboseq/Massspec_Siragusa/peptides_processed_endothel.xlsx')
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

    ####### Filtering for SO ###################################################
    # Proteins are listed with ";", separate to obtain lists
    peptides_found['Proteins'] = peptides_found['Proteins'].apply(
        lambda x: str(x).split(';') if type(x) != 'float' else x)

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
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: list(x) if type(x) != list else x)
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot.split('_') if prot.startswith('ENSG') else prot for prot in x])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot[:4] if prot[0].startswith('ENSG') else prot for prot in x])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [':'.join(prot) if prot[0].startswith('ENSG') else prot for prot in x])
    pep_found_groups['Unique_peptide'] = pep_found_groups['Proteins'].apply(
        lambda x: 1 if all([y.startswith('ENSG') for y in x]) else 0)
    unique_peptides = pep_found_groups[pep_found_groups['Unique_peptide'] == 1]
    non_unique_peptides = pep_found_groups[pep_found_groups['Unique_peptide'] == 0]

    # explode protein group information to obtain information per protein
    unique_peptides = unique_peptides.explode(['Proteins'])
    protein_counts = unique_peptides['Proteins'].value_counts()
    proteins_keep = protein_counts.index
    # aggregate sequences per protein and rename
    unique_sequences_per_prot = unique_peptides.groupby(['Proteins'])[
        ['Sequence']].agg(list).reset_index()
    unique_sequences_per_prot = unique_sequences_per_prot.rename(
        columns={'Sequence': 'Unique Sequences'})

    # extract non-unique peptide information
    non_unique_peptides = non_unique_peptides.explode(['Proteins'])
    non_unique_peptides = non_unique_peptides[non_unique_peptides['Proteins'].str.startswith(
        'ENSG')]
    non_unique_sequences_per_prot = non_unique_peptides.groupby(['Proteins'])[
        ['Sequence']].agg(list).reset_index()
    non_unique_sequences_per_prot = non_unique_sequences_per_prot[non_unique_sequences_per_prot['Proteins'].isin(
        unique_sequences_per_prot['Proteins'].to_list())]

    sequences_per_prot = unique_sequences_per_prot.copy()
    sequences_per_prot = sequences_per_prot.set_index('Proteins')
    non_unique_sequences_per_prot = non_unique_sequences_per_prot.set_index(
        'Proteins')
    sequences_per_prot['Non unique sequences'] = non_unique_sequences_per_prot['Sequence']

    # calculate number unique and non-unique peptides per protein
    sequences_per_prot['Nr unique peptides'] = sequences_per_prot['Unique Sequences'].apply(
        lambda x: len(x))
    sequences_per_prot['Nr non-unique peptides'] = sequences_per_prot['Non unique sequences'].apply(
        lambda x: len(x) if type(x) == list else 0)
    sequences_per_prot['total nr peptides'] =

    sequences_per_prot.to_csv(os.path.join(
        Cleavage_path, 'proteins_with_only_SO_peptides_17_12_24.csv'))

    ###### LATER THE CLEAVAGE STEPS, BUT ONLY ONCE I KNOW THE EXACT PARAMETERS#######
    # get possible peptides via cleavage with trypsin as dict
    # cleavage_dict = {}
    # for record in SeqIO.parse('old_SPlitORF_protein_sequences_2023/Split_ORF_protein_sequences_NMD.fasta', 'fasta'):
    #     if record.id in proteins_keep:
    #         cleavage_dict[record.id] = list(parser.cleave(str(record.seq),
    #                                                       parser.expasy_rules['trypsin'], min_length=6, missed_cleavages=2, max_length=35))
    # # make dataframe from possible peptides
    # cleavage_df = pd.DataFrame([(key, seqs, len(seqs)) for (
    #     key, seqs) in cleavage_dict.items()], columns=['Proteins', 'Sequence', 'Nr possible peptides'])
    # # get NMD SO protein names
    # NMD_proteins = cleavage_df['Proteins'].to_list()
    # # filter observed proteins for NMD (exclude RI)
    # sequences_per_prot_NMD = sequences_per_prot[sequences_per_prot['Proteins'].isin(
    #     NMD_proteins)]

    # sequences_per_prot_NMD['Nr possible cleaved peptides'] = 0


if __name__ == "__main__":
    main()
