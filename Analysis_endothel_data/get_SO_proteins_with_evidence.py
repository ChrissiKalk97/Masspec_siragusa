import csv
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
from Bio import SeqIO
from pyteomics import parser


###########################################################################################################################################################
# Start of analysis
###########################################################################################################################################################
def main():
    # proteinGroups = pd.read_csv('../proteinGroups.txt', sep='\t', header=0)
    # print(proteinGroups.head())
    # print(proteinGroups.columns)
    # proteinGroups = proteinGroups.loc[:, ['Protein IDs',
    #                                       'Peptide counts (all)',
    #                                       'Peptide counts (razor+unique)',
    #                                       'Peptide counts (unique)']]
    # print(proteinGroups.head(20))
    # proteinGroups = proteinGroups[proteinGroups['Protein IDs'].str.contains(
    #     'ENSG')]
    # print(proteinGroups.head())
    # proteinGroups['Peptide counts (unique)'] = proteinGroups['Peptide counts (unique)'].apply(
    #     lambda x: str(x).split(';') if type(x) != 'int' else [x])

    # # filter for SO proteins having at least 2 peptides found
    # proteinGroups['Peptide counts (all)'] = proteinGroups['Peptide counts (all)'].apply(
    #     lambda x: str(x).split(';') if type(x) != 'int' else [x])
    # proteinGroups['Protein IDs'] = proteinGroups['Protein IDs'].apply(
    #     lambda x: str(x).split(';') if type(x) != 'float' else [x])
    # print(proteinGroups[((proteinGroups['Protein IDs'].str.len() == 1) & (
    #     proteinGroups['Peptide counts (all)'].apply(lambda x: int(x[0])) > 1))])

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
    # Make a copy to be able to use frame later
    pep_found_groups = peptides_found.copy()
    # Filter for Groups that have a unique protein
    pep_found_groups = pep_found_groups.loc[peptides_found['Unique (Groups)']
                                            == 'yes', :]
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
    # pep_found_groups = pep_found_groups.loc[pep_found_groups['Proteins'].apply(
    #     lambda x: x[0].startswith('ENSG'))]
    # print(pep_found_groups.head())
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot.split('_') for prot in x if prot.startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot[:4] for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [':'.join(prot) for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Unique_peptide'] = pep_found_groups['Proteins'].apply(
        lambda x: 1 if len(x) == 1 else 0)
    # print(pep_found_groups)
    pep_found_groups = pep_found_groups.explode(['Proteins'])
    # pep_found_groups.groupby(['Proteins']).size()
    # pep_found_groups = pep_found_groups[pep_found_groups['Proteins'].isin()]
    protein_counts = pep_found_groups['Proteins'].value_counts()
    proteins_keep = protein_counts[protein_counts > 1].index
    pep_found_groups = pep_found_groups[pep_found_groups['Proteins'].isin(
        proteins_keep)]
    # print(pep_found_groups)
    sequences_per_prot = pep_found_groups.groupby(['Proteins'])[
        ['Sequence']].agg(list).reset_index()
    print(sequences_per_prot)
    # protein_sequences_dict = dict(
    #     zip(sequences_per_prot['Proteins'], sequences_per_prot['Sequence']))

    # print(protein_sequences_dict)

    cleavage_dict = {}
    for record in SeqIO.parse('old_SPlitORF_protein_sequences_2023/Split_ORF_protein_sequences_NMD.fasta', 'fasta'):
        if record.id in proteins_keep:
            cleavage_dict[record.id] = list(parser.cleave(str(record.seq),
                                                          parser.expasy_rules['trypsin'], min_length=6, missed_cleavages=2, max_length=35))

    cleavage_df = pd.DataFrame([(key, seqs, len(seqs)) for (
        key, seqs) in cleavage_dict.items()], columns=['Proteins', 'Sequence', 'Nr possible peptides'])

    NMD_proteins = cleavage_df['Proteins'].to_list()

    sequences_per_prot_NMD = sequences_per_prot[sequences_per_prot['Proteins'].isin(
        NMD_proteins)]
    print(sequences_per_prot_NMD)

    sequences_per_prot_NMD['Proportion found peptides'] = 0
    sequences_per_prot_NMD['Nr possible cleaved peptides'] = 0
    sequences_per_prot_NMD['Nr found peptides'] = 0
    for row in sequences_per_prot_NMD.index:
        protein = sequences_per_prot_NMD.loc[row, 'Proteins']
        sequences = sequences_per_prot_NMD.loc[row, 'Sequence']
        sequences_cleaved = cleavage_df[cleavage_df['Proteins']
                                        == protein]['Sequence']
        sequences_cleaved = sequences_cleaved.iloc[0]
        if not all(any(seq in sc for sc in sequences_cleaved) for seq in sequences):
            print(sequences)
            print(sequences_cleaved)
        else:
            sequences_per_prot_NMD.loc[row, 'Proportion found peptides'] = len(
                [seq for seq in sequences if seq in sequences_cleaved])/len(sequences_cleaved)

            sequences_per_prot_NMD.loc[row, 'Nr found peptides'] = len(
                [seq for seq in sequences if seq in sequences_cleaved])

            sequences_per_prot_NMD.loc[row, 'Nr possible cleaved peptides'] = len(
                sequences_cleaved)
    print(sequences_per_prot_NMD.head(50))
    print(sequences_per_prot_NMD.tail(50))
    # check_subset = [[seq_single[1:] if seq_single.startswith(
    #     'M') else seq_single for seq_single in seq] for seq in sequences_per_prot_NMD['Sequence']]
    # is_subset = [[seq_single if seq_single in sc else "" for seq_single in seq for sc in s_cleaved]
    #              for seq in check_subset for s_cleaved in cleavage_df['Sequence']]

    # print(f"\nAre the found sequences also found in the cleaved sequences?\
    #        Allow for the sequence being part of the cleaved sequence:  {is_subset == check_subset}")

    # print(cleavage_df)

    # run
if __name__ == "__main__":
    main()
