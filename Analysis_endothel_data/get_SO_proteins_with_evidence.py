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
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: list(x) if type(x) != list else x)
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot.split('_') for prot in x if prot.startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot[:4] for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [':'.join(prot) for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Unique_peptide'] = pep_found_groups['Proteins'].apply(
        lambda x: 1 if len(x) == 1 else 0)

    # explode protein group information to obtain information per protein
    pep_found_groups = pep_found_groups.explode(['Proteins'])
    protein_counts = pep_found_groups['Proteins'].value_counts()
    proteins_keep = protein_counts[protein_counts > 1].index
    pep_found_groups = pep_found_groups[pep_found_groups['Proteins'].isin(
        proteins_keep)]
    sequences_per_prot = pep_found_groups.groupby(['Proteins'])[
        ['Sequence']].agg(list).reset_index()

    # get possible peptides via cleavage with trypsin as dict
    cleavage_dict = {}
    for record in SeqIO.parse('old_SPlitORF_protein_sequences_2023/Split_ORF_protein_sequences_NMD.fasta', 'fasta'):
        if record.id in proteins_keep:
            cleavage_dict[record.id] = list(parser.cleave(str(record.seq),
                                                          parser.expasy_rules['trypsin'], min_length=6, missed_cleavages=2, max_length=35))
    # make dataframe from possible peptides
    cleavage_df = pd.DataFrame([(key, seqs, len(seqs)) for (
        key, seqs) in cleavage_dict.items()], columns=['Proteins', 'Sequence', 'Nr possible peptides'])
    # get NMD SO protein names
    NMD_proteins = cleavage_df['Proteins'].to_list()
    # filter observed proteins for NMD (exclude RI)
    sequences_per_prot_NMD = sequences_per_prot[sequences_per_prot['Proteins'].isin(
        NMD_proteins)]

    # get all possible peptides from Uniprot to check the uniqueness of the peptides in SO proteins observed
    peptide_list = []
    for record in SeqIO.parse('UP000005640_9606.fasta', 'fasta'):
        peptide_list = peptide_list + \
            list(parser.cleave(str(record.seq),
                 parser.expasy_rules['trypsin'], min_length=6, missed_cleavages=2, max_length=35))

    peptide_list = set(peptide_list)

    sequences_per_prot_NMD['Proportion found peptides'] = 0
    sequences_per_prot_NMD['Nr possible cleaved peptides'] = 0
    sequences_per_prot_NMD['Nr found peptides'] = 0
    sequences_per_prot_NMD['Unique peptides found'] = 0
    sequences_per_prot_NMD['Unique peptides possible'] = 0
    sequences_per_prot_NMD.loc['Unique peptide percentage found'] = 0

    for row in sequences_per_prot_NMD.index:
        go_on = True
        protein = sequences_per_prot_NMD.loc[row, 'Proteins']
        sequences = sequences_per_prot_NMD.loc[row, 'Sequence']
        sequences_cleaved = cleavage_df[cleavage_df['Proteins']
                                        == protein]['Sequence']
        if not sequences_cleaved.empty:
            sequences_cleaved = sequences_cleaved.iloc[0]
            if not all(any(seq in sc for sc in sequences_cleaved) for seq in sequences):
                if len([seq for seq in sequences if seq.startswith('M')]) > 0:
                    sequences = [seq[1:]
                                 for seq in sequences if seq.startswith('M')]
                    if not all(any(seq in sc for sc in sequences_cleaved) for seq in sequences):
                        print(sequences)
                        print(sequences_cleaved)
                        go_on = False

            if go_on == True:
                if len(sequences_cleaved) > 0:
                    sequences_per_prot_NMD.loc[row, 'Proportion found peptides'] = len(
                        [seq for seq in sequences if seq in sequences_cleaved])/len(sequences_cleaved)
                else:
                    sequences_per_prot_NMD.loc[row,
                                               'Proportion found peptides'] = 0

                sequences_per_prot_NMD.loc[row, 'Nr found peptides'] = len(
                    [seq for seq in sequences if seq in sequences_cleaved])

                sequences_per_prot_NMD.loc[row, 'Nr possible cleaved peptides'] = len(
                    sequences_cleaved)

                sequences = set(sequences)
                unique_seqs_found = {
                    seq for seq in sequences if seq not in peptide_list}
                unique_seqs_possible = {
                    seq for seq in sequences_cleaved if seq not in peptide_list}
                sequences_per_prot_NMD.loc[row, 'Unique peptides possible'] = len(
                    unique_seqs_possible)
                sequences_per_prot_NMD.loc[row,
                                           'Unique peptides found'] = len(unique_seqs_found)
                if len(unique_seqs_possible) > 0:
                    sequences_per_prot_NMD.loc[row, 'Unique peptide percentage found'] = len(
                        unique_seqs_found) / len(unique_seqs_possible)
                else:
                    sequences_per_prot_NMD.loc[row,
                                               'Unique peptide percentage found'] = 0

    print(sequences_per_prot_NMD.head(50))
    print(sequences_per_prot_NMD.tail(50))

    sequences_per_prot_NMD.to_csv(
        'unique_peptides_per_SO_protein.tsv', sep='\t')


    # run
if __name__ == "__main__":
    main()
