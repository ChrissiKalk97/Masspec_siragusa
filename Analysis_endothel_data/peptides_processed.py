import csv
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3


# helper functions
def classify_peptide_as_unique_overlapping(protein_unique_region, protein_evidence, selected_peptides, non_unique_peptides):
    if int(protein_unique_region['start']) <= protein_evidence['Start position'] <= int(protein_unique_region['end']) and\
            int(protein_unique_region['start']) <= protein_evidence['End position'] <= int(protein_unique_region['end']):
        selected_peptides.append(protein_evidence)
    # if the start position lies in the unique region, but the stop posiition lies outside
    elif int(protein_unique_region['start']) <= protein_evidence['Start position'] <= int(protein_unique_region['end']):
        # if at least 5 AA overlap with the unique region count as unique, otherwise count as non-unique
        if int(protein_unique_region['end']) - protein_evidence['Start position'] >= 5:
            selected_peptides.append(protein_evidence)
        else:
            non_unique_peptides.append(protein_evidence)
    elif int(protein_unique_region['start']) <= protein_evidence['End position'] <= int(protein_unique_region['end']):
        if protein_evidence['End position'] - int(protein_unique_region['start']) >= 5:
            selected_peptides.append(protein_evidence)
        else:
            non_unique_peptides.append(protein_evidence)
    else:
        non_unique_peptides.append(protein_evidence)
    return selected_peptides, non_unique_peptides


def compare_with_Unique_protein_regions(pep_of_interest, Unique_protein_regions, outname):
    selected_peptides = []
    non_unique_peptides = []
    pep_of_interest_in_unique_region = pep_of_interest[pep_of_interest['Proteins'].isin(
        Unique_protein_regions['ORF'])]
    for index, protein_evidence in pep_of_interest_in_unique_region.iterrows():
        unique_region_match_df = Unique_protein_regions[Unique_protein_regions['ORF'].isin(
            [protein_evidence['Proteins']])]
        if unique_region_match_df.shape[0] == 1:
            selected_peptides, non_unique_peptides = classify_peptide_as_unique_overlapping(
                unique_region_match_df.squeeze(), protein_evidence, selected_peptides, non_unique_peptides)
        else:
            for index, protein_unique_region in unique_region_match_df.iterrows():
                selected_peptides, non_unique_peptides = classify_peptide_as_unique_overlapping(
                    protein_unique_region, protein_evidence, selected_peptides, non_unique_peptides)

    protein_val_df = pd.DataFrame(selected_peptides)
    protein_val_df.drop_duplicates(ignore_index=True, inplace=True)
    protein_val_df['Peptide length'] = protein_val_df['End position'] - \
        protein_val_df['Start position']
    # protein_val_df['nr non-unique peptides'] = 0
    # protein_val_df['length non-unique peptides'] = 0

    print('Number of non-unique peptides found', len(non_unique_peptides))
    non_unique_df = pd.DataFrame(non_unique_peptides)
    non_unique_df['Peptide length'] = non_unique_df['End position'] - \
        non_unique_df['Start position']
    non_unique_df_summary = non_unique_df.groupby('Proteins').count()
    non_unique_df_summary['non-unique peptide length'] = non_unique_df.groupby(
        'Proteins')['Peptide length'].sum()
    non_unique_df_summary['nr non-unique peptides'] = non_unique_df_summary['Start position']

    # protein_val_df.loc[protein_val_df['Proteins'].isin(non_unique_df_summary.index), 'nr non-unique peptides'] =\
    # non_unique_df_summary['Start position']

    plt.hist(protein_val_df['Peptide length'])
    plt.title('Length of peptides in unique regions')
    plt.xlabel('Peptide length')
    plt.ylabel('Frequency')
    plt.savefig(f'{outname}.png')
    plt.close()

    print('Total number of peptides falling into unique regions',
          protein_val_df.shape[0])
    print('Number of unique regions with more than one peptide', sum(
        protein_val_df.groupby(['Proteins']).count()['Start position'] > 1))
    print('Number of distinct unique regions in which at least one peptide falls',
          protein_val_df.groupby(['Proteins']).count().shape[0])
    # print(protein_val_df.groupby(['Proteins'])[['Peptide length']].sum())
    length_summed_df = protein_val_df.groupby(
        ['Proteins'])[['Peptide length']].sum()
    length_summed_df['nr unique peptides'] = protein_val_df.groupby(
        ['Proteins']).count()['Start position']

    length_summed_df = length_summed_df.join(
        non_unique_df_summary[['nr non-unique peptides', 'non-unique peptide length']])
    length_summed_df.fillna(0, inplace=True)

    plt.hist(length_summed_df['Peptide length'])
    plt.title('Length of peptides in unique regions')
    plt.xlabel('Peptide length')
    plt.ylabel('Frequency')
    plt.savefig(f'{outname}_length_summed_per_unique_region.png')
    plt.close()

    return length_summed_df, Unique_protein_regions


def find_non_unique_PSMs(summary_peptide_SO_df, pep_found_groups):
    for protein_name, information in summary_peptide_SO_df.iterrows():
        # check if protein is part of the protein group and if so, then the found peptide is counted
        pep_filtered = pep_found_groups.loc[pep_found_groups['Proteins'].apply(
            lambda x: [prot for prot in x if prot == protein_name] if protein_name in x else False) != False, :].copy()
        pep_filtered['peptide length'] = pep_filtered['End position'] - \
            pep_filtered['Start position']
        summary_peptide_SO_df.loc[protein_name,
                                  'nr non-unique peptides'] += pep_filtered['Start position'].nunique()
        summary_peptide_SO_df.loc[protein_name,
                                  'non-unique peptide length'] += pep_filtered['peptide length'].sum()
    return summary_peptide_SO_df


###########################################################################################################################################################
# Start of analysis
###########################################################################################################################################################
def main():
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
    # print(pep_of_interest.head())

    # Check for start and stop positions whether this overlaps the calculated unique regions
    # groupby NMD and RI
    RI_prot_regions = \
        pd.read_csv('/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_01.11.2024-16.13.04_RI_genomic_coords/Unique_Protein_Regions_gt8_valid_filtered.bed',
                    header=None,
                    sep='\t',
                    names=['ORF', 'start', 'end'])

    NMD_prot_regions = \
        pd.read_csv('/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_01.11.2024-14.55.34_NMD_genomic_coords/Unique_Protein_Regions_gt8_valid_filtered.bed',
                    header=None,
                    sep='\t',
                    names=['ORF', 'start', 'end'])

    # ANALYSIS OF RI AND NMD UNIQUE REGIONS
    summary_peptide_SO_df_RI, RI_prot_regions = compare_with_Unique_protein_regions(
        pep_of_interest, RI_prot_regions, 'RI_peptide_length_dist')
    summary_peptide_SO_df_NMD, NMD_prot_regions = compare_with_Unique_protein_regions(
        pep_of_interest, NMD_prot_regions, 'NMD_peptide_length_dist')

    pep_found_groups = pep_found_groups.loc[pep_found_groups['Proteins'].apply(
        lambda x: len(x) > 1)]
    # this filters all proteins out that do not start with ENSG!
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot.split('_') for prot in x if prot.startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [prot[:4] for prot in x if prot[0].startswith('ENSG')])
    pep_found_groups['Proteins'] = pep_found_groups['Proteins'].apply(
        lambda x: [':'.join(prot) for prot in x if prot[0].startswith('ENSG')])

    print('pep_found_groups', pep_found_groups['Proteins'])
    summary_peptide_SO_df_RI = find_non_unique_PSMs(
        summary_peptide_SO_df_RI, pep_found_groups)
    summary_peptide_SO_df_NMD = find_non_unique_PSMs(
        summary_peptide_SO_df_NMD, pep_found_groups)
    print(summary_peptide_SO_df_RI)
    print(summary_peptide_SO_df_RI.loc[(summary_peptide_SO_df_RI['nr non-unique peptides'] > 0) | (
        summary_peptide_SO_df_RI['Peptide length'] > 7), :])
    print(summary_peptide_SO_df_NMD)
    print(summary_peptide_SO_df_NMD.loc[(summary_peptide_SO_df_NMD['nr non-unique peptides'] > 0) | (
        summary_peptide_SO_df_NMD['Peptide length'] > 7), :])

    summary_peptide_SO_df_NMD = summary_peptide_SO_df_NMD.loc[(
        summary_peptide_SO_df_NMD['nr non-unique peptides'] > 0) | (summary_peptide_SO_df_NMD['Peptide length'] > 7), :]
    summary_peptide_SO_df_RI = summary_peptide_SO_df_RI.loc[(
        summary_peptide_SO_df_RI['nr non-unique peptides'] > 0) | (summary_peptide_SO_df_RI['Peptide length'] > 7), :]

    # Check which ones were also found in Riboseq
    # NMD regions

    def get_same_ORF_coords(orf):
        orf = orf.split(':')
        orf[2] = str(int(orf[2]) + 1)
        orf = ':'.join(orf)
        return orf

    Riboseq_control_NMD_regions = pd.read_csv(
        '/Users/christina/Documents/Riboseq/Riboseq_new_08_10_24/endothel_control_NMD_unique_regions.csv', header=0, index_col=0)
    Riboseq_control_NMD_regions = Riboseq_control_NMD_regions[
        Riboseq_control_NMD_regions['signficant'] == 1]
    Riboseq_control_NMD_regions['ORF'] = Riboseq_control_NMD_regions['ORF'].apply(
        lambda x: get_same_ORF_coords(x))
    Riboseq_control_NMD_regions['Protein'] = Riboseq_control_NMD_regions['ID'] + \
        ':' + Riboseq_control_NMD_regions['ORF']
    Riboseq_control_NMD_proteomics_intersect = Riboseq_control_NMD_regions[Riboseq_control_NMD_regions['Protein'].isin(
        summary_peptide_SO_df_NMD.index)]
    Proteomics_regions_NMD = summary_peptide_SO_df_NMD.index.to_list()

    venn3([set(Riboseq_control_NMD_regions['Protein'].to_list()), set(Proteomics_regions_NMD), set(
        NMD_prot_regions['ORF'].to_list())], ('Riboseq control NMD', 'Masspec NMD regions', 'all protein NMD unique regions'))
    plt.savefig(f'Venn_NMD_control.png')
    plt.close()

    Riboseq_treatment_NMD_regions = pd.read_csv(
        '/Users/christina/Documents/Riboseq/Riboseq_new_08_10_24/endothel_treatment_NMD_unique_regions.csv', header=0, index_col=0)
    Riboseq_treatment_NMD_regions = Riboseq_treatment_NMD_regions[
        Riboseq_treatment_NMD_regions['signficant'] == 1]
    Riboseq_treatment_NMD_regions['ORF'] = Riboseq_treatment_NMD_regions['ORF'].apply(
        lambda x: get_same_ORF_coords(x))
    Riboseq_treatment_NMD_regions['Protein'] = Riboseq_treatment_NMD_regions['ID'] + \
        ':' + Riboseq_treatment_NMD_regions['ORF']
    Riboseq_treatment_NMD_proteomics_intersect = Riboseq_treatment_NMD_regions[Riboseq_treatment_NMD_regions['Protein'].isin(
        summary_peptide_SO_df_NMD.index)]
    Riboseq_treatment_NMD_proteomics_intersect.to_csv(
        'NMD_proteomics_interesct.csv')
    summary_peptide_SO_df_NMD.to_csv('NMD_proteins_found.csv')
    venn3([set(Riboseq_treatment_NMD_regions['Protein'].to_list()), set(Proteomics_regions_NMD), set(
        NMD_prot_regions['ORF'].to_list())], ('Riboseq treatment NMD', 'Masspec NMD regions', 'all protein NMD unique regions'))
    plt.savefig(f'Venn_NMD_treatment.png')
    plt.close()

    # RI regions
    Riboseq_control_RI_regions = pd.read_csv(
        '/Users/christina/Documents/Riboseq/Riboseq_new_08_10_24/endothel_control_RI_unique_regions.csv', header=0, index_col=0)
    Riboseq_control_RI_regions = Riboseq_control_RI_regions[
        Riboseq_control_RI_regions['signficant'] == 1]
    Riboseq_control_RI_regions['ORF'] = Riboseq_control_RI_regions['ORF'].apply(
        lambda x: get_same_ORF_coords(x))
    Riboseq_control_RI_regions['Protein'] = Riboseq_control_RI_regions['ID'] + \
        ':' + Riboseq_control_RI_regions['ORF']
    Riboseq_control_RI_proteomics_intersect = Riboseq_control_RI_regions[Riboseq_control_RI_regions['Protein'].isin(
        summary_peptide_SO_df_RI.index)]
    Proteomics_regions_RI = summary_peptide_SO_df_RI.index.to_list()
    venn3([set(Riboseq_control_RI_regions['Protein'].to_list()), set(Proteomics_regions_RI), set(
        RI_prot_regions['ORF'].to_list())], ('Riboseq control RI', 'Masspec RI regions', 'all protein RI unique regions'))
    plt.savefig(f'Venn_RI_control.png')
    plt.close()

    Riboseq_treatment_RI_regions = pd.read_csv(
        '/Users/christina/Documents/Riboseq/Riboseq_new_08_10_24/endothel_treatment_RI_unique_regions.csv', header=0, index_col=0)
    Riboseq_treatment_RI_regions = Riboseq_treatment_RI_regions[
        Riboseq_treatment_RI_regions['signficant'] == 1]
    Riboseq_treatment_RI_regions['ORF'] = Riboseq_treatment_RI_regions['ORF'].apply(
        lambda x: get_same_ORF_coords(x))
    Riboseq_treatment_RI_regions['Protein'] = Riboseq_treatment_RI_regions['ID'] + \
        ':' + Riboseq_treatment_RI_regions['ORF']
    Riboseq_treatment_RI_proteomics_intersect = Riboseq_treatment_RI_regions[Riboseq_treatment_RI_regions['Protein'].isin(
        summary_peptide_SO_df_RI.index)]
    Riboseq_treatment_RI_proteomics_intersect.to_csv(
        'RI_proteomics_interesct.csv')
    summary_peptide_SO_df_RI.to_csv('RI_proteins_found.csv')
    venn3([set(Riboseq_treatment_RI_regions['Protein'].to_list()), set(Proteomics_regions_RI), set(
        RI_prot_regions['ORF'].to_list())], ('Riboseq treatment RI', 'Masspec RI regions', 'all protein RI unique regions'))
    plt.savefig(f'Venn_RI_treatment.png')
    plt.close()

    # Check which regions can be found with Riboseq at all (most Protein unique regions cannot be found with Riboseq as the
    # DNA unique regions are a lot less)
    #######################################################################################################################
    # NMD
    #######################################################################################################################
    unique_region_intersection_NMD = \
        pd.read_csv('/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_01.11.2024-14.55.34_NMD_genomic_coords/Unique_Regions_Overlap_transcriptomic.bed',
                    header=None,
                    sep='\t',
                    names=['tid', 'start', 'end', 'ORF'])
    unique_region_intersection_NMD['ID'] = unique_region_intersection_NMD['tid'] + \
        ':' + unique_region_intersection_NMD['ORF']

    # Proteomics_regions_NMD = [get_same_ORF_coords(
    #     region) for region in Proteomics_regions_NMD]
    print('Number of NMD regions found with Proteomics that also have a DNA unique region and could be validated with Riboseq:')
    print(len([prot_reg for prot_reg in Proteomics_regions_NMD if prot_reg in unique_region_intersection_NMD['ID'].to_list()]))
    print('total number of protein NMD regions found:',
          len(Proteomics_regions_NMD))

    possible_riboseq_regions_NMD = [
        prot_reg for prot_reg in Proteomics_regions_NMD if prot_reg in unique_region_intersection_NMD['ID'].to_list()]
    with open("possible_riboseq_regions_protein_validated_NMD.csv", mode="w", newline="") as file:
        writer = csv.writer(file)
        for item in possible_riboseq_regions_NMD:
            writer.writerow([item])

    #######################################################################################################################
    # RI
    #######################################################################################################################
    unique_region_intersection_RI = \
        pd.read_csv('/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_01.11.2024-16.13.04_RI_genomic_coords/Unique_Regions_Overlap_transcriptomic.bed',
                    header=None,
                    sep='\t',
                    names=['tid', 'start', 'end', 'ORF'])
    unique_region_intersection_RI['ID'] = unique_region_intersection_RI['tid'] + \
        ':' + unique_region_intersection_RI['ORF']

    # Proteomics_regions_RI = [get_same_ORF_coords(
    #    region) for region in Proteomics_regions_RI]
    print('Number of RI regions found with Proteomics that also have a DNA unique region and could be validated with Riboseq:')
    print(len([prot_reg for prot_reg in Proteomics_regions_RI if prot_reg in unique_region_intersection_RI['ID'].to_list()]))
    possible_riboseq_regions_RI = [
        prot_reg for prot_reg in Proteomics_regions_RI if prot_reg in unique_region_intersection_RI['ID'].to_list()]
    with open("possible_riboseq_regions_protein_validated_RI.csv", mode="w", newline="") as file:
        writer = csv.writer(file)
        for item in possible_riboseq_regions_RI:
            writer.writerow([item])
    print('total number of protein RI regions found:', len(Proteomics_regions_RI))


# run
if __name__ == "__main__":
    main()
