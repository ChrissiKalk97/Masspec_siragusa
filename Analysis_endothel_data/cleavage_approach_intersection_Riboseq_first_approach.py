import pandas as pd
import sys
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
import os.path
import seaborn as sbn

SO_path = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_28.11.2024-10.28.44_NMD_genomic_regions_with_unique_trans_coords'
Cleavage_path = '/Users/christina/Documents/Riboseq/Massspec_Siragusa/Analysis_endothel_data/cleavage_approach/'
Riboseq_path = '/Users/christina/Documents/Riboseq/Riboseq_genomic_own_ServeR_02_12_24/NMD_genome'


def get_same_ORF_coords(orf):
    orf = orf.split(':')
    orf[2] = str(int(orf[2]) + 1)
    orf = ':'.join(orf)
    return orf


protein_regions_first_approach = pd.read_csv(
    '/Users/christina/Documents/Riboseq/Massspec_Siragusa/Analysis_endothel_data/peptides_processed/NMD_proteins_found.csv')

protein_regions_second_approach = pd.read_csv(
    os.path.join(Cleavage_path, 'unique_peptides_per_SO_protein.tsv'), sep='\t', index_col=0)


# plot distributions of unique and total percentage of found peptides
figure, ax = plt.subplots()
fig = sbn.histplot(
    protein_regions_second_approach['Unique peptide percentage found'], ax=ax)
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 500)
fig = fig.get_figure()
fig.savefig(os.path.join(Cleavage_path, 'Unique_peptide_percentage_dist.png'))
plt.close()


figure, ax = plt.subplots()
fig = sbn.histplot(
    protein_regions_second_approach['Proportion found peptides'], ax=ax)
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 500)
fig = fig.get_figure()
fig.savefig(os.path.join(Cleavage_path, 'proportion_found_peptides_dist.png'))
plt.close()


protein_regions_second_approach = protein_regions_second_approach[
    protein_regions_second_approach['Unique peptide percentage found'] > 0.2]
protein_regions_second_approach['Proteins_ribo'] = protein_regions_second_approach['Proteins'].str.split(
    ':')
protein_regions_second_approach['Proteins_ribo'] = protein_regions_second_approach['Proteins_ribo'].apply(
    lambda x: ':'.join(x[:2]))


proteins_first_app = protein_regions_first_approach['Proteins'].to_list()
proteins_found_with_both_approaches = protein_regions_second_approach[protein_regions_second_approach['Proteins'].isin(
    proteins_first_app)]
print('proteins found with both appraoches',
      proteins_found_with_both_approaches)


NMD_ribo_regions = \
    pd.read_csv(os.path.join(SO_path, 'Unique_DNA_Regions_for_riboseq.bed'),
                header=None,
                sep='\t',
                names=['ID', 'start', 'end', 'ORF'])
NMD_ribo_regions['ID'] = NMD_ribo_regions['ID'] + ':' + NMD_ribo_regions['ORF'].str.split(
    ':').apply(lambda x: x[0])


genomic_ribo_files = []
for file in os.listdir(Riboseq_path):
    if file.startswith('OHM') and file.endswith('csv'):
        file_name = os.path.join(Riboseq_path, file)
        genomic_ribo_files.append(file_name)

Ribo_validated_regions = []
for file in genomic_ribo_files:
    Ribo_validation_df = pd.read_csv(file, header=0, index_col=0)
    Ribo_validation_df = Ribo_validation_df[Ribo_validation_df['significant'] == 1]
    Ribo_region_list = Ribo_validation_df['ID'].str.split(
        ':').apply(lambda x: ':'.join(x[:2])).to_list()
    Ribo_validated_regions = Ribo_validated_regions + Ribo_region_list


Riboseq_treatment_NMD_proteomics_intersect = protein_regions_second_approach['Proteins'].isin(
    Ribo_validated_regions)
Riboseq_treatment_NMD_proteomics_intersect.to_csv(
    os.path.join(Cleavage_path, 'NMD_treatment_proteomics_interesct_cleavage_approach.csv'))
venn3([set(Ribo_validated_regions), set(
    protein_regions_second_approach['Proteins_ribo'].to_list()), set(NMD_ribo_regions['ID'].to_list())],
    ('Riboseq treatment NMD', 'Masspec NMD regions', 'NMD_ribo_regions'))
plt.savefig(os.path.join(Cleavage_path,
            'compare_riboseq', 'Venn_Ribo_NMD_masspec_cleavage.png'))
plt.close()


venn2([set(
    protein_regions_second_approach['Proteins'].to_list()),
    set(proteins_first_app)], ('Cleavage_MS_approach', 'First_MS_approach'))
plt.savefig(os.path.join(Cleavage_path,
            'compare_riboseq', 'first_vs_second_MS_approach.png'))
plt.close()
