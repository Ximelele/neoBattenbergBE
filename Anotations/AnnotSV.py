import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from natsort import natsorted
import json

annotsv_data = {}

acmg_class_mapping = {
    '1': 'benign',
    '2': 'likely benign',
    '3': 'variant of unknown significance',
    '4': 'likely pathogenic',
    '5': 'pathogenic',
    'full=1': 'benign',
    'full=2': 'likely benign',
    'full=3': 'variant of unknown significance',
    'full=4': 'likely pathogenic',
    'full=5': 'pathogenic'
}


def chromosome_with_most_pathogenic_mutations(data):
    pathogenic_data = data[data['ACMG_class'].isin(['4', '5', 'full=4', 'full=5'])]

    # Count the number of pathogenic genes per chromosome
    pathogenic_chromosome_counts = pathogenic_data['SV_chrom'].value_counts()

    # Sort chromosomes in natural order (1-22, X, Y, M)
    pathogenic_chromosome_counts_sorted = pathogenic_chromosome_counts.reindex(
        natsorted(pathogenic_chromosome_counts.index))

    dataAnno = {'label': pathogenic_chromosome_counts_sorted.index.tolist(),
                'data': pathogenic_chromosome_counts_sorted.values.tolist()}

    annotsv_data.update(MostPathogenic=dataAnno)


def ACMG_class_descriptive(acmg_class_counts):
    labels = [f"{label} ({count})" for label, count in zip(acmg_class_counts.index, acmg_class_counts.values)]

    plt.figure(figsize=(10, 8))
    colors = sns.color_palette("pastel")[:len(acmg_class_counts)]
    explode = [0.05] * len(acmg_class_counts)  # Slightly separate each slice for emphasis

    plt.pie(acmg_class_counts, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors, explode=explode,
            )
    centre_circle = plt.Circle((0, 0), 0.70, fc='white')
    plt.gca().add_artist(centre_circle)


def top_genes_relation(top_gene_disease_relationships, max_words=3):
    # Shorten the OMIM_phenotype to a specified number of words
    top_gene_disease_relationships['Short_OMIM_phenotype'] = top_gene_disease_relationships['OMIM_phenotype'].apply(
        lambda x: ' '.join(x.split()[:max_words]) + ('...' if len(x.split()) > max_words else '')
    )
    label = top_gene_disease_relationships.apply(lambda x: f"{x['Gene_name']} - {x['Short_OMIM_phenotype']}",
                                                 axis=1).tolist()
    values = top_gene_disease_relationships['Count'].tolist()
    dataAnno = {'label': label, 'data': values}
    annotsv_data.update(TopRelation=dataAnno)


def closest_left(data):
    grouped_data = data.groupby(['SV_type', 'Closest_left']).size().reset_index(name='Count')

    # Sort within each 'SV_type' group and select the top 5 most frequent 'Closest_left' values
    top_closest_left = grouped_data.groupby('SV_type').apply(lambda x: x.nlargest(5, 'Count')).reset_index(drop=True)

    # Prepare stored data without plotting
    stored_data = (
        top_closest_left
        .groupby('SV_type')
        .apply(lambda x: {'label': x['Closest_left'].tolist(), 'data': x['Count'].tolist()})

    )
    # print(stored_data['DEL'])
    # just for deletions
    annotsv_data.update(ClosestLeft=stored_data['DEL'])


def closest_right(data):
    grouped_data = data.groupby(['SV_type', 'Closest_right']).size().reset_index(name='Count')

    # Sort within each 'SV_type' group and select the top 5 most frequent 'Closest_left' values
    top_closest_right = grouped_data.groupby('SV_type').apply(lambda x: x.nlargest(5, 'Count')).reset_index(drop=True)

    stored_data = (
        top_closest_right
        .groupby('SV_type')
        .apply(lambda x: {'label': x['Closest_right'].tolist(), 'data': x['Count'].tolist()})

    )
    annotsv_data.update(ClosestRight=stored_data['DEL'])


def plot_gencc_disease(gencc_counts_sampled):
    dataAnno = {'label': gencc_counts_sampled.index.tolist(),
                'data': gencc_counts_sampled.values.tolist()}
    annotsv_data.update(GenCC=dataAnno)


# def plot_annot_ranking(data):
#     plt.figure(figsize=(8, 5))
#     sns.histplot(data['AnnotSV_ranking_score'].dropna(), bins=20, kde=True, color="darkgreen")
#     plt.title("Distribution of AnnotSV_ranking_score")
#     plt.xlabel("AnnotSV_ranking_score")
#     plt.ylabel("Frequency")
#     plt.show()
#     scores = data['AnnotSV_ranking_score'].dropna()
#
#     # Calculate the histogram data
#     hist_data, bin_edges = np.histogram(scores, bins=20)
#
#     # Prepare the result in the required format
#     result = {
#         'data': hist_data.tolist(),
#         'bins': bin_edges.tolist()
#     }
#     print(result)


def get_acmg_class_data(data):
    localData = data
    localData['ACMG_class'] = localData['ACMG_class'].map(acmg_class_mapping)
    acmg_class_counts = localData['ACMG_class'].value_counts()
    dataAnno = {'label': acmg_class_counts.index.tolist(),
                'data': acmg_class_counts.values.tolist()}
    annotsv_data.update(DonutPlot=dataAnno)


def plot_most_affected(chromosome_counts_sorted):
    chromosome_list = chromosome_counts_sorted.index.tolist()
    chromosome_list = chromosome_list[:-3] + chromosome_list[-2:]

    chromosome_name = chromosome_counts_sorted.values.tolist()
    chromosome_name = chromosome_name[:-3] + chromosome_name[-2:]

    data = {'label': chromosome_list,
            'data': chromosome_name}
    print(data)
    annotsv_data.update(MostAffect=data)


def most_frequent_gene(top_genes_sampled):
    dataAnno = {'label': top_genes_sampled.index.tolist(),
                'data': top_genes_sampled.values.tolist()}
    annotsv_data.update(MostFreq=dataAnno)


def generate_pdf_report(patient_data, filename="Mutation_Report.pdf"):
    data = pd.read_csv(patient_data, sep='\t', low_memory=False)
    data = remove_loc_genes(data)
    chromosome_with_most_pathogenic_mutations(data)

    pathogenic_data = data[data['ACMG_class'].isin(['4', '5', 'full=4', 'full=5'])]
    gene_disease_counts = pathogenic_data.groupby(['Gene_name', 'OMIM_phenotype']).size().reset_index(name='Count')
    top_gene_disease_relationships = gene_disease_counts.nlargest(10, 'Count')
    top_genes_relation(top_gene_disease_relationships)
    top_genes_sampled = data['Gene_name'].value_counts().nlargest(5)
    most_frequent_gene(top_genes_sampled)
    closest_left(data)
    closest_right(data)
    gencc_counts_sampled = data['GenCC_disease'].value_counts().nlargest(5)
    plot_gencc_disease(gencc_counts_sampled)
    # plot_annot_ranking(data)




    chromosome_counts = data['SV_chrom'].value_counts()
    chromosome_counts_sorted = chromosome_counts.reindex(natsorted(chromosome_counts.index))

    plot_most_affected(chromosome_counts_sorted)
    get_acmg_class_data(data)
    print(annotsv_data)
    out_file = open("../annotsv.json", "w")

    json.dump(annotsv_data, out_file, indent=4)
    out_file.close()



def remove_loc_genes(dataframe):
    dataframe = dataframe[~dataframe['Gene_name'].str.contains('LOC', na=False)]
    dataframe.loc[dataframe['Closest_left'].str.contains('LOC', na=False), 'Closest_left'] = np.nan
    dataframe.loc[dataframe['Closest_right'].str.contains('LOC', na=False), 'Closest_right'] = np.nan

    return dataframe


def start_analysis_annotsv(patient_data):
    patient_pdf_name = patient_data.split(".")[0] + "_annotsv" + ".pdf"
    generate_pdf_report(patient_data, patient_pdf_name)


start_analysis_annotsv("/Users/martindruzbacky/PycharmProjects/DP_helper/Lynch_1827_06_N.vcf_annotSV.tsv")
