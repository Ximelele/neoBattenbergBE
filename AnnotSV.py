import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from natsort import natsorted
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet

from StrVCTVRE_annot import save_plot_to_buffer


def chromosome_with_most_pathogenic_mutations(pathogenic_chromosome):
    plt.figure(figsize=(10, 6))
    sns.barplot(x=pathogenic_chromosome.index, y=pathogenic_chromosome.values,
                palette="Reds")
    plt.title("Chromosomes with Most Pathogenic Genes (ACMG classes 4 and 5)")
    plt.xlabel("Chromosome")
    plt.ylabel("Count of Pathogenic Genes")
    plt.xticks(rotation=45)


def ACMG_class_descriptive(acmg_class_counts):
    labels = [f"{label} ({count})" for label, count in zip(acmg_class_counts.index, acmg_class_counts.values)]

    # Plotting the fancy donut chart with exploded segments and rich colors
    plt.figure(figsize=(10, 8))
    colors = sns.color_palette("pastel")[:len(acmg_class_counts)]
    explode = [0.05] * len(acmg_class_counts)  # Slightly separate each slice for emphasis

    # Creating a donut plot by adding a central circle
    plt.pie(acmg_class_counts, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors, explode=explode,
            )
    # Adding a white circle in the center to create the donut effect
    centre_circle = plt.Circle((0, 0), 0.70, fc='white')
    plt.gca().add_artist(centre_circle)
    # Title and display


def top_10_pathogenic(top_pathogenic_genes):
    plt.figure(figsize=(10, 6))
    sns.barplot(y=top_pathogenic_genes.index, x=top_pathogenic_genes.values, palette="Reds")
    plt.title("Top 10 Pathogenic Genes (ACMG classes 4 and 5)")
    plt.xlabel("Count")
    plt.ylabel("Gene")


def top_genes_relation(top_gene_disease_relationships, max_words=3):
    # Shorten the OMIM_phenotype to a specified number of words
    top_gene_disease_relationships['Short_OMIM_phenotype'] = top_gene_disease_relationships['OMIM_phenotype'].apply(
        lambda x: ' '.join(x.split()[:max_words]) + ('...' if len(x.split()) > max_words else '')
    )

    plt.figure(figsize=(12, 8))
    sns.barplot(
        y=top_gene_disease_relationships.apply(lambda x: f"{x['Gene_name']} - {x['Short_OMIM_phenotype']}", axis=1),
        x=top_gene_disease_relationships['Count'], palette="viridis"
    )
    plt.title("Top 10 Pathogenic Gene-Disease Relationships")
    plt.xlabel("Count")
    plt.ylabel("Gene-Disease Pair")
    plt.tight_layout()


def most_frequent_gene(top_genes_sampled):
    plt.figure(figsize=(8, 5))
    sns.barplot(x=top_genes_sampled.values, y=top_genes_sampled.index, palette="magma")
    plt.title("Top 5 Most Frequented Genes")
    plt.xlabel("Count")
    plt.ylabel("Gene")


def closest_left(top_closest_left_sampled):
    plt.figure(figsize=(8, 5))
    sns.barplot(x=top_closest_left_sampled.values, y=top_closest_left_sampled.index, palette="cool")
    plt.title("Top 5 Most Frequented Closest_left")
    plt.xlabel("Count")
    plt.ylabel("Closest_left")


def closest_right(top_closest_right_sampled):
    plt.figure(figsize=(8, 5))
    sns.barplot(x=top_closest_right_sampled.values, y=top_closest_right_sampled.index, palette="Blues")
    plt.title("Top 5 Most Frequented Closest_right")
    plt.xlabel("Count")
    plt.ylabel("Closest_right")


def plot_gencc_disease(gencc_counts_sampled):
    plt.figure(figsize=(8, 5))

    sns.barplot(x=gencc_counts_sampled.values, y=gencc_counts_sampled.index, palette="Purples")
    plt.title("Top 5 GenCC_disease Associations")
    plt.xlabel("Count")
    plt.ylabel("GenCC_disease")


def plot_annot_ranking(data):
    plt.figure(figsize=(8, 5))
    sns.histplot(data['AnnotSV_ranking_score'].dropna(), bins=20, kde=True, color="darkgreen")
    plt.title("Distribution of AnnotSV_ranking_score")
    plt.xlabel("AnnotSV_ranking_score")
    plt.ylabel("Frequency")


def plot_acmg_distribution(acmg_class_counts):
    sns.barplot(x=acmg_class_counts.index, y=acmg_class_counts.values, palette="coolwarm")
    plt.title("Distribution of ACMG_class")
    plt.xlabel("ACMG_class")
    plt.ylabel("Count")
    plt.xticks(rotation=45)


def plot_most_affected(chromosome_counts_sorted):
    plt.figure(figsize=(10, 6))
    sns.barplot(x=chromosome_counts_sorted.index, y=chromosome_counts_sorted.values, palette="viridis")
    plt.title("Most Affected Chromosomes")
    plt.xlabel("Chromosome")
    plt.ylabel("Count")
    plt.xticks(rotation=45)


def generate_pdf_report(patient_data, filename="Mutation_Report.pdf"):
    data = pd.read_csv(patient_data, sep='\t', low_memory=False)
    data = remove_loc_genes(data)
    # Check unique values in ACMG_class to verify that we are filtering correctly
    doc = SimpleDocTemplate(filename, pagesize=letter)
    elements = []

    styles = getSampleStyleSheet()
    elements.append(Paragraph("AnnotSV report analysis", styles['Title']))
    elements.append(Spacer(1, 20))

    # Filter data for pathogenic classes based on observed values in ACMG_class
    pathogenic_data = data[data['ACMG_class'].isin(['4', '5', 'full=4', 'full=5'])]

    # Count the number of pathogenic genes per chromosome
    pathogenic_chromosome_counts = pathogenic_data['SV_chrom'].value_counts()

    # Sort chromosomes in natural order (1-22, X, Y, M)
    pathogenic_chromosome_counts_sorted = pathogenic_chromosome_counts.reindex(
        natsorted(pathogenic_chromosome_counts.index))

    # Plot the number of pathogenic genes per chromosome, sorted naturally

    elements.append(Paragraph("Chromosomes with Most Pathogenic Genes (ACMG classes 4 and 5)", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: chromosome_with_most_pathogenic_mutations(pathogenic_chromosome_counts_sorted),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    # Identify the top 10 most frequent pathogenic genes (based on ACMG classes 4 and 5)
    # Sample data
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

    # Example input data (replace with your actual DataFrame)
    # Map the ACMG classes to descriptive names and calculate counts
    data['ACMG_class_descriptive'] = data['ACMG_class'].map(acmg_class_mapping)
    acmg_class_counts = data['ACMG_class_descriptive'].value_counts()

    elements.append(Paragraph("ACMG_class_descriptive", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: ACMG_class_descriptive(acmg_class_counts),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    top_pathogenic_genes = pathogenic_data['Gene_name'].value_counts().nlargest(10)

    elements.append(Paragraph("ACMG_class_descriptive", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: top_10_pathogenic(top_pathogenic_genes),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())
    # Plot the top 10 pathogenic genes with horizontal bars

    gene_disease_counts = pathogenic_data.groupby(['Gene_name', 'OMIM_phenotype']).size().reset_index(name='Count')
    top_gene_disease_relationships = gene_disease_counts.nlargest(10, 'Count')

    elements.append(Paragraph("Top 10 Gene-Disease Relationships", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: top_genes_relation(top_gene_disease_relationships))
    elements.append(Image(buffer, width=500, height=400))  # Adjust dimensions as needed
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    # Step 2: Plot most frequented genes
    top_genes_sampled = data['Gene_name'].value_counts().nlargest(5)
    elements.append(Paragraph("Top 5 Most Frequented Genes", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: most_frequent_gene(top_genes_sampled),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    # Step 3: Plot most frequented Closest_left
    top_closest_left_sampled = data['Closest_left'].value_counts().nlargest(5)
    elements.append(Paragraph("Top 5 Most Frequented Closest_left", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: closest_left(top_closest_left_sampled),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    # Step 4: Plot most frequented Closest_right
    top_closest_right_sampled = data['Closest_right'].value_counts().nlargest(5)
    top_closest_left_sampled = data['Closest_left'].value_counts().nlargest(5)
    elements.append(Paragraph("Top 5 Most Frequented Closest_right", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: closest_right(top_closest_right_sampled),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    # Step 7: Plot GenCC_disease (top 5 values for simplicity)
    gencc_counts_sampled = data['GenCC_disease'].value_counts().nlargest(5)

    elements.append(Paragraph("Top 5 GenCC_disease Associations", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: plot_gencc_disease(gencc_counts_sampled),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())
    # Step 8: Analyze AnnotSV_ranking_score (using histogram for distribution)
    elements.append(Paragraph("Distribution of AnnotSV_ranking_score", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: plot_annot_ranking(data),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    plt.figure(figsize=(10, 6))
    data['ACMG_class'] = data['ACMG_class'].map(acmg_class_mapping)
    acmg_class_counts = data['ACMG_class'].value_counts()

    elements.append(Paragraph("Distribution of ACMG_class", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: plot_acmg_distribution(acmg_class_counts),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    chromosome_counts = data['SV_chrom'].value_counts()
    chromosome_counts_sorted = chromosome_counts.reindex(natsorted(chromosome_counts.index))

    # Plot the sorted chromosome counts
    elements.append(Paragraph("Most Affected Chromosomes", styles['Heading2']))
    buffer = save_plot_to_buffer(lambda: plot_most_affected(chromosome_counts_sorted),
                                 width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    doc.build(elements)


def remove_loc_genes(dataframe):
    dataframe = dataframe[~dataframe['Gene_name'].str.contains('LOC', na=False)]
    dataframe.loc[dataframe['Closest_left'].str.contains('LOC', na=False), 'Closest_left'] = np.nan
    dataframe.loc[dataframe['Closest_right'].str.contains('LOC', na=False), 'Closest_right'] = np.nan

    return dataframe

def start_analysis_annotsv(patient_data):
    patient_pdf_name = patient_data.split(".")[0] + "_annotsv" + ".pdf"
    # data = pd.read_csv(patient_data, sep='\t', low_memory=False)
    # data = remove_loc_genes(data)
    # print(data['Gene_name'].value_counts())
    # print(data['Closest_left'].value_counts())
    # print(data['Closest_right'].value_counts())
    generate_pdf_report(patient_data, patient_pdf_name)

start_analysis_annotsv("/Users/martindruzbacky/PycharmProjects/DP_helper/Lynch_1827_06_N.vcf_annotSV.tsv")