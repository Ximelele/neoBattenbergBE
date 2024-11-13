from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import io


def save_plot_to_buffer(plt_func, width=6, height=6):
    buffer = io.BytesIO()
    plt_func()
    plt.gcf().set_size_inches(width, height)
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    return buffer


def filter_numeric_info(dataframe):
    filtered_infoData = dataframe['INFO'].apply(lambda x: x.split(';')[-1] if pd.notna(x) else None).apply(
        lambda x: x.split('=')[-1] if pd.notna(x) else None)
    return pd.to_numeric(filtered_infoData, errors='coerce').dropna()


# Polar Plot Function
def plot_mutation_polar_plot(dataframe):
    allowed_chromosomes = ["chr" + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df = dataframe[dataframe['#CHROM'].isin(allowed_chromosomes)]
    filtered_infoData = filter_numeric_info(filtered_df)
    chromosome_names = sorted(filtered_df['#CHROM'].unique(),
                              key=lambda x: (int(x[3:]) if x[3:].isdigit() else float('inf'), x))
    theta = np.linspace(0.0, 2 * np.pi, len(chromosome_names), endpoint=False)
    radii = [filtered_infoData[filtered_df['#CHROM'] == chrom].count() for chrom in chromosome_names]
    radii = np.array(radii) / max(radii) * 10 if radii and max(radii) > 0 else radii
    width = 2 * np.pi / len(chromosome_names)
    plt.figure()
    ax = plt.subplot(projection='polar')
    for i, (angle, radius) in enumerate(zip(theta, radii)):
        color = 'black' if radius == min(radii) else 'blue' if radius == max(radii) else plt.cm.spring(radius / 10.)
        ax.bar(angle, radius, width=width, bottom=0, color=color, edgecolor='black')
    ax.set_xticks(theta)
    ax.set_xticklabels(chromosome_names, fontsize=8, fontweight='bold')
    plt.title("Exonic Mutation Distribution Across Chromosomes")


# Donut Plot for Exonic vs Non-Exonic Mutation Distribution
def plot_mutation_donut_chart(dataframe):
    filtered_infoData = dataframe['INFO'].apply(lambda x: x.split(';')[-1] if pd.notna(x) else None).apply(
        lambda x: x.split('=')[-1] if pd.notna(x) else None)
    non_exonic_count = filtered_infoData.value_counts()['not_exonic']
    exonic_info_data = pd.to_numeric(filtered_infoData, errors='coerce').dropna()
    exonic_count = len(exonic_info_data)
    counts = pd.DataFrame({'Mutation Type': ['Non-exonic', 'Exonic'], 'Count': [non_exonic_count, exonic_count]})
    data = counts['Count'].values
    labels = counts['Mutation Type'].values
    colors = ['#ff9999', '#66b3ff']
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw=dict(aspect="equal"))
    wedges, texts = ax.pie(data, wedgeprops=dict(width=0.3, edgecolor='w'), startangle=-40, colors=colors)
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y, x = np.sin(np.deg2rad(ang)), np.cos(np.deg2rad(ang))
        alignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        ax.annotate(f"{labels[i]}: {data[i]}", xy=(x, y), xytext=(1.4 * np.sign(x), 1.4 * y),
                    horizontalalignment=alignment, arrowprops=dict(arrowstyle="-"))


def filter_info(dataframe):
    filtered_infoData = dataframe['INFO'].apply(lambda x: x.split(';')[-1] if pd.notna(x) else None).apply(
        lambda x: x.split('=')[-1] if pd.notna(x) else None)
    exonic_df = dataframe.copy()
    exonic_df['INFO'] = pd.to_numeric(filtered_infoData, errors='coerce').dropna()
    exonic_df = exonic_df.dropna(subset=['INFO'])

    # Sort chromosomes
    sorted_chromosomes = sorted(exonic_df['#CHROM'].unique(),
                                key=lambda x: (int(x[3:]) if x[3:].isdigit() else float('inf'), x))
    exonic_df['#CHROM'] = pd.Categorical(exonic_df['#CHROM'], categories=sorted_chromosomes, ordered=True)

    return sorted_chromosomes, exonic_df


def box_plot_stats(dataframe):
    sorted_chromosomes, exonic_df = filter_info(dataframe)

    # Initialize a dictionary to hold key statistics
    chromosome_stats = {}

    # Calculate statistics for each chromosome
    for chrom in sorted_chromosomes:
        chrom_data = exonic_df[exonic_df['#CHROM'] == chrom]['INFO']
        if not chrom_data.empty:
            median = chrom_data.median()
            q1 = chrom_data.quantile(0.25)
            q3 = chrom_data.quantile(0.75)
            iqr = q3 - q1
            outliers = chrom_data[(chrom_data < q1 - 1.5 * iqr) | (chrom_data > q3 + 1.5 * iqr)].count()
            chromosome_stats[chrom] = {
                'median': median,
                'iqr': iqr,
                'outliers': outliers
            }
    return generate_dynamic_observations(chromosome_stats)


def generate_dynamic_observations(chromosome_stats):
    highest_median_chrom = max(chromosome_stats, key=lambda x: chromosome_stats[x]['median'])
    lowest_median_chrom = min(chromosome_stats, key=lambda x: chromosome_stats[x]['median'])

    # Find chromosomes with the highest variability (top 3 by IQR)
    high_variability_chromosomes = sorted(
        chromosome_stats.items(), key=lambda x: x[1]['iqr'], reverse=True)[:3]

    # Find chromosomes with the most outliers (top 3 by outlier count)
    high_outlier_chromosomes = sorted(
        chromosome_stats.items(), key=lambda x: x[1]['outliers'], reverse=True)[:3]

    # Generate observations
    observations = [
        f"<b>{highest_median_chrom}</b> has the highest median mutation count of {chromosome_stats[highest_median_chrom]['median']:.2f}.",
        f"<b>{lowest_median_chrom}</b> has the lowest median mutation count of {chromosome_stats[lowest_median_chrom]['median']:.2f}."]

    # Add highest and lowest median observations

    # Sort and add high variability chromosomes
    high_variability_chromosomes = sorted(high_variability_chromosomes, key=lambda x: x[1]['iqr'], reverse=True)
    for chrom, stats in high_variability_chromosomes:
        observations.append(
            f"<b>{chrom}</b> shows high variability, with an interquartile range (IQR) of {stats['iqr']:.2f}.")

    high_outlier_chromosomes = sorted(high_outlier_chromosomes, key=lambda x: x[1]['outliers'], reverse=True)
    for chrom, stats in high_outlier_chromosomes:
        if stats['outliers'] > 0:
            observations.append(f"<b>{chrom}</b> has a significant number of outliers ({stats['outliers']}).")

    return "<br />".join(observations)


def plot_mutation_box_plot(dataframe):
    sorted_chromosomes, exonic_df = filter_info(dataframe)
    plt.figure(figsize=(14, 8))
    sns.boxplot(data=exonic_df, x='#CHROM', y='INFO', palette='coolwarm')
    plt.title("Distribution of Exonic Mutation Counts by Chromosome", fontsize=16, fontweight='bold')
    plt.xlabel("Chromosome")
    plt.ylabel("Exonic Mutation Count")
    plt.xticks(rotation=45)


# Bar Plot for Mutation Classification
def plot_mutation_classification_bar(dataframe, threshold=0.37):
    # Filter and classify data
    exonic_info_data = filter_numeric_info(dataframe)
    exonic_info_data_classified = exonic_info_data.apply(
        lambda x: 'Likely Pathogenic/Pathogenic' if x > threshold else 'Likely Benign')
    classification_counts = exonic_info_data_classified.value_counts()

    # Set theme and create plot
    sns.set_theme(style="white")
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x=classification_counts.index, y=classification_counts.values,
                     palette=['#e78ac3', '#8da0cb'])

    # Add count labels inside each bar
    for i, value in enumerate(classification_counts.values):
        ax.text(i, value - 0.1 * value,
                str(value), ha='center', va='top', color='black', fontsize=12, fontweight='bold')

    # Customize plot appearance
    plt.title("Mutation Classifications", fontsize=16, fontweight='bold')
    plt.xlabel("Classification")
    plt.ylabel("Count")
    sns.despine(left=True, bottom=True)


# Generate the PDF
# Generate the PDF with detailed explanations
def generate_pdf_report(dataframe, filename="Mutation_Report.pdf"):
    # Calculate counts
    # Assuming '0' indicates non-exonic
    exonic_info_data = pd.to_numeric(
        dataframe['INFO'].apply(lambda x: x.split(';')[-1].split('=')[-1] if pd.notna(x) else None),
        errors='coerce').dropna()
    non_exonic_count = \
        dataframe['INFO'].apply(lambda x: x.split(';')[-1].split('=')[-1] if pd.notna(x) else None).value_counts()[
            'not_exonic']
    exonic_count = len(exonic_info_data)
    pathogenic_count = (exonic_info_data > 0.37).sum()
    benign_count = exonic_count - pathogenic_count

    # Initialize PDF document
    doc = SimpleDocTemplate(filename, pagesize=letter)
    elements = []
    styles = getSampleStyleSheet()
    elements.append(Paragraph("StrVCTVRE report analysis", styles['Title']))
    elements.append(Spacer(1, 20))

    # Polar Plot Section
    elements.append(Paragraph("Exonic Mutation Distribution Across Chromosomes", styles['Heading2']))
    elements.append(Paragraph(
        f"This polar plot visualizes the distribution of exonic mutations across various chromosomes. "
        f"Out of a total of {exonic_count + non_exonic_count} mutations, {exonic_count} are exonic, spread across chromosomes. "
        f"The largest concentration of mutations is observed in chromosomes with the tallest bars, indicating higher mutation densities.",
        styles['BodyText']))
    buffer = save_plot_to_buffer(lambda: plot_mutation_polar_plot(dataframe), width=8, height=8)
    elements.append(Image(buffer, width=400, height=400))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())
    # Donut Chart Section
    elements.append(Paragraph("Exonic vs Non-Exonic Mutation Distribution", styles['Heading2']))
    elements.append(Paragraph(
        f"This donut chart shows the fraction of exonic versus non-exonic mutations in the dataset. "
        f"Out of {exonic_count + non_exonic_count} total mutations, {exonic_count} are exonic "
        f"while {non_exonic_count} are non-exonic. This chart provides a clear breakdown of mutation types present.",
        styles['BodyText']))
    buffer = save_plot_to_buffer(lambda: plot_mutation_donut_chart(dataframe), width=8, height=8)
    elements.append(Image(buffer, width=400, height=300))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    elements.append(Paragraph("Exonic Mutation Counts by Chromosome", styles['Heading2']))
    elements.append(Paragraph(
        "This box plot illustrates the distribution of exonic mutation counts across various chromosomes, allowing for a comparison of mutation frequencies. "
        "The x-axis represents chromosomes (chr1 to chrY), and the y-axis shows the normalized count of exonic mutations. Each box plot provides insights into "
        "the range, median, quartiles, and outliers of mutation counts per chromosome. Key observations include the following:",
        styles['BodyText']))
    elements.append(Paragraph(box_plot_stats(dataframe), styles['BodyText']))
    buffer = save_plot_to_buffer(lambda: plot_mutation_box_plot(dataframe), width=10, height=6)
    elements.append(Image(buffer, width=450, height=300))
    elements.append(Spacer(1, 20))
    elements.append(PageBreak())

    elements.append(Paragraph("Mutation Classification", styles['Heading2']))
    elements.append(Paragraph(
        f"This bar plot categorizes mutations into 'Likely Pathogenic/Pathogenic' and 'Likely Benign' "
        f"based on a threshold value of 0.37. Out of {exonic_count} exonic mutations, "
        f" <b>{pathogenic_count}</b> are classified as 'Likely Pathogenic/Pathogenic', while <b>{benign_count}</b> "
        f"are classified as 'Likely Benign'.",
        styles['BodyText']))
    buffer = save_plot_to_buffer(lambda: plot_mutation_classification_bar(dataframe), width=8, height=6)
    elements.append(Image(buffer, width=450, height=300))
    elements.append(Spacer(1, 20))

    doc.build(elements)


def start_analysis(patient_data):
    patient_pdf_name = patient_data.split(".")[0] + "_StrVCTVRE" + ".pdf"

    with open(patient_data, 'r') as file:
        for line_num, line in enumerate(file):
            if line.startswith('#CHROM'):
                header_line = line_num
                break

    patient_df = pd.read_csv(patient_data, sep='\t', skiprows=header_line)
    generate_pdf_report(patient_df, patient_pdf_name)
