import json
import seaborn as sns
import pandas as pd

strvctr_data = {}


def filter_numeric_info(dataframe):
    filtered_infoData = dataframe['INFO'].apply(lambda x: x.split(';')[-1] if pd.notna(x) else None).apply(
        lambda x: x.split('=')[-1] if pd.notna(x) else None)
    return pd.to_numeric(filtered_infoData, errors='coerce').dropna()


# Polar Plot Function
def get_polar_plot_data(dataframe):
    allowed_chromosomes = ["chr" + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df = dataframe[dataframe['#CHROM'].isin(allowed_chromosomes)]
    filtered_infoData = filter_numeric_info(filtered_df)
    chromosome_names = sorted(filtered_df['#CHROM'].unique(),
                              key=lambda x: (int(x[3:]) if x[3:].isdigit() else float('inf'), x))
    radii = [int(filtered_infoData[filtered_df['#CHROM'] == chrom].count()) for chrom in chromosome_names]

    data = {'label': chromosome_names,
            'data': radii}
    strvctr_data.update(PolarPlot=data)
    return


# Donut Plot for Exonic vs Non-Exonic Mutation Distribution
def get_donut_plot_data(dataframe):
    filtered_infoData = dataframe['INFO'].apply(lambda x: x.split(';')[-1] if pd.notna(x) else None).apply(
        lambda x: x.split('=')[-1] if pd.notna(x) else None)
    non_exonic_count = filtered_infoData.value_counts()['not_exonic']
    exonic_info_data = pd.to_numeric(filtered_infoData, errors='coerce').dropna()
    exonic_count = len(exonic_info_data)
    counts = pd.DataFrame({'Mutation Type': ['Non-exonic', 'Exonic'], 'Count': [non_exonic_count, exonic_count]})
    data = counts['Count'].values
    labels = counts['Mutation Type'].values

    data = {'label': labels.tolist(),
            'data': data.tolist()}
    strvctr_data.update(Donut=data)
    return data


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


def get_box_plot_data(dataframe):
    sorted_chromosomes, exonic_df = filter_info(dataframe)

    grouped = exonic_df.groupby('#CHROM')

    chromosome_names = list(grouped.groups.keys())

    info_data = [group['INFO'].tolist() for _, group in grouped]
    data = {'label': chromosome_names,
            'data': info_data}
    strvctr_data.update(BoxPlot=data)
    return data


def get_bar_plot_data(dataframe, threshold=0.37):
    exonic_info_data = filter_numeric_info(dataframe)
    exonic_info_data_classified = exonic_info_data.apply(
        lambda x: 'Likely Pathogenic' if x > threshold else 'Likely Benign')
    classification_counts = exonic_info_data_classified.value_counts()

    data = {'label': classification_counts.index.tolist(),
            'data': classification_counts.values.tolist()}
    dataBP = {'label': ['Data classification'],
            'data': [exonic_info_data.tolist()]}

    strvctr_data.update(BarPlot=data)
    strvctr_data.update(BarPlotBox=dataBP)
    return [data, dataBP]

patient_df = None
def start_analysis(patient_data,threshold=0.37):
    global patient_df
    with open(patient_data, 'r') as file:
        for line_num, line in enumerate(file):
            if line.startswith('#CHROM'):
                header_line = line_num
                break

    patient_df = pd.read_csv(patient_data, sep='\t', skiprows=header_line)
    # get_polar_plot_data(patient_df)
    get_donut_plot_data(patient_df)
    get_box_plot_data(patient_df)
    get_bar_plot_data(patient_df, threshold=threshold)

    out_file = open("../strvctrvre.json", "w")

    json.dump(strvctr_data, out_file, indent=4)
    out_file.close()
    return strvctr_data



start_analysis("/Users/martindruzbacky/PycharmProjects/DP_helper/Data/patient1_annotated.vcf")