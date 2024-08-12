## Script 3
import pandas as pd

# match samples to obtain corresponding pais and get the frequency differences
def calculate_frequency_differences(metadata_pairs_path, filtered_freq_path, frequency_diff_B_v_4M):
    metadata_pairs_df = pd.read_csv(metadata_pairs_path, sep="\t")
    filtered_metadata_pairs_df = metadata_pairs_df[metadata_pairs_df['map'] == 'B_v_4M'] # change depending on desired comparison
    
    # dictionary to hold frequency differences
    frequency_diff_dict = {
        'site_id': [],
        'gene_id': []
    }
    
    # read the filtered frequency data
    with open(filtered_freq_path, 'r') as file:
        header = file.readline().strip().split('\t')
        sample_columns = header[2:]  # sample start from the third column

        # process line by line
        for line in file:
            line_data = line.strip().split('\t')
            site_id = line_data[0]
            gene_id = line_data[1]
            frequencies = list(map(float, line_data[2:]))
            
            frequency_diff_dict['site_id'].append(site_id)
            frequency_diff_dict['gene_id'].append(gene_id)
            
            # iterate through the filtered metadata pairs
            for _, row in filtered_metadata_pairs_df.iterrows():
                sample1 = row['accession1']
                sample2 = row['accession2']
                sample_pair = f'{sample1}, {sample2}'
    
                if sample1 in sample_columns and sample2 in sample_columns:
                    index1 = sample_columns.index(sample1)
                    index2 = sample_columns.index(sample2)
                    freq_diff = abs(frequencies[index2] - frequencies[index1])
                    if sample_pair not in frequency_diff_dict:
                        frequency_diff_dict[sample_pair] = []
                    frequency_diff_dict[sample_pair].append(freq_diff)
    
    # make dictionary into a DataFrame to be able to write output file
    frequency_diff_df = pd.DataFrame(frequency_diff_dict)
    frequency_diff_df.to_csv(frequency_diff_B_v_4M, index=False, sep="\t")

# paths
metadata_pairs_path = 'metadata_pairs.txt'
filtered_freq_path = 'B_4M_samples_filtered_freqs.txt'
frequency_diff_B_v_4M = 'B_Vulgatus_57955_freq_diff_B_v_4M.txt'

calculate_frequency_differences(metadata_pairs_path, filtered_freq_path, frequency_diff_B_v_4M)
