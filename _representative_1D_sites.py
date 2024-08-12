## Script 4
import pandas as pd
import numpy as np

# takes the frequency differences file and finds the highest frequency change per host for each unique gene
def representative_site_per_host(frequency_diff_B_v_4M, representative_sites_file):
    frequency_diff = pd.read_csv(frequency_diff_B_v_4M, sep="\t")
    
    output_data = []

    # make sure 'site_id' and 'gene_id' columns exist
    required_columns = ['site_id', 'gene_id']
    if not all(col in frequency_diff.columns for col in required_columns):
        raise ValueError(f"Missing required columns: {', '.join(required_columns)}")

    # identify unique gene_ids to process the sites for each gene, not the entire column
    unique_genes = frequency_diff['gene_id'].unique()
    for gene_id in unique_genes:
        # filter df for the current gene_id
        gene_ids_df = frequency_diff[frequency_diff['gene_id'] == gene_id]
        
        # process each column except 'site_id' and 'gene_id'
        for col in gene_ids_df.columns[2:]:
            # replace empty strings with NaN and convert to float
            column_data = gene_ids_df[col].replace('', np.nan).astype(float).abs()
            if not column_data.isnull().all():
                # this finds the index of the maximum value in the current column (the site_id)
                max_index = column_data.idxmax()
                representative_1D = gene_ids_df.at[max_index, 'site_id']
                # append processed data to output_data
                output_data.append({
                    'gene_id': gene_id,
                    'host(dyad)': col,
                    '1D_site_biggest': representative_1D
                })

    # convert to df to then write the data to a file
    representative_sites = pd.DataFrame(output_data)
    representative_sites.to_csv(representative_sites_file, index=False, sep="\t")
    print(f"Processed data written to '{representative_sites_file}'")

# paths
frequency_diff_B_v_4M = 'B_Vulgatus_57955_freq_diff_B_v_4M.txt'
representative_sites_file = 'B_Vulgatus_57955_rep_sites_B_v_4M.txt'

# call the function
representative_site_per_host(frequency_diff_B_v_4M, representative_sites_file)
