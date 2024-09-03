###  Script 4 Update  ###
## Use site_datapoints_file to obtain list of good coverage pairs for that site 
## and get the representative 1D site based on the listed pairs
import pandas as pd
import numpy as np
import csv
import ast

# paths
site_datapoints_file = '/u/project/ngarud/amzurita/ParallelEvolution_GLMM/Analysis/SpeciesData/MostPrevalent_MinReads50/Bacteroides_vulgatus_57955/site_datapoints_Bacteroides_vulgatus_57955_B_v_4M.txt'
freq_diff_file = 'B_Vulgatus_57955_freq_diff_B_v_4M.txt'
output_file = 'NEW_B_Vulgatus_57955_rep_sites_B_v_4M.txt'

# use site_datapoints_file, create a dictionary for SITE_ID and LIST_PAIRS_PASS
def load_site_pairs(site_datapoints_file):
    site_pairs_dict = {}
    with open(site_datapoints_file, 'r') as site_file:
        site_reader = csv.reader(site_file, delimiter='\t')
        next(site_reader)  # skip header
        for line in site_reader:
            site_id = line[0]
            list_pairs_pass = ast.literal_eval(line[3])
            site_pairs_dict[site_id] = list_pairs_pass
    return site_pairs_dict

# process the frequency difference file in chunks, filter based on site_pairs_dict
def process_frequency_diff(freq_diff_file, site_pairs_dict, output_file):
    with open(freq_diff_file, 'r') as freq_file, open(output_file, 'w', newline='') as out_file:
        freq_reader = csv.reader(freq_file, delimiter='\t')
        writer = csv.writer(out_file, delimiter='\t')
        
        # write the header over to the output file
        freq_header = next(freq_reader)
        writer.writerow(freq_header)
        
        output_data = []

        for freq_line in freq_reader:
            site_id = freq_line[0]
            filtered_line = freq_line[:2]  # keep site_id and gene_id
            
            # check if site_id is in site_pairs_dict
            if site_id in site_pairs_dict:
                valid_pairs = site_pairs_dict[site_id]
                
                # process each pair and check if they're valid
                for col_index, col in enumerate(freq_line[2:], start=2):
                    pair = freq_header[col_index].split(', ')
                    if list(pair) in valid_pairs:
                        filtered_line.append(col)
                    else:
                        filtered_line.append(np.nan)
            else:
                # if site_id not found, set all pairs to nan
                filtered_line.extend([np.nan] * (len(freq_line) - 2))
        
            output_data.append(filtered_line)

        return output_data, freq_header

# identify the representative site for each host
def representative_site_per_host(output_data, freq_header, output_file):
    frequency_diff = pd.DataFrame(output_data, columns=freq_header)
    output_data_final = []

    # find unique gene_ids to process the sites for each gene, not the entire column
    unique_genes = frequency_diff['gene_id'].unique()
    for gene_id in unique_genes:
        # filter df for the current gene_id
        gene_ids_df = frequency_diff[frequency_diff['gene_id'] == gene_id]

        # process each column except 'site_id' and 'gene_id'
        for col in gene_ids_df.columns[2:]:
            # replace empty strings with NaN and convert to float
            column_data = gene_ids_df[col].replace('', np.nan).astype(float).abs()
            if not column_data.isnull().all():
                # find the site_id of the maximum value in the current column
                max_index = column_data.idxmax()
                representative_1D = gene_ids_df.at[max_index, 'site_id']
                # append data to output_data as it gets processed
                output_data_final.append({
                    'gene_id': gene_id,
                    'host(dyad)': col,
                    '1D_site_biggest': representative_1D
                })

    # write df to a file
    representative_sites = pd.DataFrame(output_data_final)
    representative_sites.to_csv(output_file, index=False, sep="\t")
    print(f"Processed data written to '{output_file}'")

# main, run process
def main():
    # step 1
    site_pairs_dict = load_site_pairs(site_datapoints_file)
    # step 2
    output_data, freq_header = process_frequency_diff(freq_diff_file, site_pairs_dict, output_file)
    # step 3
    representative_site_per_host(output_data, freq_header, output_file)

# call the main function
if __name__ == "__main__":
    main()
