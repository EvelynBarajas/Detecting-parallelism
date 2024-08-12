## Script 1
import pandas as pd
import bz2

# function to get unique gene IDs from snps_info.txt.bz2
def get_unique_gene_ids(file_path):
    with bz2.open(file_path, 'rt') as file:
        df = pd.read_csv(file, sep="\t")
    unique_gene_ids = df['gene_id'].unique()
    return unique_gene_ids, df

# function to filter SNPs info by gene_id and site_type
def filter_snps_info(df, gene_id, site_type):
    filtered_sites = df[(df['gene_id'] == gene_id) & (df['site_type'] == site_type)]
    return filtered_sites['site_id']

# function to filter site_datapoints by site_ids
def filter_site_datapoints(site_datapoints_path, site_ids):
    site_datapoints_df = pd.read_csv(site_datapoints_path, sep="\t")
    filtered_site_datapoints_df = site_datapoints_df[site_datapoints_df['SITE_ID'].isin(site_ids)]
    return filtered_site_datapoints_df['SITE_ID']

# paths to files
snps_info_path = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/MegaMI/data/snps/Bacteroides_vulgatus_57955/snps_info.txt.bz2'
site_datapoints_path = '/u/project/ngarud/amzurita/ParallelEvolution_GLMM/Analysis/SpeciesData/MostPrevalent_MinReads50/Bacteroides_vulgatus_57955/site_datapoints_Bacteroides_vulgatus_57955_B_v_4M.txt'

# get all unique gene IDs and the DataFrame from snps_info.txt.bz2
unique_gene_ids, snps_info_df = get_unique_gene_ids(snps_info_path)
print(f"Number of unique genes: {len(unique_gene_ids)}")

# output file and write the header, write genes as the process is running
output_file = 'good_coverage_1Dsites_B_vulgatus_57955.txt'
with open(output_file, 'w') as file:
    file.write('site_id\tgene_id\n')

# apply filters to each unique genes
for gene_id in unique_gene_ids:
    df_site_ids = filter_snps_info(snps_info_df, gene_id, site_type='1D')
    good_coverage_sites = filter_site_datapoints(site_datapoints_path, df_site_ids)
    
    # append results to the output file
    with open(output_file, 'a') as file:
        for site_id in good_coverage_sites:
            file.write(f"{site_id}\t{gene_id}\n")
print(f"Results saved to {output_file}")
