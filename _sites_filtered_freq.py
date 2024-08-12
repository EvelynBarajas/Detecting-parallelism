## Script 2
import pandas as pd
import bz2

# get unique samples for B and 4M timepoints
def extract_unique_samples(metadata_pairs_path):
    metadata_pairs_df = pd.read_csv(metadata_pairs_path, sep="\t")
    metadata_pairs_df.columns = metadata_pairs_df.columns.str.strip()
    
    # keep rows for 'B' and '4M' samples
    filtered_df = metadata_pairs_df[(metadata_pairs_df['sample1'] == 'B') | (metadata_pairs_df['sample2'] == '4M')]

    unique_samples = pd.concat([
        filtered_df['accession1'],
        filtered_df['accession2']
    ]).unique()
    
    return set(unique_samples)

# processing by chunks and writing to output file once the chunks are finalized as it goes
def filter_snps_ref_freq(snps_ref_freq_path, good_coverage_sites_path, unique_samples, filtered_sample_freq, chunk_size=10000):
    with open(good_coverage_sites_path, 'r') as file:
        good_coverage_sites_df = pd.read_csv(file, sep="\t")
    good_coverage_sites = set(good_coverage_sites_df['site_id'])
 
    # process by chunks
    with bz2.open(snps_ref_freq_path, 'rt') as file:
        snps_ref_freq_iter = pd.read_csv(file, sep="\t", chunksize=chunk_size)
        
        with open(filtered_sample_freq, 'w') as outfile:
            outfile.write('\t'.join(['site_id', 'gene_id'] + [col for col in snps_ref_freq_iter.get_chunk(0).columns if col in unique_samples]) + '\n')
            
            for snps_ref_freq_df in snps_ref_freq_iter:
                # filter snps_ref_freq based on sites
                filtered_freq_df = snps_ref_freq_df[snps_ref_freq_df['site_id'].isin(good_coverage_sites)]
                
                # merge with good_coverage_sites_df to get gene_id to the final output
                merged_df = pd.merge(filtered_freq_df, good_coverage_sites_df[['site_id', 'gene_id']], on='site_id', how='left')
                
                # filter columns based on unique samples
                filtered_columns = ['site_id', 'gene_id'] + [col for col in snps_ref_freq_df.columns if col in unique_samples]
                filtered_freq_df = merged_df[filtered_columns]
                filtered_freq_df.to_csv(outfile, index=False, header=False, sep="\t")

# paths
snps_ref_freq_path = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/MegaMI/data/snps/Bacteroides_vulgatus_57955/snps_ref_freq.txt.bz2'
good_coverage_sites_path = 'good_coverage_1Dsites_B_vulgatus_57955.txt'
metadata_pairs_path = 'metadata_pairs.txt'
filtered_sample_freq = 'B_4M_samples_filtered_freqs.txt'

# get unique samples from metadata pairs
unique_samples = extract_unique_samples(metadata_pairs_path)
filter_snps_ref_freq(snps_ref_freq_path, good_coverage_sites_path, unique_samples, filtered_sample_freq)
