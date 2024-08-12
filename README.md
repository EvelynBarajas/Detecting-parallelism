# Detecting-parallelism

### Script 1: _get_gene_sites.py
- Identifies unique gene IDs from snps_info.txt.bz2
- Filters SNPs info by gene_id and site_type (1D only)
- Filters sites by site_ids in site_datapoints_Bacteroides_vulgatus_57955_B_v_4M.txt
- Outputs file good_coverage_1Dsites_B_vulgatus_57955.txt in the format: site_id		gene_id

### Script 2: _sites_filtered_freq.py
- Identifies unique samples from the metada_pairs.txt file. Since the previous script filters sites depending on a single comparison (B_v_4M), this function uses samples that are B and 4M.
- Filters snps_ref_fre.txt.bz2 file based on the good coverage sites identified in the previous script (good_coverage_1Dsites_B_vulgatus_57955.txt) and B and 4M samples
- Outputs B_v_4M_samples_filtered_freqs.txt file in the format: site_id	gene_id	sample1	sample2

### Script 3: _match_freq_diff.py
- Again using metadata_pairs.txt file, it matches the corresponding B and 4M samples when the ‘map’ column has rows with ‘B_v_4M’. 
- Using the output file (B_4M_samples_filtered_freqs.txt) and the matched sample pairs, it calculates the frequency difference for each pair at each site_id and formats the column header as '{sample1}, {sample2}'
- Outputs file B_Vulgatus_57955_freq_diff_B_v_4M.txt in the format: site_id	gene_id	ERR525910, ERR525909

### Script 4: _representative_1D_sites.py
- Takes file B_Vulgatus_57955_freq_diff_B_v_4M.txt to identify the site_id with the highest frequency change in each host on all genes.
- It first identifies unique genes to then identify the site_ids in the current gene it is working with
- Outputs file B_Vulgatus_57955_rep_sites_B_v_4M.txt with this format: gene_id	host(dyad)	1D_site_biggest
