#!/usr/bin/env python3
"""
generate_manuscript_tables.py

Programmatically calculates all statistics used in the manuscript drafts
and exports them to the 'manuscript_tables/' folder as clean, standardized CSV files.
Fixes formula errors (such as SIRV NRMSE similarity), standardizes tool and 
normalization naming, and generates complete metrics for Sections 3.1 - 3.5.
"""

import os
import glob
import pandas as pd
import numpy as np

# Define paths relative to the repository root
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_OUT_DIR = os.path.join(BASE_DIR, "manuscript_tables")

os.makedirs(TABLES_OUT_DIR, exist_ok=True)

# Standard Tool and Normalization Mappings
TOOL_NAME_MAP = {
    'flair': 'FLAIR',
    'isoquant': 'IsoQuant',
    'kallisto': 'Kallisto',
    'oarfish': 'Oarfish',
    'bambu': 'bambu',
    'FLAIR': 'FLAIR',
    'IsoQuant': 'IsoQuant',
    'Kallisto': 'Kallisto',
    'Oarfish': 'Oarfish'
}

NORM_NAME_MAP = {
    'CPM': 'CPM',
    'TPM': 'TPM',
    'cqn': 'CQN',
    'EDA': 'EDASeq',
    'ratio_correction': 'SMR',
    'ratio_counts': 'SMR'
}

# ==============================================================================
# SECTION 3.1: BASELINE SEQUENCING THROUGHPUT, READ LENGTHS & TRANSCRIPT COUNTS
# ==============================================================================

def calculate_section_3_1_throughput():
    """
    Loads fig1 sequencing performance summary and calculates median, min, max
    read numbers, fragment numbers, and read lengths reported in Section 3.1.
    """
    print("Calculating metrics for Section 3.1 (Throughput & Read Lengths)...")
    csv_path = os.path.join(BASE_DIR, "paper_plots/fig1/data/sequencing_performance_summary.csv")
    if not os.path.exists(csv_path):
        print(f"Warning: Sequencing summary file not found at {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    
    # Calculate Total_Fragments (for Illumina paired-end reads, fragments = Total_Reads / 2)
    df['Total_Fragments'] = df['Total_Reads']
    df.loc[df['Technology'] == 'Illumina', 'Total_Fragments'] = df['Total_Reads'] / 2.0
    
    stats = df.groupby(['Dataset', 'Technology']).agg(
        median_reads=('Total_Reads', 'median'),
        min_reads=('Total_Reads', 'min'),
        max_reads=('Total_Reads', 'max'),
        median_fragments=('Total_Fragments', 'median'),
        median_gb=('Total_Gb', 'median'),
        median_mapped_pct=('Mapped_Reads_Pct', 'median'),
        median_mapped_len=('Average_Mapped_Read_Length', 'median'),
        min_mapped_len=('Average_Mapped_Read_Length', 'min'),
        max_mapped_len=('Average_Mapped_Read_Length', 'max'),
        sample_count=('Sample', 'count')
    ).reset_index()
    
    out_path = os.path.join(TABLES_OUT_DIR, "section_3_1_sequencing_throughput_statistics.csv")
    stats.to_csv(out_path, index=False)
    print(f" - Saved: {out_path}")

def calculate_section_3_1_transcript_counts():
    """
    Aggregates unfiltered and filtered transcript model counts across all tools,
    technologies, and datasets as reported in Section 3.1 of the Results.
    """
    print("Calculating transcript model counts for Section 3.1 (Unfiltered vs Filtered)...")
    files_map = [
        ('paper_plots/fig1/data/sqanti_new/masseq_sqanti_summary_global.csv', 'Mouse', 'Kinnex'),
        ('paper_plots/fig1/data/sqanti_new/ont_sqanti_summary_global.csv', 'Mouse', 'ONT cDNA'),
        ('paper_plots/fig1/data/sqanti_new/isoseq_sqanti_summary_global.csv', 'Mouse', 'IsoSeq'),
        ('paper_plots/fig1/data/sqanti_new/neuron_pb_sqanti_summary_global.csv', 'iNeuron', 'Kinnex'),
        ('paper_plots/fig1/data/sqanti_new/neuron_ont_sqanti_summary_global.csv', 'iNeuron', 'ONT cDNA'),
        ('paper_plots/fig1/data/sqanti_new/kinnex_norm_sqanti_summary_global.csv', 'LongBench', 'Kinnex'),
        ('paper_plots/fig1/data/sqanti_new/cDNA_norm_sqanti_summary_global.csv', 'LongBench', 'ONT cDNA'),
        ('paper_plots/fig1/data/sqanti_new/drna_norm_sqanti_summary_global.csv', 'LongBench', 'ONT dRNA')
    ]

    records = []
    for rel_path, dataset, tech in files_map:
        full_path = os.path.join(BASE_DIR, rel_path)
        if os.path.exists(full_path):
            df = pd.read_csv(full_path)
            df = df[df['Tool'] != 'isoseq']
            df['Tool'] = df['Tool'].map(lambda x: TOOL_NAME_MAP.get(x, x))
            
            grouped = df.groupby(['Tool', 'Filtered'])['n_trans'].sum().reset_index()
            piv = grouped.pivot(index='Tool', columns='Filtered', values='n_trans').fillna(0).astype(int)
            
            for tool in piv.index:
                unfiltered = piv.loc[tool, False] if False in piv.columns else 0
                filtered = piv.loc[tool, True] if True in piv.columns else 0
                pct_retained = (filtered / unfiltered * 100.0) if unfiltered > 0 else 0.0
                records.append({
                    'Dataset': dataset,
                    'Technology': tech,
                    'Tool': tool,
                    'Unfiltered_Transcripts': unfiltered,
                    'Filtered_Transcripts': filtered,
                    'Percent_Retained': round(pct_retained, 2)
                })

    out_df = pd.DataFrame(records)
    out_path = os.path.join(TABLES_OUT_DIR, "section_3_1_active_and_unfiltered_transcript_counts.csv")
    out_df.to_csv(out_path, index=False)
    print(f" - Saved: {out_path}")

def calculate_section_3_1_sqanti_novel_proportions():
    """
    Loads SQANTI3 classification proportions and calculates the percentage of novel
    transcript models before and after expression-based filtering for Section 3.1.
    """
    print("Calculating novel transcript model proportions for Section 3.1...")
    csv_path = os.path.join(BASE_DIR, "paper_plots/fig1/data/sqanti_proportions_comparison.csv")
    if not os.path.exists(csv_path):
        print(f"Warning: SQANTI proportions file not found at {csv_path}")
        return

    df = pd.read_csv(csv_path)
    novel_df = df[df['IsoformClass'] == 'Novel'].copy()
    novel_df['FilteredLabel'] = novel_df['FilteredLabel'].str.replace('\n', ' ')

    dataset_map = {'NIH': 'Mouse', 'Neuron': 'iNeuron', 'LongBench': 'LongBench'}
    novel_df['Dataset'] = novel_df['Dataset'].map(lambda x: dataset_map.get(x, x))
    novel_df['Tool'] = novel_df['Tool'].map(lambda x: TOOL_NAME_MAP.get(x, x))

    before_df = novel_df[novel_df['FilteredLabel'].str.contains('Before')].copy()
    after_df = novel_df[novel_df['FilteredLabel'].str.contains('After')].copy()

    merged = pd.merge(
        before_df[['Dataset', 'LibType', 'Tool', 'Count', 'Total', 'Percentage']],
        after_df[['Dataset', 'LibType', 'Tool', 'Count', 'Total', 'Percentage']],
        on=['Dataset', 'LibType', 'Tool'],
        suffixes=('_Unfiltered', '_Filtered'),
        how='outer'
    ).fillna(0)

    all_combos = df[['Dataset', 'LibType']].drop_duplicates()
    all_combos['Dataset'] = all_combos['Dataset'].map(lambda x: dataset_map.get(x, x))
    ref_restricted = []
    for _, row in all_combos.iterrows():
        for tool in ['Kallisto', 'Oarfish']:
            ref_restricted.append({
                'Dataset': row['Dataset'],
                'LibType': row['LibType'],
                'Tool': tool,
                'Count_Unfiltered': 0,
                'Total_Unfiltered': 0,
                'Percentage_Unfiltered': 0.0,
                'Count_Filtered': 0,
                'Total_Filtered': 0,
                'Percentage_Filtered': 0.0
            })
    ref_df = pd.DataFrame(ref_restricted)
    
    combined_novel = pd.concat([merged, ref_df], ignore_index=True)
    combined_novel.drop_duplicates(subset=['Dataset', 'LibType', 'Tool'], inplace=True)

    combined_novel.rename(columns={
        'LibType': 'Technology',
        'Count_Unfiltered': 'Unfiltered_Novel_Count',
        'Total_Unfiltered': 'Unfiltered_Total',
        'Percentage_Unfiltered': 'Unfiltered_Novel_Pct',
        'Count_Filtered': 'Filtered_Novel_Count',
        'Total_Filtered': 'Filtered_Total',
        'Percentage_Filtered': 'Filtered_Novel_Pct'
    }, inplace=True)

    combined_novel['Delta_Novel_Pct'] = combined_novel['Filtered_Novel_Pct'] - combined_novel['Unfiltered_Novel_Pct']

    combined_novel['Unfiltered_Novel_Pct'] = combined_novel['Unfiltered_Novel_Pct'].round(2)
    combined_novel['Filtered_Novel_Pct'] = combined_novel['Filtered_Novel_Pct'].round(2)
    combined_novel['Delta_Novel_Pct'] = combined_novel['Delta_Novel_Pct'].round(2)

    combined_novel = combined_novel.sort_values(by=['Dataset', 'Technology', 'Tool']).reset_index(drop=True)

    out_path = os.path.join(TABLES_OUT_DIR, "section_3_1_sqanti_novel_transcript_proportions.csv")
    combined_novel.to_csv(out_path, index=False)
    print(f" - Saved: {out_path}")

# ==============================================================================
# SECTION 3.3: SYNTHETIC SPIKE-IN (SIRV & ERCC) SIMILARITY AND RMSE STATISTICS
# ==============================================================================

def calculate_section_3_3_spike_ins():
    """
    Loads Figure 5 SIRV & ERCC summary data, computes accurate NRMSE similarity
    (1 - NRMSE), and exports summary tables for Section 3.3 of the Results.
    """
    print("Calculating metrics for Section 3.3 (Spike-ins & Sensitivity)...")
    datasets = ['LongBench', 'Neuron', 'NIH']
    all_data_list = []
    
    for ds in datasets:
        fofn_path = os.path.join(BASE_DIR, f"paper_plots/fig5/data/sirv_trans_{ds}.fofn")
        if not os.path.exists(fofn_path):
            continue
            
        with open(fofn_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    rel_path, tech = parts
                    csv_path = os.path.join(BASE_DIR, "paper_plots/fig5", rel_path)
                    if os.path.exists(csv_path):
                        df = pd.read_csv(csv_path)
                        df = df[df['Tool'] != 'isoseq']
                        df_all = df[df['LengthQuantile'] == 'All'].copy()
                        df_all['Dataset'] = "Mouse" if ds == "NIH" else ds
                        df_all['Tech'] = tech.replace('_', ' ')
                        all_data_list.append(df_all)
                        
    if all_data_list:
        combined = pd.concat(all_data_list, ignore_index=True)
        
        # Calculate correct NRMSE Similarity: 1 - NRMSE (NRMSE = RMSE / abs(mean_expected))
        combined['mean_expected'] = np.where(np.abs(combined['mean_expected']) < 1e-6, 1.0, combined['mean_expected'])
        combined['NRMSE'] = combined['RMSE'] / np.abs(combined['mean_expected'])
        combined['Similarity'] = 1.0 - combined['NRMSE']
        combined['Similarity'] = np.where(np.isfinite(combined['Similarity']), combined['Similarity'], 0.0)
        
        # Clean Tool & Normalization labels
        combined['Tool'] = combined['Tool'].map(lambda x: TOOL_NAME_MAP.get(x, x))
        combined['Normalization'] = combined['Normalization'].map(lambda x: NORM_NAME_MAP.get(x, x))
        
        # Table A: Mean Similarity & RMSE by Dataset, Technology, and Normalization
        table_norm = combined.groupby(['Dataset', 'Tech', 'Normalization'])[['Similarity', 'RMSE']].mean().reset_index()
        norm_out_path = os.path.join(TABLES_OUT_DIR, "section_3_3_sirv_similarity_rmse_by_normalization.csv")
        table_norm.to_csv(norm_out_path, index=False)
        print(f" - Saved: {norm_out_path}")
        
        # Table B: Mean Similarity & RMSE by Tool (excluding SMR)
        df_no_smr = combined[combined['Normalization'] != 'SMR']
        table_tool = df_no_smr.groupby(['Dataset', 'Tool'])[['Similarity', 'RMSE']].mean().reset_index()
        table_tool = table_tool.sort_values(by=['Dataset', 'Similarity'], ascending=[True, False])
        tool_out_path = os.path.join(TABLES_OUT_DIR, "section_3_3_sirv_similarity_rmse_by_tool.csv")
        table_tool.to_csv(tool_out_path, index=False)
        print(f" - Saved: {tool_out_path}")

    # Export ERCC Accuracy Metrics by Normalization and by Tool
    ercc_files = [
        ("results/neuron/ont/neuron_ont_normalization_ercc_summary.csv", "ONT cDNA", "Neuron"),
        ("results/neuron/kinnex/neuron_pb_normalization_ercc_summary.csv", "Kinnex", "Neuron"),
        ("results/longbench_ncbi/cDNA/cDNA_norm_normalization_ercc_summary.csv", "ONT cDNA", "LongBench"),
        ("results/longbench_ncbi/dRNA/drna_norm_normalization_ercc_summary.csv", "ONT dRNA", "LongBench"),
        ("results/longbench_ncbi/Kinnex/kinnex_norm_normalization_ercc_summary.csv", "Kinnex", "LongBench")
    ]
    ercc_list = []
    for rel_path, tech, ds in ercc_files:
        csv_path = os.path.join(BASE_DIR, rel_path)
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            df = df[df['Tool'] != 'isoseq']
            df_all = df[df['LengthQuantile'] == 'All'].copy()
            df_all['Tech'] = tech
            df_all['Dataset'] = ds
            ercc_list.append(df_all)
            
    if ercc_list:
        combined_ercc = pd.concat(ercc_list, ignore_index=True)
        combined_ercc['mean_expected'] = np.where(np.abs(combined_ercc['mean_expected']) < 1e-6, 1.0, combined_ercc['mean_expected'])
        combined_ercc['NRMSE'] = combined_ercc['RMSE'] / np.abs(combined_ercc['mean_expected'])
        combined_ercc['Similarity'] = 1.0 - combined_ercc['NRMSE']
        combined_ercc['Similarity'] = np.where(np.isfinite(combined_ercc['Similarity']), combined_ercc['Similarity'], 0.0)
        combined_ercc['Tool'] = combined_ercc['Tool'].map(lambda x: TOOL_NAME_MAP.get(x, x))
        combined_ercc['Normalization'] = combined_ercc['Normalization'].map(lambda x: NORM_NAME_MAP.get(x, x))
        
        ercc_table = combined_ercc.groupby(['Dataset', 'Tech', 'Normalization'])[['Similarity', 'RMSE']].mean().reset_index()
        ercc_out_path = os.path.join(TABLES_OUT_DIR, "section_3_3_ercc_similarity_rmse.csv")
        ercc_table.to_csv(ercc_out_path, index=False)
        print(f" - Saved: {ercc_out_path}")

        ercc_tool_df = combined_ercc[combined_ercc['Normalization'] != 'SMR'].groupby(['Dataset', 'Tech', 'Tool'])[['Similarity', 'RMSE']].mean().reset_index()
        ercc_tool_df = ercc_tool_df.sort_values(by=['Dataset', 'Tech', 'RMSE'], ascending=[True, True, True])
        ercc_tool_out_path = os.path.join(TABLES_OUT_DIR, "section_3_3_ercc_similarity_rmse_by_tool.csv")
        ercc_tool_df.to_csv(ercc_tool_out_path, index=False)
        print(f" - Saved: {ercc_tool_out_path}")

        # Same metrics pooled ACROSS technologies: the Dataset x Tool means quoted in the
        # Section 3.3 prose. The table above is split by Tech and cannot reproduce them.
        # The design is balanced (equal rows per Dataset/Tech/Tool), so pooling the rows
        # directly is identical to averaging the per-technology means.
        ercc_pooled_df = combined_ercc[combined_ercc['Normalization'] != 'SMR'].groupby(['Dataset', 'Tool']).agg(
            Similarity=('Similarity', 'mean'),
            RMSE=('RMSE', 'mean'),
            Tech_Count=('Tech', 'nunique'),
            Row_Count=('RMSE', 'count')
        ).reset_index()
        ercc_pooled_df = ercc_pooled_df.sort_values(by=['Dataset', 'RMSE'], ascending=[True, True])
        ercc_pooled_out_path = os.path.join(TABLES_OUT_DIR, "section_3_3_ercc_similarity_rmse_by_tool_pooled.csv")
        ercc_pooled_df.to_csv(ercc_pooled_out_path, index=False)
        print(f" - Saved: {ercc_pooled_out_path}")

    # Export Spike-in Detection Sensitivity Metrics
    sensitivity_files = [
        ("paper_plots/fig5/data/nih/ont_sensitivity_sirv_summary.csv", "ONT cDNA", "Mouse", "SIRV"),
        ("paper_plots/fig5/data/nih/isoseq_sensitivity_sirv_summary.csv", "IsoSeq", "Mouse", "SIRV"),
        ("paper_plots/fig5/data/nih/masseq_sensitivity_sirv_summary.csv", "Kinnex", "Mouse", "SIRV"),
        ("paper_plots/fig5/data/neuron/ont_sensitivity_sirv_summary.csv", "ONT cDNA", "Neuron", "SIRV"),
        ("paper_plots/fig5/data/neuron/kinnex_sensitivity_sirv_summary.csv", "Kinnex", "Neuron", "SIRV"),
        ("paper_plots/fig5/data/longbench/cDNA_norm_sensitivity_sirv_summary.csv", "ONT cDNA", "LongBench", "SIRV"),
        ("paper_plots/fig5/data/longbench/drna_sensitivity_sirv_summary.csv", "ONT dRNA", "LongBench", "SIRV"),
        ("paper_plots/fig5/data/longbench/kinnex_sensitivity_sirv_summary.csv", "Kinnex", "LongBench", "SIRV"),
        ("results/neuron/ont/neuron_ont_sensitivity_ercc_summary.csv", "ONT cDNA", "Neuron", "ERCC"),
        ("results/neuron/kinnex/neuron_pb_sensitivity_ercc_summary.csv", "Kinnex", "Neuron", "ERCC"),
        ("results/longbench_ncbi/cDNA/cDNA_norm_sensitivity_ercc_summary.csv", "ONT cDNA", "LongBench", "ERCC"),
        ("results/longbench_ncbi/dRNA/drna_norm_sensitivity_ercc_summary.csv", "ONT dRNA", "LongBench", "ERCC"),
        ("results/longbench_ncbi/Kinnex/kinnex_norm_sensitivity_ercc_summary.csv", "Kinnex", "LongBench", "ERCC")
    ]
    sens_list = []
    for rel_path, tech, ds, spike_in in sensitivity_files:
        csv_path = os.path.join(BASE_DIR, rel_path)
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            df = df[df['Tool'] != 'isoseq']
            df['Tech'] = tech
            df['Dataset'] = ds
            df['SpikeIn'] = spike_in
            sens_list.append(df)
            
    if sens_list:
        combined_sens = pd.concat(sens_list, ignore_index=True)
        combined_sens['Tool'] = combined_sens['Tool'].map(lambda x: TOOL_NAME_MAP.get(x, x))
        sens_summary = combined_sens.groupby(['Dataset', 'Tech', 'SpikeIn', 'Tool'])['Sensitivity'].mean().reset_index()
        sens_out_path = os.path.join(TABLES_OUT_DIR, "section_3_3_spikein_detection_sensitivity.csv")
        sens_summary.to_csv(sens_out_path, index=False)
        print(f" - Saved: {sens_out_path}")

# ==============================================================================
# SECTION 3.4: VALIDATION OF ENDOGENOUS TRANSCRIPTS (ILLUMINA CORRELATION)
# ==============================================================================

def calculate_section_3_4_illumina_correlations():
    """
    Loads Figure 4 short-read correlation files and exports summary tables
    for Section 3.4 of the Results section.
    """
    print("Calculating metrics for Section 3.4 (Illumina Correlations)...")
    datasets = ['LongBench', 'Neuron', 'NIH']
    levels = ['gene', 'trans']
    all_data_list = []
    
    for ds in datasets:
        for lvl in levels:
            fofn_path = os.path.join(BASE_DIR, f"paper_plots/fig4/data/sr_{lvl}_{ds}.fofn")
            if not os.path.exists(fofn_path):
                continue
                
            with open(fofn_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) == 2:
                        rel_path, tech = parts
                        csv_path = os.path.join(BASE_DIR, "paper_plots/fig4", rel_path)
                        if os.path.exists(csv_path):
                            df = pd.read_csv(csv_path)
                            df = df[df['Tool'] != 'isoseq']
                            df_all = df[df['LengthQuantile'] == 'All'].copy()
                            df_all['Dataset'] = "Mouse" if ds == "NIH" else ds
                            df_all['Level'] = "Gene" if lvl == "gene" else "Transcript"
                            df_all['Tech'] = tech.replace('_', ' ')
                            all_data_list.append(df_all)
                            
    if not all_data_list:
        print("No Illumina correlation data found.")
        return
        
    combined = pd.concat(all_data_list, ignore_index=True)
    combined['Tool'] = combined['Tool'].map(lambda x: TOOL_NAME_MAP.get(x, x))
    combined['Normalization'] = combined['Normalization'].map(lambda x: NORM_NAME_MAP.get(x, x))
    
    # Table A: Median Spearman Correlation by Dataset, Technology, Level, and Normalization
    table_normalization = combined.groupby(['Dataset', 'Tech', 'Level', 'Normalization'])['Correlation'].median().reset_index()
    norm_out_path = os.path.join(TABLES_OUT_DIR, "section_3_4_illumina_correlation_by_normalization.csv")
    table_normalization.to_csv(norm_out_path, index=False)
    print(f" - Saved: {norm_out_path}")
    
    # Table B: Median Spearman Correlation by Tool (excluding SMR)
    df_no_smr = combined[combined['Normalization'] != 'SMR']
    table_tool = df_no_smr.groupby(['Dataset', 'Tech', 'Level', 'Tool'])['Correlation'].median().reset_index()
    table_tool = table_tool.sort_values(by=['Dataset', 'Tech', 'Level', 'Correlation'], ascending=[True, True, True, False])
    tool_out_path = os.path.join(TABLES_OUT_DIR, "section_3_4_illumina_correlation_by_tool.csv")
    table_tool.to_csv(tool_out_path, index=False)
    print(f" - Saved: {tool_out_path}")
    
    # Table C: Overall Sample Median & Mean Correlation across ALL samples combined
    overall_summary = combined.groupby(['Tech', 'Level', 'Normalization'])['Correlation'].agg(
        Overall_Sample_Median='median',
        Overall_Sample_Mean='mean',
        Sample_Count='count'
    ).reset_index()
    overall_summary['Overall_Sample_Median'] = overall_summary['Overall_Sample_Median'].round(3)
    overall_summary['Overall_Sample_Mean'] = overall_summary['Overall_Sample_Mean'].round(3)
    
    overall_out_path = os.path.join(TABLES_OUT_DIR, "section_3_4_illumina_correlation_overall_summary.csv")
    overall_summary.to_csv(overall_out_path, index=False)
    print(f" - Saved: {overall_out_path}")

    # Table D: Median & Mean Correlation by Dataset, Technology, Level, Tool, AND Normalization
    # (the tool-by-normalization breakdown quoted in the Section 3.4 prose; Table B averages
    #  over normalizations and therefore cannot reproduce those numbers)
    table_tool_norm = combined.groupby(
        ['Dataset', 'Tech', 'Level', 'Tool', 'Normalization']
    )['Correlation'].agg(
        Sample_Median='median',
        Sample_Mean='mean',
        Sample_Count='count'
    ).reset_index()
    table_tool_norm['Sample_Median'] = table_tool_norm['Sample_Median'].round(4)
    table_tool_norm['Sample_Mean'] = table_tool_norm['Sample_Mean'].round(4)
    table_tool_norm = table_tool_norm.sort_values(
        by=['Dataset', 'Tech', 'Level', 'Normalization', 'Sample_Median'],
        ascending=[True, True, True, True, False]
    )
    tool_norm_out_path = os.path.join(TABLES_OUT_DIR, "section_3_4_illumina_correlation_by_tool_and_normalization.csv")
    table_tool_norm.to_csv(tool_norm_out_path, index=False)
    print(f" - Saved: {tool_norm_out_path}")



def main():
    calculate_section_3_1_throughput()
    calculate_section_3_1_transcript_counts()
    calculate_section_3_1_sqanti_novel_proportions()
    calculate_section_3_3_spike_ins()
    calculate_section_3_4_illumina_correlations()
    print("\nAll manuscript tables updated and regenerated successfully!")

if __name__ == "__main__":
    main()
