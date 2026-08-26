import os
import pandas as pd
import glob

def calculate_n50(lengths, counts):
    """Calculate N50 from a frequency distribution."""
    df = pd.DataFrame({'len': lengths, 'count': counts})
    df = df.sort_values('len', ascending=False)
    df['total_bp'] = df['len'] * df['count']
    total_sum = df['total_bp'].sum()
    cumsum = df['total_bp'].cumsum()
    n50_threshold = total_sum / 2
    n50 = df[cumsum >= n50_threshold]['len'].iloc[0]
    return n50

def parse_stats(file_path):
    # Extract metadata from directory structure: read_stats/Dataset/Technology/Sample.stats
    parts = file_path.split(os.sep)
    dataset = parts[-3]
    tech = parts[-2]
    name = parts[-1].replace(".stats", "")
    
    stats = {}
    rl_hist_len = []
    rl_hist_count = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("SN"):
                parts = line.split("\t")
                key = parts[1].strip(":")
                val = float(parts[2])
                stats[key] = val
            elif line.startswith("RL"):
                parts = line.split("\t")
                rl_hist_len.append(int(parts[1]))
                rl_hist_count.append(int(parts[2]))
    
    # Extract specific metrics
    res = {
        "Dataset": dataset,
        "Technology": tech,
        "Sample": name,
        "Total_Reads": stats.get("raw total sequences", 0),
        "Total_Gb": stats.get("total length", 0) / 1e9,
        "Mapped_Reads_Pct": (stats.get("reads mapped", 0) / stats.get("raw total sequences", 1)) * 100,
        "Avg_Read_Length": stats.get("average length", 0),
        "Max_Read_Length": stats.get("maximum length", 0),
    }
    
    if rl_hist_len:
        res["Read_N50"] = calculate_n50(rl_hist_len, rl_hist_count)
    else:
        res["Read_N50"] = 0
        
    # Look for matching mapped_stats file
    mapped_stats_path = file_path.replace(".stats", ".mapped_stats")
    if os.path.exists(mapped_stats_path):
        with open(mapped_stats_path, 'r') as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    key = parts[0]
                    val = parts[1]
                    try:
                        if "." in val:
                            res[key] = float(val)
                        else:
                            res[key] = int(val)
                    except ValueError:
                        res[key] = val
        
    return res

def main():
    # Recursive search for .stats files
    stat_files = glob.glob("read_stats/**/*.stats", recursive=True)
    if not stat_files:
        print("No .stats files found in read_stats/")
        return
        
    all_results = [parse_stats(f) for f in stat_files]
    df = pd.DataFrame(all_results)
    
    # Reorder columns for readability
    cols = ["Dataset", "Technology", "Sample", "Total_Reads", "Total_Gb", "Mapped_Reads_Pct", "Avg_Read_Length", "Read_N50"]
    extra_cols = ["Average_Mapped_Read_Length", "N50_Mapped_Read_Length", "Percent_Bases_Mapped"]
    for ec in extra_cols:
        if ec in df.columns:
            cols.append(ec)
            
    df = df[cols]
    os.makedirs("analysis_results", exist_ok = True)    
    df.to_csv("analysis_results/sequencing_performance_summary.csv", index=False)
    print("Summary saved to analysis_results/sequencing_performance_summary.csv")
    print(df)

if __name__ == "__main__":
    main()
