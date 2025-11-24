import os
import glob
import pandas as pd
import argparse
import gc

__version__ = "1.0"

def get_file_attributes(filename):
    """
    Parses filename to extract attributes.
    Assumes format: SampleID.Reference.Strand.Mod...
    """
    basename = os.path.basename(filename)
    parts = basename.split('.')
    
    # 1. Extract Reference (Assumed to be second part)
    # Example: ENCBS048EXC.GRCh38.plus... -> GRCh38
    ref = parts[1] if len(parts) > 1 else "unknown"
    
    # 2. Extract Sample ID (Assumed to be the first part)
    # This becomes the column header (e.g., "Pct_ENCBS048EXC")
    sample_id = parts[0] if len(parts) > 0 else "unknown_sample"
    
    # 3. Extract Modification
    known_mods = ['m5C', 'm6A', 'Nm', 'pseU', 'inosine', '5hmCG', '6mA']
    mod = "unknown"
    for m in known_mods:
        if f".{m}." in basename or f"_{m}_" in basename:
            mod = m
            break
    
    # 4. Strand
    strand = "unknown"
    if ".plus" in basename or "_plus" in basename:
        strand = "plus"
    elif ".minus" in basename or "_minus" in basename:
        strand = "minus"
        
    return sample_id, ref, mod, strand

def process_replicates(input_folders, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if isinstance(input_folders, str):
        input_folders = [input_folders]

    # Structure: inventory[ref][mod][sample_id] = {'plus': path, 'minus': path}
    inventory = {}

    print(f"Scanning {len(input_folders)} folders...")

    # 1. Scan ALL folders and pool files together
    total_files = 0
    for folder in input_folders:
        files = glob.glob(os.path.join(folder, "*.bed"))
        total_files += len(files)
        
        for f in files:
            sample_id, ref, mod, strand = get_file_attributes(f)
            
            if mod == "unknown" or strand == "unknown":
                continue
                
            # Initialize nested structure
            if ref not in inventory: inventory[ref] = {}
            if mod not in inventory[ref]: inventory[ref][mod] = {}
            if sample_id not in inventory[ref][mod]: inventory[ref][mod][sample_id] = {}
            
            # Store file path
            inventory[ref][mod][sample_id][strand] = f

    print(f"Found {total_files} files. Grouping by Reference and Modification...")

    # 2. Process each Ref+Mod Group
    for ref in inventory:
        for mod in inventory[ref]:
            
            samples_dict = inventory[ref][mod]
            sample_list = sorted(samples_dict.keys()) # List of sample IDs (e.g. ENCBS1, ENCBS2...)
            
            if len(sample_list) < 2:
                print(f"\nSkipping {ref}-{mod}: Only found 1 sample ({sample_list[0]}). Need >1 to compare.")
                continue

            print(f"\nProcessing Group: {ref} - {mod}")
            print(f"  found {len(sample_list)} samples: {', '.join(sample_list)}")
            
            merged_df = None
            
            # 3. Merge logic
            for sample_id in sample_list:
                print(f"  Loading {sample_id}...")
                
                file_paths = samples_dict[sample_id]
                
                # Load Plus and Minus for this specific Sample ID
                sample_dfs = []
                for strand in ['plus', 'minus']:
                    if strand in file_paths:
                        fpath = file_paths[strand]
                        try:
                            # Read BED file
                            # Col 0=Chrom, 1=Start, 2=End, 5=Strand, 10=Percentage
                            df = pd.read_csv(fpath, sep='\t', header=None, comment='#', 
                                           usecols=[0, 1, 2, 5, 10],
                                           names=['Chrom', 'Start', 'End', 'Strand', f'Pct_{sample_id}'],
                                           dtype={'Chrom': str, 'Strand': str})
                            sample_dfs.append(df)
                        except Exception as e:
                            print(f"    Error reading {os.path.basename(fpath)}: {e}")
                
                if not sample_dfs:
                    continue

                # Combine Plus/Minus for this sample
                current_sample_df = pd.concat(sample_dfs, ignore_index=True)
                
                # Merge into Master Matrix
                if merged_df is None:
                    merged_df = current_sample_df
                else:
                    # Outer merge ensures we keep sites even if they only appear in one sample
                    merged_df = pd.merge(merged_df, current_sample_df, 
                                       on=['Chrom', 'Start', 'End', 'Strand'], 
                                       how='outer')

            # 4. Save Result
            if merged_df is not None:
                out_name = f"Compare_{ref}_{mod}.csv"
                out_path = os.path.join(output_folder, out_name)
                
                merged_df.to_csv(out_path, index=False)
                print(f"  -> Saved: {out_name} ({len(merged_df)} rows)")
                
                del merged_df
                gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Samples/Replicates (Folder Independent)")
    parser.add_argument("input_folders", nargs='+', help="List of folders to scan for files")
    parser.add_argument("--output", help="Output folder", default="comparison_results")
    print(f"report_bed_rep.py version {__version__}")
    
    args = parser.parse_args()
    process_replicates(args.input_folders, args.output)