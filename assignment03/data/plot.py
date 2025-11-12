import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import io

def create_combined_scaling_plots(file_list, output_dir="combined_plots"):
    """
    Creates a single strong scaling plot (MFlops vs. Ranks) per input file,
    with different problem sizes (NITER and N) plotted as separate series.
    """

    # 1. Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    generated_plots = []

    # Custom colors to ensure consistency across plots (matching the previous generated plots)
    # This is done by mapping the sorted problem size groups to fixed colors.
    custom_colors = ['tab:red', 'tab:green', 'tab:orange', 'tab:blue']

    for full_filepath in file_list:
        try:
            # --- FIX: Use a regex separator for robust reading of space-separated data ---
            df = pd.read_csv(full_filepath, sep=r'[,\s]+', engine='python')
            
            # Drop the first row if it was read as a header but contains NaN or is empty due to separators
            if df.iloc[0].isnull().all():
                df = df.iloc[1:].copy()

            required_cols = ['Ranks', 'NITER', 'N', 'MFlops']
            if not all(col in df.columns for col in required_cols):
                print(f"Skipping {full_filepath}: Missing one of the required columns: {required_cols} after reading.")
                continue

            # Ensure Ranks, NITER, N, and MFlops are clean and numeric
            for col in required_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            df.dropna(subset=required_cols, inplace=True)
            df.sort_values(by='Ranks', inplace=True)

            # 2. Group by NITER and N (fixed problem size)
            # Sort the groups to ensure a consistent color assignment (e.g., N=1000 is always the first group)
            problem_sizes_groups = df.groupby(['NITER', 'N']).apply(lambda x: x.sort_values('Ranks')).reset_index(drop=True).groupby(['NITER', 'N'])
            
            # Use os.path.basename to get just the filename for titles/output
            base_filename = os.path.basename(full_filepath)
            title_part = base_filename.replace('result_bench_', '').replace('.csv', '').replace(' ', '_').replace('.', '_')

            # 3. Create a single figure for the file
            plt.figure(figsize=(12, 8))

            # Iterate through each fixed problem size (series)
            i = 0
            for (niter, n), group in problem_sizes_groups:

                # Get the base performance (minimum Ranks in the group, which should be 1)
                min_ranks_row = group.iloc[0]
                base_ranks = min_ranks_row['Ranks']
                base_mflops = min_ranks_row['MFlops']

                # Calculate Ideal Scaling
                ideal_mflops_series = pd.Series(dtype='float64')
                if base_ranks > 0:
                    ideal_mflops_series = base_mflops * (group['Ranks'] / base_ranks)

                # Get the color for this series
                color = custom_colors[i % len(custom_colors)]

                # Plot Measured MFlops (Solid line)
                plt.plot(group['Ranks'], group['MFlops'], marker='o', linestyle='-',
                            label=f'Measured (NITER={int(niter)}, N={int(n)})',
                            color=color)

                # Plot Ideal Scaling (Dashed line)
                if not ideal_mflops_series.empty:
                    plt.plot(group['Ranks'], ideal_mflops_series, linestyle='--',
                                label=f'Ideal (NITER={int(niter)}, N={int(n)})',
                                color=color, alpha=0.6)
                
                i += 1 # Increment color index

            # 4. Finalize plot settings
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Number of Ranks (Log Scale)')
            plt.ylabel('Performance (MFlops) (Log Scale)')

            display_title = title_part.replace('_', ' ')
            title = f'Strong Scaling: {display_title}'
            plt.title(title)

            # Place legend outside the plot to prevent cluttering the data area
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.grid(True, which="both", ls="--")
            plt.tight_layout(rect=[0, 0, 0.85, 1])

            # 5. Save the plot
            plot_filename = os.path.join(
                output_dir,
                f"combined_scaling_{title_part}.png"
            )
            plt.savefig(plot_filename)
            plt.close()
            generated_plots.append(plot_filename)
            print(f"Successfully generated plot: {plot_filename}")

        except Exception as e:
            print(f"An error occurred while processing {full_filepath}: {e}")

    return generated_plots

# Example execution list, assuming a directory structure with 'blocking' and 'non-blocking' folders
files_to_process = [
    # Memdomain (Ranks 1-18)
    "./blocking/result_bench_memdomain_blocking.csv",
    "./non-blocking/result_bench_memdomain_non-blocking.csv",
    
    # Intranode (Ranks 1-72, single node)
    "./blocking/result_bench_intranode_blocking.csv",
    "./non-blocking/result_bench_intranode_non-blocking.csv",
    
    # Internode (Ranks 72-288, cluster)
    "./blocking/result_bench_internode_blocking.csv",
    "./non-blocking/result_bench_internode_non-blocking.csv",
]

# --- Main Execution ---
if __name__ == "__main__":
    created_files = create_combined_scaling_plots(files_to_process)
    print("\nPlotting complete. Generated files:", created_files)
