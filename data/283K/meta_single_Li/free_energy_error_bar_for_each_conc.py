import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Define results path
results_base = "/global/cfs/cdirs/m4248/xiaoxusr/solvation_283K/meta_single_Li/results_bootstrap"
plt.rcParams.update({
    "font.size": 22,  # Increased overall font size while keeping the original style
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})
# Define concentrations
concentrations = [0.5, 1, 2, 5, 7, 10, 12, 15, 17, 19, 21, 23]


for conc in concentrations:
    free_energy_data = {}
    
    # Loop over SKIP_FRAMES
    for skip in range(10):
        file_path = f"{results_base}/{conc}M/SKIP_{skip}/corrected_free_energy.csv"
        try:
            df = pd.read_csv(file_path)
            for index, row in df.iterrows():
                cluster = row["formula"]
                energy = row["energy_normalized"]
                if cluster not in free_energy_data:
                    free_energy_data[cluster] = []
                free_energy_data[cluster].append(energy)
        except FileNotFoundError:
            pass

    # Ensure at least some data is collected
    if not free_energy_data:
        print(f"No valid data found for {conc}M, skipping.")
        continue

    # Compute mean and standard deviation
    mean_free_energy = {cluster: np.mean(energies) for cluster, energies in free_energy_data.items()}
    std_free_energy = {cluster: np.std(energies) for cluster, energies in free_energy_data.items()}

    # Rank clusters by their mean free energy
    sorted_clusters = sorted(mean_free_energy.keys(), key=lambda x: mean_free_energy[x])

    # Extract sorted values
    sorted_means = [mean_free_energy[cluster] for cluster in sorted_clusters]
    sorted_stds = [std_free_energy[cluster] for cluster in sorted_clusters]

    # Plot free energy with error bars
    plt.figure(figsize=(10, 6))
    plt.errorbar(sorted_clusters[:5], sorted_means[:5], yerr=sorted_stds[:5], fmt='o-', capsize=5, label=f"{conc}M")
    plt.xlabel("Cluster Formula")
    plt.ylabel("Free Energy (kbT)")
    plt.xticks(rotation=45, ha="right")
    plt.title(f"Free Energy vs Cluster (Conc: {conc}M)")
    plt.legend()
    plt.tight_layout()
    
    output_path = f"{results_base}/{conc}M/free_energy_error_plot.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Saved error bar plot for {conc}M concentration at {output_path}.")