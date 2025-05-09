import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Define results path
results_base = "/global/cfs/cdirs/m4248/xiaoxusr/solvation_313K/meta_single_Li/results_bootstrap"
# Re-import necessary libraries since execution state was reset
import matplotlib.pyplot as plt
import numpy as np

# Increase font size without changing the original settings' structure
plt.rcParams.update({
    "font.size": 22,  # Increased overall font size while keeping the original style
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})


# Define concentrations
concentrations = [0.5, 1, 2, 5, 7, 10, 12, 15, 17, 19, 21, 23, 26]
cluster_labels = ["H8LiO4", "H6ClLiO3", "H4ClLiO2", "Cl3Li", "H2Cl2LiO", "H4Cl2LiO2"]  # You will dynamically update this list

for cluster in cluster_labels:
    cluster_data = {}

    if cluster == "Cl3Li":
        cluster_only = "Cl3Li"
        cluster_combined = ["Cl3Li", "Cl5Li2"]
    
    else:
        cluster_only = cluster
        cluster_combined = [cluster]

    free_energy_only = {}
    free_energy_sum = {}

    for conc in concentrations:
        free_energy_only[conc] = []
        free_energy_sum[conc] = []
        
        for skip in range(10):
            file_path = f"{results_base}/{conc}M/SKIP_{skip}/corrected_free_energy.csv"
            try:
                df = pd.read_csv(file_path)
                
                # Energy for Cl3Li alone
                energy_only = df[df["formula"] == cluster_only]["energy_normalized"].values
                if len(energy_only) > 0:
                    free_energy_only[conc].append(np.sum(energy_only))  # Sum of Cl3Li alone

                # Energy sum of Cl3Li + Cl5Li2
                energy_combined = df[df["formula"].isin(cluster_combined)]["energy_normalized"].values
                probability = 0
                for e in energy_combined:
                    probability += np.exp(-e)
                energy = -np.log(probability)
                if len(energy_combined) > 0:
                    free_energy_sum[conc].append(energy)  # Sum of both

            except FileNotFoundError:
                pass

    # Ensure at least some data is collected
    if not free_energy_only or not free_energy_sum:
        print(f"No valid data found for {cluster}, skipping.")
        continue

    # Compute means and standard deviations
    mean_free_energy_only = [np.mean(energies) if len(energies) > 0 else np.nan for energies in free_energy_only.values()]
    std_free_energy_only = [np.std(energies) if len(energies) > 0 else np.nan for energies in free_energy_only.values()]
    
    mean_free_energy_sum = [np.mean(energies) if len(energies) > 0 else np.nan for energies in free_energy_sum.values()]
    std_free_energy_sum = [np.std(energies) if len(energies) > 0 else np.nan for energies in free_energy_sum.values()]

    # Plot free energy with error bars
    plt.figure(figsize=(10, 6))

    if cluster == "Cl3Li":
        # Line for Cl3Li only
        plt.errorbar(concentrations, mean_free_energy_only, yerr=std_free_energy_only, fmt='o-', capsize=5, label="Cl3Li")
        # Line for Cl3Li + Cl5Li2
        plt.errorbar(concentrations, mean_free_energy_sum, yerr=std_free_energy_sum, fmt='s--', capsize=5, label="Cl3Li + Cl5Li2")

    else:
        plt.errorbar(concentrations, mean_free_energy_only, yerr=std_free_energy_only, fmt='o-', capsize=5, label=cluster)

    plt.xlabel("Concentrations (M)")
    plt.ylabel("Free Energy (kbT)")
    plt.xticks(rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    
    output_path = f"{results_base}/free_energy_error_plot_{cluster}.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"Saved error bar plot for {cluster} at {output_path}.")