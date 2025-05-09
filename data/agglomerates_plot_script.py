import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import os


# Define temperatures and relevant clusters
temperatures = [283, 298, 313]
result_base_list = ["/global/cfs/cdirs/m4248/xiaoxusr/solvation_283K/meta_single_Li/results_bootstrap2", 
"/global/cfs/cdirs/m4248/xiaoxusr/solvation_298K/meta_single_Li/results_bootstrap",
"/global/cfs/cdirs/m4248/xiaoxusr/solvation_313K/meta_single_Li/results_bootstrap2"]
cluster_combinations = ["Cl3Li", "Cl5Li2", "Cl2Li","H2Cl2LiO", "H4Cl2LiO2", "H4Cl4Li2O2", "H2Cl4Li2O", "H2Cl5Li2O"]
title="Agglomerates (LiCl3, Li2Cl5, LiCl2) \n Free Energy at Different Temperatures and Concentrations"
output_path = f"/global/cfs/cdirs/m4248/xiaoxusr/free_energy_Agg_temp_conc2.png"

def plot_agglomerates_free_energy_over_temps(
    temperatures, 
    result_base_list, 
    cluster_combinations, 
    high_conc=False,
    title="Agglomerates Free Energy Plot",
    output_path="./free_energy_Agg_temp_conc_.png",
):
    """
    conc_list: list of concentrations to plot. If False, it will extract from directory names. default: False 
    """
    # Increase font size for better visualization
    plt.rcParams.update({
        "font.size": 22,
        "axes.titlesize": 22,
        "axes.labelsize": 22,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "legend.fontsize": 20
    })

    # Initialize plot
    plt.figure(figsize=(10, 6))

    for i in range(len(temperatures)):
        temp = temperatures[i]
        results_base = result_base_list[i]
        search_base = os.path.join(results_base, "*M")
        concentration_dirs = sorted(glob(search_base), key=lambda x: float(os.path.basename(x).replace("M", "").replace(".0","")))
        
        if high_conc:
            concentration_dirs = [d for d in concentration_dirs if float(os.path.basename(d).replace("M", "").replace(".0","")) >= 15]
            
            # Extract concentration values from directory names

        concentrations = [os.path.basename(d).replace("M", "").replace(".0","") + "M" for d in concentration_dirs]
        
        free_energy_data = {conc: [] for conc in concentrations}
        
        for conc in concentrations:
            numeric_conc = conc.replace("M", "")
            for skip in range(10):
                file_path = f"{results_base}/{numeric_conc}M/SKIP_{skip}/corrected_free_energy.csv"
                try:
                    df = pd.read_csv(file_path)
                    
                    energy_combined = df[df["formula"].isin(cluster_combinations)]["energy_normalized"].values
                        
                    if len(energy_combined) > 0:
                        probability = sum(np.exp(-energy_combined))
                        energy = -np.log(probability)
                        free_energy_data[conc].append(energy)
            
                except FileNotFoundError:
                    pass
        
        # Compute means and standard deviations
        means = [np.mean(free_energy_data[conc]) if len(free_energy_data[conc]) > 0 else np.nan for conc in concentrations]
        stds = [np.std(free_energy_data[conc]) if len(free_energy_data[conc]) > 0 else np.nan for conc in concentrations]

        plt.errorbar([float(c.replace("M", "")) for c in concentrations], means, yerr=stds, fmt='o-', capsize=5, label=f"{temp}K")

    # Final plot formatting
    plt.title(title)
    plt.xlabel("Concentration (M)")
    plt.ylabel("Free Energy (kbT)")
    plt.xticks(rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Saved combined free energy plot at {output_path}.")
