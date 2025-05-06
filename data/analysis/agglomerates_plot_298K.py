import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import os

# Increase font size for better visualization
plt.rcParams.update({
    "font.size": 22,
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})

# Define temperatures and relevant clusters
temperatures = [298]
result_base_list = ["/global/cfs/cdirs/m4248/xiaoxusr/solvation_298K/meta_single_Li/results_bootstrap"]
cluster_combinations = ["Cl3Li", "Cl5Li2", "Cl2Li"]
concentrations = [0.5, 2, 5, 7, 10, 12, 15, 18, 20, 22, 25]
# Initialize plot
plt.figure(figsize=(10, 6))

for i in range(len(temperatures)):
    temp = temperatures[i]
    results_base = result_base_list[i]
    search_base = os.path.join(results_base, "*M")
    concentration_dirs = sorted(glob(search_base), key=lambda x: float(os.path.basename(x).replace("M", "").replace(".0", "")))

    # Extract concentration values from directory names
    concentrations = [os.path.basename(d).replace("M", "").replace(".0", "") + "M" for d in concentration_dirs]
    print(concentrations)
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
                # print(free_energy_data)
            except FileNotFoundError:
                pass
    
    # Compute means and standard deviations
    means = [np.mean(free_energy_data[conc]) if len(free_energy_data[conc]) > 0 else np.nan for conc in concentrations]
    stds = [np.std(free_energy_data[conc]) if len(free_energy_data[conc]) > 0 else np.nan for conc in concentrations]

    plt.errorbar([float(c.replace("M", "")) for c in concentrations], means, yerr=stds, fmt='o-', capsize=5, label="Cl3Li + Cl5Li2 + Cl2Li")
    print(means)
# Final plot formatting
# plt.title("Agglomerates (LiCl3, Li2Cl5, LiCl2) \n Free Energy at Different Temperatures and Concentrations")
plt.xlabel("Concentration (M)")
plt.ylabel("Free Energy (kbT)")
plt.xticks(rotation=45, ha="right")
plt.legend()
plt.tight_layout()
output_path = f"./free_energy_Agg_temp_conc.png"
plt.savefig(output_path, dpi=300)
plt.close()

print(f"Saved combined free energy plot at {output_path}.")
