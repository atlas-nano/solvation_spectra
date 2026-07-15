import sys
sys.path.append("/global/cfs/cdirs/m4248/xiaoxusr/solvation_2024/solvation_scripts/python")

import argparse
import os
import matplotlib.pyplot as plt
import glob
import pandas as pd
import pickle

from sea_urchin import SeaUrchin
from sea_urchin.plotting.rendering import plot_structures
from free_energy_tool import EnergyCorrectionAnalyzer, load_bias_potential_data, plot_structures

import MDAnalysis as mda

# Custom unpickler to handle module path changes
class CompatibilityUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        # Map old module paths to new ones
        if module == "sea_urchin.sea_urchin":
            module = "sea_urchin"
        return super().find_class(module, name)
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Cluster analysis for solvation structures.")
    parser.add_argument("--curr_dir", type=str, required=True, help="Current working directory for analysis")
    parser.add_argument("--base_path", type=str, required=True, help="Base directory path")
    parser.add_argument("--nstrides", type=int, required=True, help="Number of strides")
    parser.add_argument("--Li_O_radii", type=float, required=True, help="Radius for lithium-oxygen")
    parser.add_argument("--Cl_H_radii", type=float, required=True, help="Radius for chlorine-hydrogen")
    parser.add_argument("--Li_Cl_radii", type=float, required=True, help="Radius for lithium-chlorine")
    parser.add_argument("--Li_index", type=int, required=True, help="Index of lithium")
    parser.add_argument("--T", type=int, required=True, help="Temperature in Kelvin")
    parser.add_argument("--conc", type=float, required=True, help="Concentration in M/L")
    parser.add_argument("--data_stride_dir", type=str, required=True, help="Directory for data stride files")
    parser.add_argument("--skip_water_calc", action="store_true", help="Skip water fraction calculation (load from cache)")
    args = parser.parse_args()

    # Define parameters
    base_path = args.base_path
    curr_dir = args.curr_dir
    nstrides = args.nstrides
    data_stride_dir = args.data_stride_dir
    Li_O_radii = args.Li_O_radii
    Cl_H_radii = args.Cl_H_radii
    # Li_Cl_radii = args.Li_Cl_radii
    Li_index = args.Li_index
    T = args.T
    conc = args.conc
    Li_id = Li_index+1

    # File paths
    lmp_file = f"{base_path}/01_IDNR/lammps.298K.prod.mtd.lammpstrj"
    data_file = f"{base_path}/01_IDNR/lammps.data"
    traj_list = glob.glob(f"{base_path}/*_IDNR/*lammpsdump")
    
    
    # Initialize analyzer
    analyzer = EnergyCorrectionAnalyzer(base_path, nstrides, data_file, traj_list, T)

    # Load SeaUrchin object with compatibility handling for module path changes
    pkl_path = f"{data_stride_dir}/urchin_LiClOH_{nstrides}.pkl"
    try:
        with open(pkl_path, 'rb') as f:
            obj = CompatibilityUnpickler(f).load()
    except Exception as e:
        print(f"Error loading pickle: {e}")
        raise

    # Load and process data
    df = pd.read_csv(f"{data_stride_dir}/clu_analysis_sorted.csv", index_col=0).reset_index(drop=True)
    u_list = [mda.Universe(data_file, traj) for traj in traj_list]

    # Calculate or load free water mole fraction (only once per concentration)
    water_cache_file = f"{curr_dir}/x_free_water_all_list.pkl"
    if args.skip_water_calc and os.path.exists(water_cache_file):
        # Load cached water fraction
        with open(water_cache_file, 'rb') as f:
            x_free_water_all_list = pickle.load(f)
        print(f"Loaded cached water fraction from {water_cache_file}")
    else:
        # Calculate water fraction and cache it
        x_free_water_all_list = analyzer.calculate_free_water_fraction(u_list, distance_range=range(12,13), Li_id=Li_id, O_radii=Li_O_radii, H_radii=Cl_H_radii)
        with open(water_cache_file, 'wb') as f:
            pickle.dump(x_free_water_all_list, f)
        print(f"Calculated and cached water fraction to {water_cache_file}")

    # Correct free energy
    df_corrected_sorted = analyzer.correct_free_energy_from_mole_fraction(df, x_free_water_all_list)
    df_corrected_sorted.to_csv("./corrected_free_energy.csv")
    # Plot corrected free energy
    analyzer.plot_corrected_free_energy(df_corrected_sorted)

    # Visualize top clusters
    analyzer.visualize_top_clusters(df_corrected_sorted, obj)

if __name__ == "__main__":
    main()