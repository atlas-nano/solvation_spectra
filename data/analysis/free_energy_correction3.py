import sys
sys.path.append("/global/cfs/cdirs/m4248/xiaoxusr/solvation_scripts/python")

import argparse
import os
import matplotlib.pyplot as plt
import glob
import pandas as pd

from sea_urchin.sea_urchin import SeaUrchin
from sea_urchin.plotting.rendering import plot_structures
from free_energy_tool import EnergyCorrectionAnalyzer, load_bias_potential_data, plot_structures

import MDAnalysis as mda
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Cluster analysis for solvation structures.")
    parser.add_argument("--base_path", type=str, required=True, help="Base directory path")
    parser.add_argument("--nstrides", type=int, required=True, help="Number of strides")
    parser.add_argument("--O_radii", type=float, required=True, help="Radius for oxygen")
    parser.add_argument("--H_radii", type=float, required=True, help="Radius for hydrogen")
    parser.add_argument("--Cl_radii", type=float, required=True, help="Radius for chlorine")
    parser.add_argument("--Li_index", type=int, required=True, help="Index of lithium")
    parser.add_argument("--T", type=int, required=True, help="Temperature in Kelvin")
    parser.add_argument("--conc", type=float, required=True, help="Concentration in M/L")
    args = parser.parse_args()

    # Define parameters
    base_path = args.base_path
    nstrides = args.nstrides
    O_radii = args.O_radii
    H_radii = args.H_radii
    Cl_radii = args.Cl_radii
    Li_index = args.Li_index
    T = args.T
    conc = args.conc
    Li_id = Li_index+1

    # File paths
    lmp_file = f"{base_path}/01_IDNR/lammps.{T}K.prod.mtd.lammpstrj"
    data_file = f"{base_path}/01_IDNR/lammps.data"
    traj_list = glob.glob("*_IDNR/*lammpsdump")
    
    
    # Initialize analyzer
    analyzer = EnergyCorrectionAnalyzer(base_path, nstrides, data_file, traj_list, T)

    # Load SeaUrchin object
    obj = SeaUrchin(f"{base_path}/urchin_LiClOH_{nstrides}.pkl")

    # Load and process data
    df = pd.read_csv("./clu_analysis_sorted.csv", index_col=0).reset_index(drop=True)
    u_list = [mda.Universe(data_file, traj) for traj in traj_list]

    # Calculate free water mole fraction
    x_free_water_all_list = analyzer.calculate_free_water_fraction(u_list, distance_range=range(12,13), Li_id=Li_id, O_radii=O_radii, H_radii=H_radii)

    # Correct free energy
    df_corrected_sorted = analyzer.correct_free_energy(df, x_free_water_all_list, conc=conc)
    df_corrected_sorted.to_csv("./corrected_free_energy.csv")
    # Plot corrected free energy
    analyzer.plot_corrected_free_energy(df_corrected_sorted)

    # Visualize top clusters
    analyzer.visualize_top_clusters(df_corrected_sorted, obj)

if __name__ == "__main__":
    main()