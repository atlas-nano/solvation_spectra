import sys
sys.path.append("/global/cfs/cdirs/m4248/xiaoxusr/solvation_scripts/python")

import argparse
import os
import matplotlib.pyplot as plt
from sea_urchin.sea_urchin import SeaUrchin
from sea_urchin.plotting.rendering import plot_structures
from free_energy_tool import ClusterAnalyzer, get_multiple_replica_files, load_bias_potential_data, plot_structures

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Cluster analysis for solvation structures.")
    parser.add_argument("--base_path", type=str, required=True, help="Base directory path")
    parser.add_argument("--skip_frames", type=int, required=True, help="Number of frames skipped")
    parser.add_argument("--nstrides", type=int, required=True, help="Number of strides")
    parser.add_argument("--O_radii", type=float, required=True, help="Radius for oxygen")
    parser.add_argument("--Cl_radii", type=float, required=True, help="Radius for chlorine")
    parser.add_argument("--Li_index", type=int, required=True, help="Index of lithium")
    parser.add_argument("--T", type=int, required=True, help="Temperature in Kelvin")
    args = parser.parse_args()

    # Define parameters
    base_path = args.base_path
    skip_frames = args.skip_frames
    nstrides = args.nstrides
    O_radii = args.O_radii
    Cl_radii = args.Cl_radii
    Li_index = args.Li_index
    T = args.T

    # File paths
    lmp_file = f"{base_path}/01_IDNR/lammps.{T}K.prod.mtd.lammpstrj"
    bgf_file = "data.lammps"
    cluster_file = f"{base_path}/urchin_LiClOH_{nstrides}.pkl"

    if not os.path.isfile(cluster_file):
        # Cutoff distances
        cutoff_dist = {
            "O": O_radii,
            "Cl": Cl_radii,
        }

        # Reconstruction settings
        reconstruct = {
            "type": {1: "molecules"},
            "depth": 1,
            "inverse": False,
            "merge": True,
        }

        # Initialize SeaUrchin object
        obj = SeaUrchin(
            lmp_file,
            bgf_file=bgf_file,
            ref_atom="Li",
            coord_env=["O", "H", "Cl", "Li"],
            cutoff_dist=cutoff_dist,
            parallel=True,
            multiple_replica=base_path,
            skip_frames=skip_frames,
            nstrides=nstrides,
            save_name=f"urchin_LiClOH_{nstrides}.pkl",
            reconstruct=reconstruct,
            free_energy=False,
            save=True,
        )

    # Reload SeaUrchin object
    obj = SeaUrchin(f"{base_path}/urchin_LiClOH_{nstrides}.pkl")

    # Filter clusters containing the target Li atom
    clu_list = [clu for clu in obj.clusters if Li_index in clu.info["ori_idx"]]

    # Load bias potential data
    dir_path_sorted, trj_files, bias_potential_files, partial_pmf_files, full_pmf_files = get_multiple_replica_files(base_path)
    bias_pot_df = load_bias_potential_data(bias_potential_files)

    # Perform cluster analysis
    clu_analysis = ClusterAnalyzer(
        obj=obj,
        bias_pot_df=bias_pot_df,
        target_atoms_ix=[Li_index],
        geometry_pickle_data=None,
        temperature=T,
    )
    clu_analysis.calculate_total_cluster_data()
    clu_analysis.save_total_cluster_data(filename=f"total_cluster_data_{nstrides}.pkl", output_directory="./")
    clu_analysis.load_total_cluster_data(f"./total_cluster_data_{nstrides}.pkl")
    clu_analysis_df = clu_analysis.get_formula_data()

    # Sort by probability
    clu_analysis_sorted = clu_analysis_df.sort_values(["probability"], ascending=False)
    clu_analysis_sorted.to_csv("clu_analysis_sorted.csv")

    # Plot probability
    plt.figure()
    plt.plot(clu_analysis_sorted["formula"][:10], clu_analysis_sorted["probability"][:10])
    plt.xlabel("Formula")
    plt.ylabel("Probability")
    plt.xticks(rotation=90)
    plt.savefig("probability.png", bbox_inches="tight")
    plt.close()

    # Plot free energy surface
    plt.figure()
    plt.plot(clu_analysis_sorted["formula"][:5], clu_analysis_sorted["energy"][:5])
    plt.xlabel("Formula")
    plt.ylabel("Free Energy Surface (kbT)")
    plt.xticks(rotation=90)
    plt.savefig("free_energy_surface.png", bbox_inches="tight")
    plt.close()

    # Plot top cluster structures
    for i, f in enumerate(clu_analysis_sorted["formula"][:5]):
        clusters = obj.get_cluster_with_formula(f)
        try:
            plot_structures(clusters[:6])
            plt.savefig(f"formula_{i}.png", bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"Error plotting structure for formula {f}: {e}")

if __name__ == "__main__":
    main()