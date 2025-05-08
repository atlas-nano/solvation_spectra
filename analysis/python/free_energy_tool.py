from sea_urchin.sea_urchin import SeaUrchin
from sea_urchin.plotting.rendering import plot_structures

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import os
import glob
import pickle
import re

def get_multiple_replica_files(base_path):
    # Regex patterns to extract the number from the filenames
    # partial_pattern = re.compile(r'colvar\.out\.partial\.(\d+)\.pmf')
    # full_pattern = re.compile(r'colvar\.out\.(\d+)\.pmf')
    trj_files = []
    partial_pmf_files = []
    bias_potential_files = []
    full_pmf_files=[]
    
    dir_path_list = glob.glob(os.path.join(base_path, '*_IDNR'))
    dir_path_sorted = sorted(dir_path_list)
    # Iterate through each *_IDNR directory
    for dir_path in dir_path_sorted:
        # Separate matching files into partial and full categories
        trj = glob.glob(os.path.join(dir_path, "lammps.*K.prod.*.mtd.lammp*"))
        trj = sorted(trj)
        bias_potential_file = glob.glob(os.path.join(dir_path, "colvar.out.colvars.traj"))
        partial_pmf_file = glob.glob(os.path.join(dir_path, 'colvar.out.partial.pmf'))
        full_pmf_file = glob.glob(os.path.join(dir_path, 'colvar.out.pmf'))

            
        trj_files.append(trj)
        bias_potential_files.append(bias_potential_file)
        partial_pmf_files.append(partial_pmf_file)
        full_pmf_files.append(full_pmf_file)
    # Find the latest file for both categories
    return dir_path_sorted, trj_files, bias_potential_files, partial_pmf_files, full_pmf_files 

def load_bias_potential_data(bias_potential_files):
    """Load bias potential data from multiple files into a single DataFrame."""
    bias_pot_df = pd.DataFrame()
    df0 = pd.read_csv(bias_potential_files[0][0], sep="\s+", comment="#")
    print(df0.columns)
    if len(df0.columns)==4:
        col_names = ["step", "coord_Li_O", "coord_Li_Cl", "E_metadyn_d1"]
    else:
        col_names = ["step", "coord_Li_O", "coord_Li_Cl", "coord_Li", "E_metadyn_d1"]

    for i, file in enumerate(bias_potential_files):
        df = pd.read_csv(file[0], sep="\s+", comment="#", names=col_names)
        df["replica_id"]=i
        bias_pot_df = pd.concat([bias_pot_df, df], ignore_index=True)
    return bias_pot_df



class ClusterAnalyzer:
    def __init__(self, 
                 obj, 
                 bias_pot_df, 
                 target_atoms_ix=None,
                 formula_list=None, 
                 num_formulas=None, 
                 geometry_pickle_data=None, 
                 temperature=298):
        """
        Initialize the ClusterAnalysis class.
        
        Parameters:
        - obj: An object containing the cluster data.
        - bias_pot_df: A DataFrame containing the bias potential data.
        - : all the cluster data including the replica, timestep and biased potential applied
        - target_atom_ix: If provided a list of index of atom, the clusters including these atoms will be analyzed. e.g. When we control 1 Li+ in metadynamics.
        - geometry_pickle_data: Pickle data for geometry-level cluster analysis.
        - formula_list: A list of formulas to analyze.
        - num_formulas: Number of formulas to analyze (if no list is provided).
        - temperature: Temperature of the system.
        """
    
        self.obj = obj
        self.bias_pot_df = bias_pot_df
        self.target_atoms_ix = target_atoms_ix
        self.formula_list = formula_list
        self.geometry_pickle_data = geometry_pickle_data
        self.temperature = temperature
        self.num_formulas = num_formulas
        self.k_b = 1.380649e-23  # Boltzmann constant (in kbT units)
        self.kbT = self.k_b * self.temperature
        self.convert_kcal_per_mol_J = 4184 / 6.02e23
        self.convert_kcal_per_mol_kbT = self.convert_kcal_per_mol_J / self.kbT
        if self.target_atoms_ix:
            target_atoms_clusters = []
            for clu in obj.clusters:
                if any(idx in clu.info["ori_idx"] for idx in self.target_atoms_ix):
                    target_atoms_clusters.append(clu)
        self.clusters = obj.clusters if target_atoms_ix is None else target_atoms_clusters
    
    def calculate_total_cluster_data(self):
        """Calculate the total_cluster_data by extracting bias potentials for all formulas."""
        total_cluster_data = {
            "clu_type": [], 
            "replica": [], 
            "timestep": [], 
            "bias_potential": []
        }
        
        formulas = self.get_all_formulas()
        
        for formula in formulas:
            replica_id_list, timestep_list, bias_potential_list = self.extract_biased_potential_for_formula(formula)
            total_cluster_data["clu_type"].append(formula)
            total_cluster_data["replica"].append(replica_id_list)
            total_cluster_data["timestep"].append(timestep_list)
            total_cluster_data["bias_potential"].append(bias_potential_list)
        
        self.total_cluster_data = pd.DataFrame(total_cluster_data)
        return self.total_cluster_data

    def save_total_cluster_data(self, filename, output_directory="./"):
        """Save the total_cluster_data to a pickle file."""
        with open(os.path.join(output_directory, filename), 'wb') as f:
            pickle.dump(self.total_cluster_data, f)
        print(f"Total cluster data saved to {self.total_cluster_data}.")

    def load_total_cluster_data(self, path_to_file):
        """Load the total_cluster_data from a pickle file."""
        with open(path_to_file, 'rb') as f:
            self.total_cluster_data = pickle.load(f)
        print(f"Total cluster data loaded from {path_to_file}.")
        return self.total_cluster_data
        
    ### Method to get all formulas ###
    
    def get_all_formulas(self):
        """Return all possible formulas from the object."""
        return list(set([cluster.get_chemical_formula() for cluster in self.clusters]))
   
    ### Case1: Geometry-level Analysis ###
    def get_geometry_data(self):
        """Perform geometry-level analysis if the geometry pickle data is provided."""
        
        def extract_biased_potential_for_geometry(geometry_label):
            """Extract bias potential for a specific geometry from the bias potential DataFrame."""
            replica_id_list, timestep_list, bias_potential_list = [], [], []
            clusters = [clu for idx, clu in enumerate(geometry_pickle_data.clusters) if geometry_pickle_data.labels[idx] == geometry_label]
            formula = clusters[0].get_chemical_formula()

            for idx, clu in enumerate(clusters):
                timestep = clu.info["timestamp"]
                replica_id = clu.info["replica_id"]

                if timestep is not None and replica_id is not None:
                    timestep_list.append(timestep)
                    replica_id_list.append(replica_id)

                    bias_potential = self.get_bias_potential_for_step(replica_id, timestep, idx, formula)
                    bias_potential_list.append(bias_potential)
            return replica_id_list, timestep_list, bias_potential_list        
        
        if self.geometry_pickle_data is None:
            raise ValueError("Geometry pickle data is not provided.")
            
        geometry_data = {"geometry_label": [], 
                         "replica": [], 
                         "timestep": [], 
                         "bias_potential": [], 
                         "probability": [], 
                         "energy": []}
        geometry_labels = list(set(self.geometry_pickle_data.labels))
        
        for geometry_label in geometry_labels:
            replica_id_list, timestep_list, bias_potential_list = extract_biased_potential_for_geometry(geometry_label)
            probability_list = self.compute_probability(bias_potential_list)
            energy_list = self.get_relative_energy_from_prob(probability_list)
            
            geometry_data["geometry_label"].append(geometry_label)
            geometry_data["replica"].append(replica_id_list)
            geometry_data["timestep"].append(timestep_list)
            geometry_data["bias_potential"].append(bias_potential_list)
            geometry_data["probability"].append(probability_list)
            geometry_data["energy"].append(energy_list)

            with open(os.path.join(output_directory, geometry_data.pkl), 'wb') as f:
                pickle.dump(geometry_data, f)
            print(f"Total cluster data saved to geometry_data.pkl.")

        return pd.DataFrame(geometry_data)
    

        
    ### Case 2: Formula-level Analysis ### 
    
    def get_formula_data(self, output_directory="./"):
        """Perform formula-level analysis if only formula list or formula number is provided."""
        if self.formula_list is None and self.num_formulas is None and self.target_atoms_ix is None:
            raise ValueError("Neither a formula list nor a number of formulas is provided.")
        
        formula_data = {"formula": [], "replica": [], "timestep": [], "bias_potential": [], "probability": [], "energy": []}
        
        if self.formula_list is not None:
            desired_formulas = self.formula_list
        elif self.target_atoms_ix is not None:
            desired_formulas = self.get_all_formulas()
        else:
            cluster_sorted = sorted(self.obj.cluster_types.items(), key=lambda item: item[1], reverse=True)
            desired_formulas = list(dict(cluster_sorted).keys())[:self.num_formulas]
        
        for formula in desired_formulas:
            replica_id_list, timestep_list, bias_potential_list = self.extract_biased_potential_for_formula(formula)
            probability_list = self.compute_probability(bias_potential_list)
            energy_list = self.get_relative_energy_from_prob(probability_list)
            
            formula_data["formula"].append(formula)
            formula_data["replica"].append(replica_id_list)
            formula_data["timestep"].append(timestep_list)
            formula_data["bias_potential"].append(bias_potential_list)
            formula_data["probability"].append(probability_list)
            formula_data["energy"].append(energy_list)
        
            with open(os.path.join(output_directory, "formula_data.pkl"), 'wb') as f:
                pickle.dump(formula_data, f)
            print(f"Total cluster data saved to formula_data.pkl.")
        return pd.DataFrame(formula_data)
    
 
    def extract_biased_potential_for_formula(self, formula):
        """Extract bias potential for a specific formula from the bias potential DataFrame."""
        replica_id_list, timestep_list, bias_potential_list = [], [], []
        for idx, cluster in enumerate(self.clusters):
            if cluster.get_chemical_formula() == formula:
                timestep = cluster.info["timestamp"]
                replica_id = cluster.info["replica_id"]
            
                if timestep is not None and replica_id is not None:
                    timestep_list.append(timestep)
                    replica_id_list.append(replica_id)

                    bias_potential = self.get_bias_potential_for_step(replica_id, timestep, idx, formula)
                    bias_potential_list.append(bias_potential)

        return replica_id_list, timestep_list, bias_potential_list
    
    ### Common Methods for Both Cases ###    
    
    def compute_probability(self, bias_potential_list):
        count_total_clusters = len(self.clusters)
        # print(max(bias_potential_list))
        biased_probabilities = np.ones(len(bias_potential_list)) / count_total_clusters  
        # Apply bias potential correction
        bias_potential = np.array(bias_potential_list) * self.convert_kcal_per_mol_kbT # Convert from kcal/mol to kbT
        logsumexp_trick = bias_potential+np.log(biased_probabilities) # q = p*e^(bias_potential/kbT), P = sum(q)/Q

        total_sum_bias = self.get_total_sum_bias()

        return np.exp(logsumexp(logsumexp_trick)-total_sum_bias)

    
    def get_total_sum_bias(self):
        # total_sum_bias = 0
        logsumexp_trick_list = []
        count_total_clusters = len(self.clusters)
        
        for idx, formula in enumerate(self.total_cluster_data["clu_type"]):
            bias_potential_list = np.array(self.total_cluster_data.iloc[idx]["bias_potential"])* self.convert_kcal_per_mol_kbT # convert to J/mol
            biased_probabilities = np.ones(len(bias_potential_list)) / count_total_clusters
            bias_potential = np.array(bias_potential_list)  # Convert from kcal/mol to kbT
            logsumexp_trick = bias_potential+np.log(biased_probabilities)
            logsumexp_trick_list.append(logsumexp_trick)
        return logsumexp(np.concatenate(logsumexp_trick_list)) # return a logsumexp value = log(sum(P*e^(beta*bias_potential)))

            
    def get_bias_potential_for_step(self, replica_id, timestep, cluster_idx, label):
        """Get bias potential for a specific replica and timestep."""
        try:
            bias_potential = self.bias_pot_df[
                (self.bias_pot_df["step"] == timestep) & 
                (self.bias_pot_df["replica_id"] == replica_id)
            ]["E_metadyn_d1"].values[0]
        except (IndexError, KeyError):
            print(f"No bias potential recorded for {label} at replica {replica_id}, timestep {timestep}")
            bias_potential = 0
        return bias_potential
    
    def get_relative_energy_from_prob(self, prob):
        """Calculate relative free energy from probabilities."""
        return -np.log(prob)

def logsumexp(x):
    # first reduce max value c among all number
    # take exponential, sum and log, 
    # and eventually add c to all value
    c = x.max()
    return c + np.log(np.sum(np.exp(x-c)))


class EnergyCorrectionAnalyzer():
    def __init__(self, base_path, nstrides, data_file, traj_list, T):
        self.base_path = base_path
        self.nstrides = nstrides
        self.data_file = data_file
        self.traj_list = traj_list
        self.T = T

    @staticmethod
    def count_oxygen_atoms(formula):
        """Count the number of oxygen atoms in a chemical formula."""
        oxygen_matches = re.findall(r'O(\d*)', formula)
        oxygen_count = 0
        for match in oxygen_matches:
            oxygen_count += int(match) if match else 1
        return oxygen_count

    def calculate_free_water_fraction(self, u_list, distance_range=range(12, 13), Li_id=2947, O_radii=2.65, H_radii=2.95):
        """Calculate the local free water mole fraction."""
        x_free_water_all_list = []
        for x in tqdm(distance_range):
            x_free_water_mean_list = []
            for u in u_list:
                x_free_water_list = []
                water_selection = f"byres ((around {x} (id {Li_id})) and ((type 1) or (type 2)))"
                salt_selection = f"byres ((around {x} (id {Li_id})) and ((type 3) or (type 4)))"
                non_free_water_selection = (
                    f"byres ((around {x} (id {Li_id})) and (((type 1) and (around {O_radii} (type 3))) or ((type 2) and (around {H_radii} (type 4)))))"
                )
                for ts in u.trajectory[::20]:
                    water_atoms = u.select_atoms(water_selection)
                    salt_atoms = u.select_atoms(salt_selection)
                    non_free_water_atoms = u.select_atoms(non_free_water_selection)
                    non_free_water = len(non_free_water_atoms) / 3
                    water = len(water_atoms) / 3
                    salt = len(salt_atoms) / 2
                    try:
                        x_free_water = (water - non_free_water) / (water + salt)
                        x_free_water_list.append(x_free_water)
                    except: 
                        pass
                
                x_free_water_mean_list.append(np.mean(x_free_water_list))
                print(f"{water=}, {non_free_water=}, {salt=}")
            x_free_water_all_list.append(np.mean(x_free_water_mean_list))

        return x_free_water_all_list
    
    def get_activity_from_conc(self, conc):
        T = self.T
        def get_activity(conc):
            if T==298:
                return -0.0444*conc + 1.0014
            elif T==283:
                return -0.0507*conc + 1
            elif T==313:
                return -0.0422*conc + 1
        if T==298:
            solubility = 20
        elif T==283:
            solubility = 17.5
        elif T==313: 
            solubility = 21
        if conc <= solubility:
            return get_activity(conc)
        else:
            return get_activity(conc=solubility)
        
    def correct_free_energy(self, df, x_free_water_all_list, conc):
        """Apply a correction to the free energy based on the local free water fraction."""
        x_bulk = 1
        activity = self.get_activity_from_conc(conc=conc)
        print(f"{activity=}")
        print(f"free_water_mole_fraction={x_free_water_all_list[-1]}")
        delta_mu = np.log(x_free_water_all_list[-1] * activity / x_bulk)
        df_corrected = df.copy()
        df_corrected["N_oxygen"] = [self.count_oxygen_atoms(f) for f in df_corrected["formula"]]
        df_corrected["energy_corrected"] = (
            df_corrected["energy"] - np.array(df_corrected["N_oxygen"]) * delta_mu
        )
        # Normalize energy after energy correction
        energy_corrected = df_corrected["energy_corrected"]
        probability_corrected = np.exp(-1 * energy_corrected)
        df_corrected["probability_normalized"] = probability_corrected / np.sum(probability_corrected)
        df_corrected["energy_normalized"] = -np.log(df_corrected["probability_normalized"])
        return df_corrected.sort_values(["energy_normalized"], ascending=True)     

    def plot_corrected_free_energy(self, df_corrected_sorted):
        """Plot the corrected free energy."""
        plt.plot(df_corrected_sorted["formula"][:5], df_corrected_sorted["energy_normalized"][:5])
        plt.ylabel("Corrected Free Energy (kbT)")
        plt.xlabel("Formula")
        plt.savefig("corrected_free_energy.png", bbox_inches="tight")
        plt.close()

        plt.plot(df_corrected_sorted["formula"][:5], df_corrected_sorted["probability_normalized"][:5])
        plt.ylabel("Corrected Probability")
        plt.xlabel("Formula")
        plt.savefig("corrected_probability.png", bbox_inches="tight")
        plt.close()

    def visualize_top_clusters(self, df_corrected_sorted, obj, n_top=5):
        """Visualize the top clusters based on corrected free energy."""
        for i, f in enumerate(df_corrected_sorted["formula"][:n_top]):
            clusters = obj.get_cluster_with_formula(f)
            try:
                plot_structures(clusters[:6])
                plt.savefig(f"formula_corrected_{i}.png", bbox_inches="tight")
                plt.close()
            except Exception as e:
                print(f"Error plotting structure for formula {f}: {e}")

