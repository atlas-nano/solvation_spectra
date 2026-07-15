from sea_urchin import SeaUrchin
from sea_urchin.plotting.rendering import plot_structures

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import os
import glob
import pickle

import MDAnalysis as mda
import re

def get_multiple_replica_files(base_path,
                               prefix="colvar.out",
                               results_subdir=None):
    """
    Collect per-replica file lists from a multi-replica IDNR run.

    Parameters
    ----------
    base_path : str
        Root directory containing *_IDNR subdirectories.
    prefix : str
        Colvar file prefix. Use "colvar.fixedbias" for fixedbias runs,
        "colvar.out" (default) for plain metadynamics runs.
    results_subdir : str or None
        If set (e.g. "results"), colvar/PMF files are searched inside
        {replica_dir}/{results_subdir}/. Trajectory files are still
        searched in the replica root first, then the subdir.
    """
    trj_files = []
    partial_pmf_files = []
    bias_potential_files = []
    full_pmf_files = []

    dir_path_list = glob.glob(os.path.join(base_path, '*_IDNR'))
    dir_path_sorted = sorted(dir_path_list)

    for dir_path in dir_path_sorted:
        colvar_dir = (os.path.join(dir_path, results_subdir)
                      if results_subdir else dir_path)

        # Trajectory (.lammpstrj): check replica root first, then results subdir
        trj = glob.glob(os.path.join(dir_path, "lammps.*K.prod.*.mtd.lammpstrj"))
        if not trj and results_subdir:
            trj = glob.glob(os.path.join(colvar_dir, "lammps.*K.prod.*.mtd.lammpstrj"))
        trj = sorted(trj)

        bias_potential_file = glob.glob(
            os.path.join(colvar_dir, f"{prefix}.colvars.traj"))
        partial_pmf_file = glob.glob(
            os.path.join(colvar_dir, f"{prefix}.partial.pmf"))
        full_pmf_file = glob.glob(
            os.path.join(colvar_dir, f"{prefix}.pmf"))

        trj_files.append(trj)
        bias_potential_files.append(bias_potential_file)
        partial_pmf_files.append(partial_pmf_file)
        full_pmf_files.append(full_pmf_file)

    return dir_path_sorted, trj_files, bias_potential_files, partial_pmf_files, full_pmf_files


def build_colvar_va_table(meta_path, n_replicas=10, colvar_filename="colvar.out.colvars.traj",
                          results_subdir=None, max_step=None):
    """
    Build a V_a(t) lookup table from the full plain-MTD COLVAR trajectories.

    For multiple-walker metadynamics all replicas share one bias, so at each
    timestep t the spatial average is estimated by averaging V across all
    replicas that have data at t.  V_a(t) is the cumulative mean of those
    per-step averages from step 0 up to t.

    Parameters
    ----------
    meta_path : str
        Root directory containing XX_IDNR subdirectories.
    n_replicas : int
        Number of replicas (default 10).
    colvar_filename : str
        Name of the COLVAR trajectory file inside each replica directory.
    results_subdir : str or None
        If set, look for the COLVAR file inside {replica_dir}/{results_subdir}/.
    max_step : int or None
        Truncate all trajectories at this step (optional).

    Returns
    -------
    dict {int -> float}
        Maps integer timestep → V_a (kcal/mol).  Only steps present in the
        COLVAR data are included; use nearest-earlier lookup for gaps.
    """
    col_names = ["step", "coord_Li_O", "coord_Li_Cl", "coord_Li", "E_metadyn_d1"]
    frames = []
    for i in range(1, n_replicas + 1):
        rep_dir = os.path.join(meta_path, f"{i:02d}_IDNR")
        if results_subdir:
            rep_dir = os.path.join(rep_dir, results_subdir)
        fpath = os.path.join(rep_dir, colvar_filename)
        if not os.path.exists(fpath):
            print(f"  Warning: {fpath} not found, skipping replica {i}")
            continue
        try:
            df = pd.read_csv(fpath, sep=r"\s+", comment="#", header=None,
                             names=col_names, usecols=["step", "E_metadyn_d1"])
            df["step"] = df["step"].astype(np.int64)
            if max_step is not None:
                df = df[df["step"] <= max_step]
            frames.append(df)
        except Exception as e:
            print(f"  Warning: could not load {fpath}: {e}")

    if not frames:
        raise RuntimeError("No COLVAR files loaded — check meta_path and colvar_filename.")

    full_df = pd.concat(frames, ignore_index=True)

    # At each step, average V across all replicas that have data there
    step_mean = (
        full_df.groupby("step")["E_metadyn_d1"]
        .mean()
        .sort_index()
        .reset_index()
    )
    step_mean.columns = ["step", "V_mean"]

    # Cumulative mean over time (V_a grows as bias accumulates)
    step_mean["V_a"] = step_mean["V_mean"].expanding().mean()

    va_table = dict(zip(step_mean["step"].astype(int), step_mean["V_a"]))
    print(f"  build_colvar_va_table: {len(frames)} replicas, "
          f"{len(step_mean):,} unique steps, "
          f"step range [{step_mean['step'].iloc[0]}, {step_mean['step'].iloc[-1]}]")
    return va_table


def load_bias_potential_data(bias_potential_files, extra_files=None):
    """
    Load bias potential data from multiple replica files into a single DataFrame.

    Parameters
    ----------
    bias_potential_files : list of list of str
        One inner list per replica (output of get_multiple_replica_files).
    extra_files : list of list of str or None
        Optional extra traj files per replica (e.g. restart run files) that
        will be concatenated after the main file for each replica.
    """
    col_names_4 = ["step", "coord_Li_O", "coord_Li_Cl", "E_metadyn_d1"]
    col_names_5 = ["step", "coord_Li_O", "coord_Li_Cl", "coord_Li", "E_metadyn_d1"]

    # Detect column count from first available file
    first_file = next((f[0] for f in bias_potential_files if f), None)
    if first_file is None:
        raise FileNotFoundError("No bias potential files found.")
    df0 = pd.read_csv(first_file, sep=r"\s+", comment="#", nrows=2)
    col_names = col_names_4 if len(df0.columns) == 4 else col_names_5

    frames = []
    for i, file_list in enumerate(bias_potential_files):
        replica_frames = []
        files_to_load = list(file_list)
        if extra_files and extra_files[i]:
            files_to_load += extra_files[i]
        for fpath in files_to_load:
            try:
                df = pd.read_csv(fpath, sep=r"\s+", comment="#", names=col_names)
                df["replica_id"] = i
                replica_frames.append(df)
            except Exception as e:
                print(f"Warning: could not load {fpath}: {e}")
        if replica_frames:
            frames.append(pd.concat(replica_frames, ignore_index=True))

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()



class ClusterAnalyzer:
    def __init__(self, 
                 obj, 
                 bias_pot_df, 
                 target_atoms_ix=None,
                 formula_list=None, 
                 num_formulas=None, 
                 geometry_pickle_data=None, 
                 temperature=298,
                 filtered_clusters=False
                ):
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
        self.filtered_clusters=filtered_clusters # only analyze these clusters 

        # if there is filtered_clusters, only use these clusters as dataset
        clusters_dataset = obj.clusters if not self.filtered_clusters else self.filtered_clusters
        # if there is target_atoms_ix, only use these clusters within dataset
        if self.target_atoms_ix:
            target_atoms_clusters = []
            for clu in clusters_dataset:
                if any(idx in clu.info["ori_idx"] for idx in self.target_atoms_ix):
                    target_atoms_clusters.append(clu)
        
        self.clusters = clusters_dataset if target_atoms_ix is None else target_atoms_clusters

        
            
    ### Methods for total_cluster_data Handling ###
    
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
        # Drop frames where colvar data was not recorded (beyond end of colvar file)
        valid = np.array([v for v in bias_potential_list if not np.isnan(v)])
        if len(valid) == 0:
            return 0.0
        biased_probabilities = np.ones(len(valid)) / count_total_clusters
        bias_potential = valid * self.convert_kcal_per_mol_kbT  # Convert from kcal/mol to kbT
        logsumexp_trick = bias_potential + np.log(biased_probabilities)

        total_sum_bias = self.get_total_sum_bias()

        return np.exp(logsumexp(logsumexp_trick) - total_sum_bias)


    def get_total_sum_bias(self):
        logsumexp_trick_list = []
        count_total_clusters = len(self.clusters)

        for idx, formula in enumerate(self.total_cluster_data["clu_type"]):
            raw = self.total_cluster_data.iloc[idx]["bias_potential"]
            # Drop NaN frames (dump steps beyond colvar coverage)
            valid = np.array([v for v in raw if not np.isnan(v)])
            if len(valid) == 0:
                continue
            bias_potential = valid * self.convert_kcal_per_mol_kbT
            biased_probabilities = np.ones(len(valid)) / count_total_clusters
            logsumexp_trick = bias_potential + np.log(biased_probabilities)
            logsumexp_trick_list.append(logsumexp_trick)
        return logsumexp(np.concatenate(logsumexp_trick_list))


    def get_bias_potential_for_step(self, replica_id, timestep, cluster_idx, label):
        """Get bias potential for a specific replica and timestep."""
        try:
            bias_potential = self.bias_pot_df[
                (self.bias_pot_df["step"] == timestep) &
                (self.bias_pot_df["replica_id"] == replica_id)
            ]["E_metadyn_d1"].values[0]
        except (IndexError, KeyError):
            # Dump frame falls beyond the colvar file's last recorded step
            # (SLURM wall-time kills COLVARS output before the final dump flush).
            # Return NaN so this frame is excluded from reweighting.
            bias_potential = np.nan
        return bias_potential
    
    def get_relative_energy_from_prob(self, prob):
        """Calculate relative free energy from probabilities."""
        return -np.log(prob)

    def get_formula_data_global_reweight(self, window_steps=None, equil_steps=0, max_step=None, va_table=None):
        """
        Global reweighting with cumulative time-dependent bias correction.

        For multiple-walker metadynamics all replicas share the same bias, so
        V_a(t) is estimated as the cross-replica cumulative mean of V up to t.

        If va_table is provided (built from the full COLVAR via build_colvar_va_table),
        it is used directly for V_a lookup at each frame's step.  Otherwise V_a is
        estimated from the subsampled cluster frames only (fallback).

        Each frame's weight: w_i = exp(β * (V_i - V_a(t_i)))
        All frames pooled:   p(formula) = Σ_{i ∈ formula} w_i / Σ_i w_i

        Parameters
        ----------
        window_steps : int, optional
            Unused — kept for backward compatibility.
        equil_steps : int
            Discard frames with step < equil_steps.
        max_step : int, optional
            Discard frames with step > max_step (used for convergence checks).
        va_table : dict {int -> float}, optional
            Pre-computed V_a(t) from the full COLVAR trajectory.  Keys are
            timestep integers; values are the cumulative cross-replica mean bias
            up to that step.  Build with build_colvar_va_table().

        Returns
        -------
        DataFrame
            Columns: formula, probability, energy  (sorted by energy).
        """
        if not hasattr(self, "total_cluster_data") or self.total_cluster_data is None:
            raise RuntimeError("Call calculate_total_cluster_data() first.")

        # --- build flat (formula, replica_id, step, V) table ---
        records = []
        for idx in range(len(self.total_cluster_data)):
            row     = self.total_cluster_data.iloc[idx]
            formula = row["clu_type"]
            for r, s, v in zip(row["replica"], row["timestep"], row["bias_potential"]):
                if not np.isnan(float(v)) and int(s) >= equil_steps:
                    records.append({"formula": formula,
                                    "replica_id": int(r),
                                    "step":       int(s),
                                    "V":          float(v)})
        if not records:
            print("Warning: no valid frames in total_cluster_data.")
            return pd.DataFrame(columns=["formula", "probability", "energy"])

        flat_df = pd.DataFrame(records)
        if max_step is None:
            max_step = int(flat_df["step"].max())
        else:
            flat_df = flat_df[flat_df["step"] <= max_step]
            if flat_df.empty:
                print("Warning: no frames after applying max_step filter.")
                return pd.DataFrame(columns=["formula", "probability", "energy"])

        # --- assign V_a for each frame ---
        if va_table is not None:
            # Use pre-computed dense COLVAR V_a — look up by step, fall back to
            # nearest earlier step for any gaps.
            va_steps = np.array(sorted(va_table.keys()), dtype=np.int64)
            def _lookup(step):
                i = np.searchsorted(va_steps, step, side="right") - 1
                return va_table[int(va_steps[max(i, 0)])]
            flat_df["V_a"] = flat_df["step"].map(_lookup)
            va_source = "full-COLVAR"
        else:
            # Fallback: cross-replica cumulative mean from subsampled frames.
            # Multiple-walker MTD shares one bias, so average across all replicas.
            flat_df = flat_df.sort_values("step").reset_index(drop=True)
            flat_df["V_a"] = flat_df["V"].expanding().mean()
            # broadcast: frames at the same step share the last V_a at that step
            va_by_step = flat_df.groupby("step")["V_a"].last()
            flat_df["V_a"] = flat_df["step"].map(va_by_step)
            va_source = "subsampled-fallback"

        # --- compute weights and pool all frames ---
        beta_dV      = (flat_df["V"].values - flat_df["V_a"].values) * self.convert_kcal_per_mol_kbT
        flat_df["w"] = np.exp(beta_dV)
        Z            = flat_df["w"].sum()

        grp = (
            flat_df.groupby("formula")["w"]
            .sum()
            .div(Z)
            .rename("probability")
            .reset_index()
        )
        grp["energy"] = grp["probability"].apply(
            lambda p: -np.log(p) if p > 0 else np.inf
        )
        grp = grp.sort_values("energy").reset_index(drop=True)

        n_steps = flat_df["step"].nunique()
        print(f"  {va_source} | {n_steps} unique steps | {len(flat_df):,} frames | "
              f"top={grp.iloc[0]['formula']}  E={grp.iloc[0]['energy']:.3f} kbT")
        return grp

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
        # activity = activity_coefficient * mole_fraction 
        T = self.T # activity coefficient is calculated based on the concentration of LiCl in water, and the relationship is fitted based on experimental data. The activity coefficient is used to correct the free energy calculation by accounting for the non-ideal behavior of the solution at different concentrations.
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
    
    def correct_free_energy(self, 
                                          df, 
                                          x_free_water_all_list,
                                          conc):
        """Apply a correction to the free energy based on the local free water fraction.
        original method
        """
        x_bulk = 1
        activity = self.get_activity_from_conc(conc=conc)
        print(f"{activity=}")
        print(f"free_water_mole_fraction={np.mean(x_free_water_all_list)}")
        delta_mu = np.log((np.mean(x_free_water_all_list) * activity) / x_bulk) 
        # this should be corrected
        # delta_mu = np.log(
            # (np.mean(x_free_water_all_list) * water_activity_coefficient_simulated)
            # / (x_bulk*water_activity_coefficient_experimental)
            # )
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
            
    def correct_free_energy_from_activity_coefficient(self, 
                                          df, 
                                          x_free_water_all_list, 
                                          water_activity_coefficient_simulated, 
                                          conc):
        """Apply a correction to the free energy based on the local free water fraction.
        Apply when you know simulated activity coefficient, rather than just experimental activity coefficient. 
        delta_mu = kT*ln(
                        (water_activity_coefficient_simulated*free_water_mole_fraction)/
                        (water_activity_coefficient_experimental*local_bulk_water_mole_fraction)
                    )
        local_bulk_water_mole_fraction is assumed to be 1, so the delta_mu is simplified to 
        kT*ln(water_activity_coefficient_simulated*free_water_mole_fraction) - kT*ln(water_activity_coefficient_experimental)
        
        experimentally measured water activity is water activity coefficient * water mole fraction (average)
        the water activity coefficient at saturation is calculated and then fitted with linear regression with 1 at ideal solution
        """
        
        x_bulk = 1
        water_activity_coefficient_experimental = self.get_activity_from_conc(conc=conc)
        print(f"{water_activity_coefficient_experimental=}")
        print(f"free_water_mole_fraction={np.mean(x_free_water_all_list)}")
        delta_mu = np.log(
            (np.mean(x_free_water_all_list) * water_activity_coefficient_simulated)
            / (x_bulk*water_activity_coefficient_experimental)
            )
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

    def correct_free_energy_from_mole_fraction(self, 
                                               df, 
                                               x_free_water_all_list):
        """Apply a correction to the free energy based on the local free water fraction.
        Apply when you only know experimental activity coefficient, and you assume the simulated activity coefficient is the same.
        delta_mu = kT*ln(
                        (water_activity_coefficient_simulated*free_water_mole_fraction)/
                        (water_activity_coefficient_experimental*local_bulk_water_mole_fraction)
                    )
        local_bulk_water_mole_fraction is assumed to be 1, and water_activity_coefficient_simulated is assumed to be the same as water_activity_coefficient_experimental, 
        so the delta_mu is simplified to 
        kT*ln(free_water_mole_fraction)
        """
        x_bulk = 1 # local bulk water mole fraction is assumed to be 1
        print(f"free_water_mole_fraction={np.mean(x_free_water_all_list)}")
        delta_mu = np.log(np.mean(x_free_water_all_list) / x_bulk)
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

