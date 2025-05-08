import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pycolvars import pmf
from scipy.interpolate import griddata

plt.rcParams.update({
        'font.size': 16,              # Set default font size
        'axes.labelweight': 'bold',   # Bold labels
        'axes.linewidth': 1.5,        # Axis line width
        'grid.color': '#888888',      # Gray grid lines
        'grid.linestyle': '--',       # Dashed grid lines
        'grid.linewidth': 0.5,        # Grid line width
        'xtick.major.size': 5,        # Major tick size
        'ytick.major.size': 5         # Major tick size
    })

class ColvarsAnalyzer:
    def __init__(self, base_dir, number_of_cv=1, cv_labels=None):
        """
        Initialize the ColvarsAnalyzer class.

        Parameters:
        - base_dir: The base directory containing subdirectories with simulation data.
        - number_of_cv: Number of collective variables (CVs) to analyze (1, 2, or 3).
        - cv_labels: List of labels for the CVs. Must match the number of CVs.
        """
        self.base_dir = base_dir
        self.number_of_cv = number_of_cv
        self.cv_labels = cv_labels or self._default_cv_labels()
        if len(self.cv_labels) != self.number_of_cv:
            raise ValueError("The length of cv_labels must match the number_of_cv.")

        self.directories = [os.path.join(base_dir, f'{i:02}_IDNR') for i in range(1, 11)]
        self.pmf_files = [os.path.join(dir, "colvar.out.pmf") for dir in self.directories]
        self.colvar_files = [os.path.join(dir, "colvar.out.colvars.traj") for dir in self.directories]

    def _default_cv_labels(self):
        """Provide default labels for CVs based on the number of CVs."""
        if self.number_of_cv == 1:
            return ["Li-Cl CN"]
        elif self.number_of_cv == 2:
            return ["Li-O CN", "Li-Cl CN"]
        elif self.number_of_cv == 3:
            return ["Li-O CN", "Li-Cl CN", "Li-Li CN"]
        else:
            raise ValueError("Invalid number_of_cv. Only 1, 2, or 3 CVs are supported.")

    @staticmethod
    def read_data(path):
        """Read a data file and handle any trailing null rows."""
        data = pd.read_csv(path, comment='#', sep='\s+', header=None)
        if data.iloc[-1].isnull().any():
            data = data.iloc[:-1]
        return data

    def plot_pmf(self):
        """
        Plot the Potential of Mean Force (PMF) based on the number of CVs.
        Saves the resulting plot as an image.
        """
        if self.number_of_cv == 1:
            self._plot_pmf_1cv()
        elif self.number_of_cv == 2:
            self._plot_pmf_2cv()
        elif self.number_of_cv == 3:
            self._plot_pmf_3cv()
        else:
            raise ValueError("Invalid number of CVs. Only 1, 2, or 3 are supported.")

    def _plot_pmf_1cv(self):
        """Plot PMF for 1 CV."""
        labels = self.cv_labels
        data = self.read_data(self.pmf_files[0])
        plt.plot(data[0], data[1], label="Free Energy (kBT)")
        plt.xlabel(f"{labels[0]}")
        plt.ylabel("Free Energy (kBT)")
        plt.title("Potential of Mean Force (1 CV)")
        plt.legend()
        plt.grid(True)
        plt.savefig("PMF_1CV.png", bbox_inches="tight")
        plt.close()

    def _plot_pmf_2cv(self):
        """Plot PMF for 2 CVs."""
        print(self.base_dir)
        obj = pmf.MultipleReplica(pmf_path=self.base_dir, pmf_name="colvar.out.pmf")
        obj.plot_colvars_data(con_steps=10)
        plt.title("Potential of Mean Force (2 CV)")
        plt.savefig("PMF_2CV.png", bbox_inches="tight")
        plt.close()

    def _plot_pmf_3cv(self):
        """Generate 3D PMF projections for 3 CVs and save the plot."""
        labels = self.cv_labels
        print(self.pmf_files[0])
        data = np.loadtxt(self.pmf_files[0])
        x, y, z, pmf = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

        # Assuming data is already loaded as x, y, z, and pmf
        # Generate grid for imshow
        grid_x = np.linspace(min(x), max(x), 100)
        grid_y = np.linspace(min(y), max(y), 100)
        grid_z = np.linspace(min(z), max(z), 100)

        # Interpolate PMF values onto a grid
        pmf_xy = griddata((x, y), pmf, (grid_x[None, :], grid_y[:, None]), method='cubic')
        pmf_xz = griddata((x, z), pmf, (grid_x[None, :], grid_z[:, None]), method='cubic')

        # Create subplots
        fig, axes = plt.subplots(1, 2, figsize=[20, 8])
        cmap = plt.get_cmap('viridis').copy()
        # Define color map and set vmin and vmax for consistent scaling
        vmin, vmax = np.nanmin(pmf), np.nanmax(pmf)

        # XY Plane Projection
        im_xy = axes[0].imshow(pmf_xy, extent=(min(x), max(x), min(y), max(y)), origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(im_xy, ax=axes[0])
        cbar.set_label(r'$\Delta G$ [kT]', fontsize=20, fontweight="bold")
        axes[0].set_xlabel(f"{labels[0]}", fontsize=20, fontweight="bold")
        axes[0].set_ylabel(f"{labels[1]}", fontsize=20, fontweight="bold")
        axes[0].set_title('XY Plane Projection', fontsize=20, fontweight="bold")
        axes[0].grid(True)

        # Add contour lines
        contour_levels = np.arange(0, vmax, 5)
        axes[0].contour(pmf_xy, levels=contour_levels, extent=(min(x), max(x), min(y), max(y)), colors='black', linewidths=0.75)

        # XZ Plane Projection
        im_xz = axes[1].imshow(pmf_xz, extent=(min(x), max(x), min(z), max(z)), origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(im_xz, ax=axes[1])
        cbar.set_label(r'$\Delta G$ [kT]', fontsize=20, fontweight="bold")
        axes[1].set_xlabel(f"{labels[0]}", fontsize=20, fontweight="bold")
        axes[1].set_ylabel(f"{labels[2]}", fontsize=20, fontweight="bold")
        axes[1].set_title('XZ Plane Projection', fontsize=20, fontweight="bold")
        axes[1].grid(True)
        axes[1].contour(pmf_xz, levels=contour_levels, extent=(min(x), max(x), min(z), max(z)), colors='black', linewidths=0.75)

        # Adjust layout for a better fit
        plt.tight_layout()
        plt.show()
        fig.savefig("PMF_3CV", bbox_inches="tight")

    def plot_cv_histogram_full(self):
        """
        Plot histograms of all CVs in a single figure with stacked subplots.
        """
        all_data = []
        for colvar_file in self.colvar_files:
            if os.path.isfile(colvar_file):
                all_data.append(self.read_data(colvar_file))
        
        # Concatenate all data into one DataFrame
        all_data = pd.concat(all_data)

        n = self.number_of_cv  # Number of CVs
        figs, axes = plt.subplots(n, 1, figsize=(10, n * 5))  # Create subplots
        axes = axes.flatten()

        for i, label in enumerate(self.cv_labels):
            axes[i].hist(all_data.iloc[:, i + 1], bins=15, alpha=0.7, label=label)
            axes[i].set_xlabel("CN")
            axes[i].set_ylabel("Frequency")
            axes[i].legend()
            axes[i].grid(True)

        figs.suptitle("Li-O CN and Li-Cl CN Histogram", fontsize=20)
        plt.tight_layout()
        figs.savefig("CV_Histograms.png", bbox_inches="tight")
        plt.close()
    
    def plot_traj_overtime_all(self):
        colvar_data_all = [] 
        for f in self.colvar_files:
            if os.path.isfile(f): 
                colvar_data = self.read_data(f)
                colvar_data_all.append(colvar_data)
            else:
                colvar_data_all.append(None)
        n = self.number_of_cv
        figs, axes = plt.subplots(n+1, 1, figsize=(15, 15))
        axes = axes.flatten()
        # Plot the coordination number over time
        # axes[0].plot(data[0], data[1], label='Coordination Number Li-O')
        for i in range(self.number_of_cv+1):
            if i<self.number_of_cv:
                label=self.cv_labels[i]
                for j, data in enumerate(colvar_data_all):
                    axes[i].plot(data.iloc[:,0], data.iloc[:,i+1], label=f'{label} Replica {j+1}')
                    axes[i].set_label(f'{label}')
                    axes[i].set_xlabel("Time Step (fs)")
                    axes[i].legend()
            else:
                for j, data in enumerate(colvar_data_all):
                    axes[i].plot(data.iloc[:, 0], data.iloc[:, i+1], label=f'Biased potential applied Replica {i+1}')
                    axes[i].set_xlabel("Time Step (fs)")
                    axes[i].set_ylabel("Biased Potential")
                    axes[i].legend()
            axes[i].grid(True)
        figs.suptitle('Colvars and Biased potential applied Over Time', fontsize=20)
        figs.savefig("CV_trajectory.png", bbox_inches="tight")
        plt.close()

       
                
                
