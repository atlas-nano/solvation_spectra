import sys
sys.path.append("/global/cfs/cdirs/m4248/xiaoxusr/solvation_scripts/python")
import argparse
import os
from colvars_analysis_tool import ColvarsAnalyzer

# Set up the parameters for the analysis
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Analyze colvars data.")
    parser.add_argument("--base_dir", type=str, required=True, help="Base directory path for data")
    parser.add_argument("--number_of_cv", type=int, required=True, help="Number of collective variables (1, 2, or 3)")
    parser.add_argument("--cv_labels", type=str, nargs='*', default=None, help="Labels for the collective variables (optional)")

    args = parser.parse_args()

    # Extract parameters
    base_dir = args.base_dir
    print("base_dir")
    number_of_cv = args.number_of_cv
    cv_labels = args.cv_labels if args.cv_labels is not None else None

    # cv_labels = ["Li-O CN", "Li-Cl CN", "Li-Li CN"]  # Labels for the CVs

    # Initialize the analyzer
    analyzer = ColvarsAnalyzer(base_dir=base_dir, number_of_cv=number_of_cv, cv_labels=cv_labels)

    # Run various analyses
    print("Plotting PMF...")
    analyzer.plot_pmf()

    print("Plotting trajectory data over time for all replicas...")
    analyzer.plot_traj_overtime_all()

    print("Plotting histograms of CVs...")
    analyzer.plot_cv_histogram_full()

    print("CV Histogram Analysis complete! Check the output files in the current directory.")

if __name__ == "__main__":
    main()