import os

def replace_pppm_in_file(filepath):
    with open(filepath, 'r') as file:
        content = file.read()

    if 'pppm 0.001' in content:
        content = content.replace('pppm 0.001', 'pppm 1e-5')
        with open(filepath, 'w') as file:
            file.write(content)
        print(f"Updated: {filepath}")
    else:
        print(f"No change needed: {filepath}")

def replace_pppm_in_lammps_files(root_dir):
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename == 'in.lammps' or filename == "in_colvars.lammps":
                filepath = os.path.join(dirpath, filename)
                replace_pppm_in_file(filepath)

# Replace this with your target directory path
target_directory = '/Users/user/software_dev/solvation_spectra/data'
replace_pppm_in_lammps_files(target_directory)
