
# restart_ls=(5500000 6000000)
for ((i=1; i<=10; i++)); do
    # Define the directory name
    dir_name="$(printf '%02d' $i)_IDNR"
    rep_id="$(printf '%02d' $i)"
    # restart_id=${restart_ls[$((i-1))]}
    restart_id=$((2400000 + 10#$i * 300000)) # 10#$i indicates a base 10 number

    # Check if the directory exists
    slurm_file="lammps.lammps_colvars.slurm"
    colvar_file="colvar.lmp"
    if [ -d "$dir_name" ]; then
        # Change to the directory
        cd "$dir_name"
        sed -i "s/__REPID__/$rep_id/g" "$slurm_file"
        sed -i "s/__RESTARTID__/$restart_id/g" "$slurm_file"
        sed -i "s/__REPID__/$rep_id/g" "$colvar_file"
        sed -i "s/__i__/$i/g" "$slurm_file"
        # Submit the SLURM job
        sbatch lammps.lammps_colvars.slurm
        
        # Go back to the parent directory
        cd ..
    else
        echo "Directory $dir_name does not exist!"
    fi
done