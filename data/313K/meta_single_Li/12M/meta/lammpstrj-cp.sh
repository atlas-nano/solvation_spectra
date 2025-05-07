
for ((i=1; i<=10; i++)); do
   dir_name="$(printf '%02d' $i)_IDNR"
   data1="/global/cfs/cdirs/m4248/xiaoxusr/solvation_313K/meta_single_Li/12M/classical/data.lammps"
   data2="data.lammps"
   data3="lammps.data"
   trj="lammps.313K.prod.$(printf '%02d' $i).mtd.lammps"
   trj_1="lammps.313K.prod.$(printf '%02d' $i).mtd.lammpsdump"
   trj_2="lammps.313K.prod.$(printf '%02d' $i).mtd.lammpstrj"
   trj_3="lammps.313K.prod.mtd.lammpstrj"
   # dir_name="${i}_IDNR"
   cd $dir_name
   cp $trj $trj_1
   cp $trj $trj_2
   cp $trj_1 $trj_3
   rm $trj
   cp $data1 $data2
   cp $data2 $data3
   # update parameters for each replicas
   cd ..
done

