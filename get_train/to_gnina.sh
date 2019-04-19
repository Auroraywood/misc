files=/home/zdx/cnndock/get_train/tmp/gnina_active.sh
path=/home/zdx/demo/docking_pdbqt/ligand_active
find $path -name '*.pdbqt' > $files
sed -i 's/^/gninatyper /g' $files

files=/home/zdx/cnndock/get_train/tmp/gnina_decoy.sh
path=/home/zdx/demo/docking_pdbqt/ligand_decoy
find $path -name '*.pdbqt' > $files
sed -i 's/^/gninatyper /g' $files

sh ./tmp/gnina_active.sh
sh ./tmp/gnina_decoy.sh

gninatyper /home/zdx/demo/target/*.pdbqt
