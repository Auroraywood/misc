#!bin/sh
files=active_files.txt
path=/home/zdx/demo/docking_pdbqt/ligand_active
find -name '*.gninatypes' > $files
sed -i 's/^/1 /g' $files

files=decoy_files.txt
path=/home/zdx/demo/docking_pdbqt/ligand_decoy
find -name '*.gninatypes' > $files
sed -i 's/^/0 /g' $files

files=receptor.txt
path=/home/zdx/demo/target
find -name '*.gninatypes' > $files
