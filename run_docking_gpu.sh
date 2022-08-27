export CUDA_VISIBLE_DEVICES=${0}
./Vina-GPU --config config.txt --ligand /tmp/${1}.pdbqt --out /tmp/${1}_out.pdbqt --log /tmp/${1}_log.txt
date


# change path to vina

