module load gcc/7.4.0
module load boost
ulimit -s 8192
if [ ! -z "$2" ]
then
export CUDA_VISIBLE_DEVICES=${2}
fi
./Vina-GPU --config config.txt --ligand /tmp/${1}.pdbqt --out /tmp/${1}_out.pdbqt --log /tmp/${1}_log.txt
date


# change path to vina

