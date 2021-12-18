obabel -ipdb ${1}.pdb -omol2 -O ${1}.mol2 --partialcharges gasteiger
obabel -imol2 ${1}.mol2 -opdbqt -O ${1}.pdbqt

