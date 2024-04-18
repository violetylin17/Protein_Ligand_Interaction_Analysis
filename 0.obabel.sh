#enter each file to split drug.pdbqt files and convert the first drug(highest score) to drug.sdf 

#!/bin/bash
path="/home/yvonne/opt/ADFRsuite-1.1dev/bin/obabel"
#name_o="zinc151"
#name_new='ZINC151'

cd $1
for i in {01..12}; do
	cd c$i
	$path $2.pdbqt -O $2_s.pdbqt -m
	$path $2_s1.pdbqt -O $2_s1.sdf
	sed -i '/>  <TORSDO>/i \\n> <Mode>\n1\n\n\> <Score>\n0\n\n\> <RMSD_LB>\n0.0\n\n\> <RMSD_UB>\n0.0\n\n\> <Iteration>\n0\n' $2_s1.sdf
	for j in {2..9}; do
		rm $2_s${j}.pdbqt
	done
	cd ..
done

# $1 = file name
# $2 = new ligand name
# $3 = old ligand name
# command: 0.obabel.sh $1 $2 $3


#單一指令:/home/yvonne/opt/ADFRsuite-1.1dev/bin/obabel ZINC000100013130_1.pdbqt -O ZINC130_s.pdbqt -m

