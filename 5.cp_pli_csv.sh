#!/bin/bash
cd $1
mkdir $1_pli_csv
for i in {01..12}; do
	cd c$i
	cp output/*.csv ../$1_pli_csv/$1_${i}_pli.csv 
	cd ..
done

# if file name and lig name are different, set two varible: change last two $1 to $2 , than command ./.sh file_name zinc_name


#/home/yvonne/opt/ADFRsuite-1.1dev/bin/obabel ZINC000100013130_1.pdbqt -O ZINC130_s.pdbqt -m

