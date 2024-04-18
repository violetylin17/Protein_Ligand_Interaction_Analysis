#do vina docing

import re
from os import path
from vina import Vina
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

f_root = '/home/yvonne/Desktop/docking_result/4.hybrid_dockingfda/top_ligands/drug_base' #file that contain protein.pdbqt and ligand.pdbqt

file = input()
drug = "351"
pro_num = input()

v = Vina(sf_name='vina')    #set forcefield ad4 or vina or vinardo
v.set_receptor(f'{f_root}/zinc{drug}/{file}/{pro_num}trxa.pdbqt')  ## drug base
v.set_ligand_from_file(f'{f_root}/zinc{drug}/{file}/ZINC{drug}_s1.pdbqt') 
v.compute_vina_maps(center=[70.2, 99.99, 19.2], box_size=[30, 30, 30])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

energysr = str(energy[0])
print(energysr)


filename = f'{f_root}/zinc{drug}/{file}/ZINC{drug}_s1.sdf'

with open(filename) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if '> <Score>' in line:
        match = re.search('> <Score>', line)
        if match:
            line_str = str(lines[i+1])
            lines[i+1] = line_str.replace(line_str, energysr+' (kcal/mol)'+"\n")
            
with open(filename, 'w') as f:
    f.writelines(lines)
