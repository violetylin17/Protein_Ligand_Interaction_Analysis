segments = ["TRXA","MGA"]

import sys
sys.path.append('/home/yvonne/Desktop/docking_result')

from numpy import exp
import gzip
import shutil
import subprocess
import re
import pandas as pd
import time
import os
import os.path
from multiprocessing.pool import ThreadPool
from openbabel import pybel
import tarfile
import glob
import shutil
from rdkit.Chem import PandasTools
from collections import OrderedDict
from PLI import amino_acid, analysis, interaction, license, io
import warnings
warnings.filterwarnings('ignore')

lig_f = 351
frame_n = input() ## drug base
pro_f = input()

root_folder = f"/home/yvonne/Desktop/docking_result/4.hybrid_dockingfda/top_ligands/drug_base/zinc{lig_f}/{frame_n}/"

rep_file = f"{pro_f}trxa.pdb"
lig_file =f"ZINC{lig_f}_s1.sdf"
i = "input"
o = "output"
l = "log"

i_path=os.path.join(root_folder,i)
o_path=os.path.join(root_folder,o)
l_path=os.path.join(root_folder,l)
receptor=os.path.join(root_folder,rep_file)
ligand_file=os.path.join(root_folder,lig_file)

ligands = PandasTools.LoadSDF(ligand_file, molColName=None) # Do not show Mol figure 
print(ligands)

os.makedirs(i_path, exist_ok=True)   #
os.makedirs(o_path, exist_ok = True)  # 
os.makedirs(l_path, exist_ok = True)     #
ligands = PandasTools.LoadSDF(ligand_file, molColName=None) # Do not show Mol figure 
ligandlist = []
for m in pybel.readfile("sdf", ligand_file):
        t0 = m.data["Mode"]
        pdbqt_out = os.path.join(i_path, f'{m.title}_{t0}.pdbqt')
        m.write("pdbqt", pdbqt_out, overwrite=True)
        ligandlist.append(f"{m.title}_{t0}")
ligands = ligands[["ID", "Mode", "Score", "RMSD_LB", "RMSD_UB"]]


segid = segments
prot = amino_acid.read_protein_all_atom_model(receptor, segid)

intro_message = "Date: {}   Creator: {}\n".format(time.asctime(), os.getenv('USER'))
intro_message += license.intro()


t0 = []
t1 = []
t2 = []
t3 = []
t4 = []
t5 = []
for r in prot.residues:
    segid = r.segid
    resid = r.resid
    resn = r.resn
    for a in r.atoms:
        id = a.id  ## 2022/08/16 amended by Bruce
        name = a.name 
        type = a.type
        t0.append(segid)
        t1.append(resid)
        t2.append(resn)
        t3.append(name)
        t4.append(type)
        t5.append(id)
atomlist = pd.DataFrame([t5, t0, t1, t2, t3, t4], index=["id", "segid", "resid", "resn", "name", "type"]).T

d_rid_table = dict()
for row in atomlist.iterrows():
    key = (row[1]['segid'], row[1]['resid'], row[1]['resn'], row[1]['name'])
    d_rid_table[key] =  row[1]['id']
    

ligAtomList = []
interactionList = []
for ligand in ligandlist:
    input = os.path.join(root_folder, i, f'{ligand}.pdbqt')
    output = os.path.join(root_folder, o, f'{ligand}.csv')
    log = os.path.join(root_folder, l, f'{ligand}.log')

    file = open(log, 'w', encoding='UTF-8')       #amend
    lig = amino_acid.read_atom_group(input)

    main_message = ""
    all_interaction = []

    # Calculate steric interaction
    steric_main, steric_side, msg = interaction.steric_interaction(prot, lig, log=True)
    all_interaction.append(steric_main)
    all_interaction.append(steric_side)
    main_message += msg

    # Calculate hydrophobic interaction
    hydrophobic_main, hydrophobic_side, msg = interaction.hydrophobic_interaction(prot, lig, log=True)
    all_interaction.append(hydrophobic_main)
    all_interaction.append(hydrophobic_side)
    main_message += msg

    # Calculate hbond interaction
    hbond_main, hbond_side, msg = interaction.hbond_interaction(prot, lig, log=True)
    all_interaction.append(hbond_main)
    all_interaction.append(hbond_side)
    main_message += msg

    all_interaction.reverse()
    features = OrderedDict()
    
    # Generate interaction table
    for index, types in enumerate(["steric", "hydrophobic", "hbond"]):
        for pose in interaction.list_of_interaction()[types]:
            features[pose] = all_interaction.pop()
            
    table = analysis.create_feature_table(**features)
    table.to_csv(output, float_format='%.5f', index=False)
    
    main_message += "\n\nDate: {}\n".format(time.asctime())
    file.write(intro_message + main_message)
    file.close()
    
    #
    t1 = []
    t3 = []
    t4 = []
    t5 = []

    for a in lig.atoms:
        id = int(a.index) 
        name = a.name 
        type = a.type
        t3.append(name)
        t4.append(type)
        t5.append(id)
    t1 = pd.DataFrame([t5, t3, t4], index=["id", "name", "type"]).T
    ligAtomList.append(t1)
    
    #
    st = table.filter(like="Steric").sum(axis=1)
    hy = table.filter(like="Hydrophobic").sum(axis=1)
    hb = table.filter(like="H-bond").sum(axis=1)
    surface_d = table.Surface_d
    lid = ((table["PairJ"].str.extract(r'.*[a-zA-Z]([0-9]+)')).astype(int) + 1)[0]
    t = zip(list(table["Segid"]), list(table["Resid"]), list(table["Resn"]), list(table["PairI"]))
    
    #ammemded by yizao
    rid = list()
    for x in t:
        key = (x[0], x[1], x[2], x[3])
        rid.append(d_rid_table[key])
        
        #if (x[0] == 'TRXB') & (x[1]==14) & (x[2]=='ASP') & (x[3]=='OD1'):
        #    print(d_rid_table[key])
    #rid = pd.DataFrame(rid)[0]  #####Something really wrong here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    d_result = {'lid': lid, 'rid': rid, 'st': st, 'hy': hy, 'hb': hb, 'surface_d': surface_d}
    df_final = pd.DataFrame(d_result)
    interactionList.append(df_final)
    
store = pd.HDFStore(root_folder + f"score{lig_f}.h5")
store["receptor/atoms"] = atomlist

for (index, row), m, n in zip(ligands.iterrows(), ligAtomList, interactionList):
    t = "ligands/" + str(str(row["ID"]) + "," + str(row["Mode"]))
    store[t + "/atoms"] = m
    store[t + "/interaction"] = n

store.keys()
store.close()
