import sys
sys.path.append('/home/yvonne/Desktop/docking_result')

from PLI import SA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path
import warnings
warnings.filterwarnings('ignore')

lig_n = 351  
frame_n = "c"+input()
root_folder = f"/home/yvonne/Desktop/docking_result/4.hybrid_dockingfda/top_ligands/drug_base/zinc{lig_n}/{frame_n}/"

h5_file = f"{root_folder}/score{lig_n}.h5"
store = pd.HDFStore(h5_file,"r")
print(store.keys())

total_score = SA.ScoreAnalysis(store)
total_score[-1].to_csv(root_folder + f"h5_{frame_n}_zinc{lig_n}.csv",index=True)

df = pd.read_csv(f"{root_folder}h5_{frame_n}_zinc{lig_n}.csv")
df = df.set_index('resname')
df = df.loc[(df !=0).any(1)]
df = df.reset_index()
df.to_csv(root_folder + f"zeroless_{frame_n}_zinc{lig_n}.csv",index=False)

#check if there's any atom clash
df2 = pd.read_csv(root_folder + f"/zeroless_{frame_n}_zinc{lig_n}.csv")
# rename column
df2 = df2.rename(columns={df2.columns[1]: f'ZINC{lig_n}_{frame_n}'}, inplace=False)
# Save DataFrame back to CSV
df2.to_csv(root_folder + f"/zeroless_{frame_n}_zinc{lig_n}.csv", index=False)

mask = df2[f'ZINC{lig_n}_{frame_n}'] > 0
if mask.any():
    print(f"clash")
else:
    print(f"safe")
    
print(df2[mask])
