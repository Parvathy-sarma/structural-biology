######RMSD — Kabsch algorithm, rotation matrix, superposition
# RMSD from scratch with NumPy; validate against MDAnalysis built-in
# NMR ensemble — multiple models in one PDB

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser, PDBList
import json
import requests
import warnings
warnings.filterwarnings('ignore')
import MDAnalysis as mda
from MDAnalysis.analysis import rms 

####The goal in here is to compare the RMSD obtained from MD analysis package, simply RMSD based on least square analysis method and RMSD of the protein after superimposing the models one over the other.
####The superimposing of one molecule over the other is done using Kabsch algorithm.


###First goal is to calculate the rotation, translation using KAbsch algorithm###
##Steps for Kabsch algorithm are as follows###

def centroid(coords):
    return coords.mean(axis = 0)

def kabsch_rotation(P,Q):
    #Step 1 center both the proteins##
    ##Centering just includes bringing the protein to one point that will be (0,0,0),
    #To bring the new centroid to be zero, all you need is to just subtract each and every point to centroid.

    P_centered = P - centroid(P)
    Q_centered = Q - centroid(Q)

    ##Step 2 - find covariance matrix - finds how one protein is positioned relative to other
    H = P_centered.T @ Q_centered

    ###Step 3: Single value Decomposition -
    # SVD finds the best coordinate systems (axes) in which:
    # P and Q are most aligned
    # their directional correlations are maximized

    U,S,Vt = np.linalg.svd(H)

    #Step 4: handle reflection (ensure proper rotation and not mirror)
    d = np.linalg.det(Vt.T@U.T)
    D = np.diag([1,1,d]) #d = +1 or -1

    ##Rotation matrix
    R = Vt.T@D@U.T
    #####
    # print("The values of kabsch rotation, for debugging")
    # print(f"{Q_centered} : Q_centered")
    # print(f"{P_centered}: P_centered ")
    # print(f"{U}, {S}, {Vt}: SVD values")
    # print(f"{d},{D}: Diagonal matrix values for handling reflection")
    # print(f"{R}, the rotation matrix for the selected models")

    return R


#####Caluculate normal RMSD without superimpositions###

def rmsd_no_superposition(P,Q):
    difference = P-Q
    return np.sqrt((difference**2).sum(axis=1).mean())

####Kabsch_rmsd###

def rmsd_kabsch(P,Q):
    
    P_centered = P - centroid(P)
    Q_centered = Q - centroid(Q)
    R = kabsch_rotation(P_centered,Q_centered)

    P_rotated = P_centered@R.T
    ##RMSD
    diff = P_rotated - Q_centered
    # print(f"{R}, the RMSD for the selected models using rmsd kabsch algorithm")
    return np.sqrt((diff**2).sum(axis =1).mean())


###Fetch protein PARP1-Zn12 and the coordinates

###Get coords of C-alpha atoms of every model in the protein###

def get_ca_coords_all_models(pdb_id):
    protein = PDBList()
    protein.retrieve_pdb_file(pdb_id, file_format = "pdb", pdir= ".",)
    parser = PDBParser()
    structure = parser.get_structure('protein', file= f"pdb{pdb_id.lower()}.ent")

    all_models = []
    for model in structure:
        coords = []
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_vector().get_array())
        if len(coords) > 0:
            all_models.append(np.array(coords))

    return all_models

models = get_ca_coords_all_models(pdb_id="2N8A")
print(f"Number of models in ensemble: {len(models)}")
print(f"Residues per model: {len(models[0])}")

###Compute pairwise RMSD using kabasch RMSD

n_models = len(models)
rmsd_matrix = np.zeros((n_models,n_models))
nosuper = np.zeros((n_models,n_models))
for i in range(n_models):
    for j in range(n_models):
        if i != j:
            rmsd_matrix[i,j] = rmsd_kabsch(models[i], models[j])
            nosuper[i,j] = rmsd_no_superposition(models[i], models[j])
            # print(f"{i}th and {j}th model evaluation completed")
            #np.savetxt(f"./rmsd_matrix{i}_{j}.txt", rmsd_matrix)



# np.savetxt('rmsd_mtrix(kabsh)_vs first.txt', rmsd_vs_first, delimiter=',')
print(f"\nPairwise RMSD statistics: ")
print(f"\n RMSD without any superposition: {nosuper.mean():.3f} Å" )
print(f"Mean RMSD between models: {rmsd_matrix[rmsd_matrix>0].mean():.3f} Å")
print(f"Max  RMSD between models: {rmsd_matrix.max():.3f} Å")
print(f"Min  RMSD between models: {rmsd_matrix[rmsd_matrix>0].min():.3f} Å")

# ════════════════════════════════════════════════════════════════
# PART 3: Validate against MDAnalysis
# ════════════════════════════════════════════════════════════════

print("\n── Validating against MDAnalysis ───────────────────────")

# Load with MDAnalysis
u = mda.Universe('pdb2n8a.ent')
print(f"MDAnalysis loaded: {len(u.trajectory)} frames")

# Select Cα atoms
ca_selection = u.select_atoms('name CA')

# Compute RMSD of all frames vs frame 0 using MDAnalysis
R = rms.RMSD(ca_selection, ref_frame=0)
R.run()
mda_rmsds = R.results.rmsd[:, 2]  # column 2 = RMSD values


# Compare our implementation vs MDAnalysis for frames 0 vs 1
u.trajectory[0]
coords_frame0 = ca_selection.positions.copy()
u.trajectory[1]
coords_frame1 = ca_selection.positions.copy()

our_rmsd = rmsd_kabsch(coords_frame0, coords_frame1)
mda_rmsd = mda_rmsds[1]


print(f"\nFrame 0 vs Frame 1:")
print(f"  Our Kabsch RMSD:    {our_rmsd:.4f} Å")
print(f"  MDAnalysis RMSD:    {mda_rmsd:.4f} Å")
print(f"  Difference:         {abs(our_rmsd - mda_rmsd):.6f} Å")

if abs(our_rmsd - mda_rmsd) < 0.001:
    print("  ✓ VALIDATION PASSED — our implementation matches MDAnalysis")
else:
    print("  ✗ Check implementation — values differ")

# ════════════════════════════════════════════════════════════════
# PART 4: Plots
# ════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Pairwise RMSD matrix
im = axes[0].imshow(rmsd_matrix, cmap='viridis')
axes[0].set_title('Pairwise RMSD Matrix\n(all NMR models)', fontsize=11)
axes[0].set_xlabel('Model index')
axes[0].set_ylabel('Model index')
plt.colorbar(im, ax=axes[0], label='RMSD (Å)')

# Plot 2: RMSD of each model vs model 1
rmsd_vs_first = rmsd_matrix[0, :]
axes[1].plot(range(n_models), rmsd_vs_first,
             'o-', color='steelblue', markersize=5)
axes[1].set_title('RMSD vs Model 1\n(conformational spread)', fontsize=11)
axes[1].set_xlabel('Model index')
axes[1].set_ylabel('RMSD (Å)')
axes[1].axhline(y=rmsd_vs_first.mean(), color='red',
                linestyle='--', label=f'Mean = {rmsd_vs_first.mean():.2f} Å')
axes[1].legend()

# Plot 3: MDAnalysis RMSD over frames
axes[2].plot(mda_rmsds, color='coral', linewidth=1.5)
axes[2].set_title('RMSD over NMR frames\n(MDAnalysis)', fontsize=11)
axes[2].set_xlabel('Frame')
axes[2].set_ylabel('RMSD vs frame 0 (Å)')
axes[2].axhline(y=mda_rmsds.mean(), color='red',
                linestyle='--', label=f'Mean = {mda_rmsds.mean():.2f} Å')
axes[2].legend()

plt.suptitle('RMSD Analysis of 2N8A NMR Ensemble', fontsize=13)
plt.tight_layout()
plt.savefig('results_today/w02d2_rmsd_analysis.png', dpi=150)
plt.show()
print("\nSaved to w02d2_rmsd_analysis.png")
