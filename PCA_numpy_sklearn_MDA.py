####PCA analysis using sckit learn and MD analysis####

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from Bio.PDB import PDBParser, PDBList
import MDAnalysis as mda
from MDAnalysis.analysis import pca as mdapca
import urllib.request
import os
import warnings
warnings.filterwarnings('ignore')

###Download PDB files and save the CA coordinates####

def dwd_ca_all_models(pdb_file, chain_id='A'):
    pdb_test = PDBList()
    pdb_test.retrieve_pdb_file(pdb_file, pdir="./", file_format='pdb')

    parser = PDBParser()
    structure = parser.get_structure('pdb_test', f"pdb{pdb_file.lower()}.ent")

    models = []

    for model in structure:
        coords = []  # 

        for chain in model:
            if chain.id != chain_id:   # 
                continue

            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())

        if len(coords) > 0:
            models.append(coords)


    if len(models) == 0:
        raise ValueError("No CA atoms found. Check chain ID or PDB file.")
   
    min_len = min(len(m) for m in models)
    coords_array = np.array([m[:min_len] for m in models])

    return coords_array

####PCA using numpy#####
def pca_numpy(coords_array):
    
    n_models,n_residues,n_dim = coords_array.shape
    ###Flattening the array, such that each model represent each rows, and each coloumn represent residues)
    X = coords_array.reshape(n_models,n_residues*n_dim)
    print(f"Data matrix shape: {X.shape}")
    print(f" {n_models} conformations, {n_residues*n_dim} features")


    ####center your entire system###
    ####just find the X-Xmean###
    X_centered = X- X.mean(axis=0)

    ###Find the covariance matrix of the matrix X_centerd###

    X_cov = np.cov(X_centered.T)
    print(f"Covariance matrix shape = {X_cov.shape}")

    ###Find eigen values and eigen vectors####
    X_eigval, X_eigvec = np.linalg.eigh(X_cov)

    ####Sort in descending eigenvalue####
    eig_sort = np.argsort(X_eigval)[::-1]
    X_eigval = X_eigval[eig_sort]
    X_eigvec = X_eigvec[:,eig_sort]

    # Explained variance ratio
    total_variance = X_eigval.sum()
    explained_variance_ratio = X_eigval / total_variance
    
    # Project data onto principal components
    # (n_models, n_components)
    projections = X_centered @ X_eigvec
    
    return {
        'projections':              projections,
        'eigenvectors':             X_eigvec,
        'eigenvalues':              X_eigval,
        'explained_variance_ratio': explained_variance_ratio,
        'mean_conformation':        X.mean(axis =0),
        'X_centred':                X_centered
    }

####Compute PCA for the PDB####
results ={}

coords = dwd_ca_all_models('2n8a', chain_id='A')
pca_result = pca_numpy(coords)

####Using sklearn PCA for validation###
n_models, n_res, _ = coords.shape
X_flat = coords.reshape(n_models, -1)
X_centered = X_flat - X_flat.mean(axis =0)

sklearn_pca = PCA(n_components=min(10,n_models-1))
sklearn_pca.fit(X_centered)

 # Compare
our_ev   = pca_result['explained_variance_ratio'][:3] * 100
sk_ev    = sklearn_pca.explained_variance_ratio_[:3] * 100
    
print(f"\nExplained variance (first 3 PCs):")
print(f"  Our PCA:    PC1={our_ev[0]:.1f}%  PC2={our_ev[1]:.1f}%  PC3={our_ev[2]:.1f}%")
print(f"  sklearn:    PC1={sk_ev[0]:.1f}%  PC2={sk_ev[1]:.1f}%  PC3={sk_ev[2]:.1f}%")
    
cumulative = np.cumsum(pca_result['explained_variance_ratio']) * 100
n_for_90   = np.argmax(cumulative >= 90) + 1
print(f"\nPCs needed to explain 90% variance: {n_for_90}")
    
results['2n8a_pdb'] = {
        'pca':        pca_result,
        'sklearn':    sklearn_pca,
        'coords':     coords,
        'cumulative': cumulative
    }

###Validate PCS with MD analysis####

u = mda.Universe('pdb2n8a.ent')
ca = u.select_atoms('name CA')
pc = mdapca.PCA(u, select='name CA', n_components=3)
pc.run()

print(f"MDAnalysis PC1 explained variance: "
      f"{pc.results.variance[0]/pc.results.variance.sum()*100:.1f}%")
print(f"Our PC1 explained variance:        "
      f"{results['2n8a_pdb']['pca']['explained_variance_ratio'][0]*100:.1f}%")


fig = plt.figure(figsize=(18, 12))
gs  = gridspec.GridSpec(2, 4, figure=fig)

for col, (pdb_id, res) in enumerate(results.items()):
    pca_r  = res['pca']
    proj   = pca_r['projections']
    evr    = pca_r['explained_variance_ratio'] * 100
    cumul  = res['cumulative']
    n_mod  = proj.shape[0]

    # ── Plot 1: Scree plot ────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, col*2])
    n_show = min(10, len(evr))
    ax1.bar(range(1, n_show+1), evr[:n_show],
            color='steelblue', alpha=0.8, label='Individual')
    ax1.plot(range(1, n_show+1), cumul[:n_show],
             'ro-', markersize=5, label='Cumulative')
    ax1.axhline(y=90, color='gray', linestyle='--',
                alpha=0.7, label='90% threshold')
    ax1.set_xlabel('Principal Component')
    ax1.set_ylabel('Explained Variance (%)')
    ax1.set_title(f'{pdb_id}\nScree Plot', fontsize=11)
    ax1.legend(fontsize=8)
    ax1.set_xticks(range(1, n_show+1))

    # ── Plot 2: PC1 vs PC2 scatter ───────────────────────────
    ax2 = fig.add_subplot(gs[0, col*2+1])
    scatter = ax2.scatter(proj[:, 0], proj[:, 1],
                          c=range(n_mod),
                          cmap='viridis', s=80, zorder=3)
    # Label each point with model number
    for i in range(n_mod):
        ax2.annotate(str(i+1), (proj[i,0], proj[i,1]),
                     fontsize=7, ha='center', va='bottom')
    plt.colorbar(scatter, ax=ax2, label='Model index')
    ax2.set_xlabel(f'PC1 ({evr[0]:.1f}%)')
    ax2.set_ylabel(f'PC2 ({evr[1]:.1f}%)')
    ax2.set_title(f'{pdb_id} — PC1 vs PC2\nConformational landscape', fontsize=11)
    ax2.grid(True, alpha=0.3)

    # ── Plot 3: Per-residue contribution to PC1 ───────────────
    ax3 = fig.add_subplot(gs[1, col*2])
    # PC1 eigenvector reshaped to (n_residues, 3)
    n_res   = res['coords'].shape[1]
    pc1_vec = pca_r['eigenvectors'][:, 0].reshape(n_res, 3)
    # Magnitude of displacement per residue
    pc1_mag = np.sqrt((pc1_vec**2).sum(axis=1))
    ax3.bar(range(1, n_res+1), pc1_mag, color='coral', alpha=0.8)
    ax3.set_xlabel('Residue index')
    ax3.set_ylabel('PC1 displacement magnitude (Å)')
    ax3.set_title(f'{pdb_id} — Per-residue PC1 contribution\n'
                  f'(flexible regions = high values)', fontsize=11)

    # ── Plot 4: PC1 vs PC2 vs PC3 (3D) ───────────────────────
    ax4 = fig.add_subplot(gs[1, col*2+1], projection='3d')
    ax4.scatter(proj[:,0], proj[:,1], proj[:,2],
                c=range(n_mod), cmap='viridis', s=60)
    ax4.set_xlabel(f'PC1 ({evr[0]:.1f}%)', fontsize=8)
    ax4.set_ylabel(f'PC2 ({evr[1]:.1f}%)', fontsize=8)
    ax4.set_zlabel(f'PC3 ({evr[2]:.1f}%)', fontsize=8)
    ax4.set_title(f'{pdb_id} — 3D conformational space', fontsize=11)

plt.suptitle('PCA of NMR Conformational Ensembles', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('w02d3_pca_conformations.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nSaved w02d3_pca_conformations.png")
