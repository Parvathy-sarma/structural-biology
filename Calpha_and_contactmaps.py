#######w05:BioPython PDB parser — SMCRA hierarchy, Residue, Atom, coordinates#####

import numpy as np
import matplotlib.pyplot as plt
import Bio.PDB.PDBParser
from Bio.PDB import PDBList
import warnings
warnings.filterwarnings('ignore')

###Calculate pairwise distance matrix####
def pairwisedist(coords):
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
    return dist_matrix

###Get contact map####
def get_contact_map(dist_matrix, threshold=8.0):
    """Binary contact map: 1 if distance < threshold, 0 otherwise."""
    return (dist_matrix < threshold).astype(int)


####Fetch PDB files####

pdb_list = PDBList()
for id in ['1MBN','4HHB','1TIM']:
        pdb_list.retrieve_pdb_file(id, pdir='.', file_format='pdb')
        print(f"Downloded {id}")

def get_ca_atoms(pdb_id):
        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure(pdb_id,f'pdb{id}.ent')
        ca_atoms = []
        residue_ids = []
        for model in structure:
                for chain in model:
                        for residue in chain:
                                if 'CA' in residue:
                                        ca_atoms.append(residue['CA'].get_vector().get_array())
                                        residue_ids.append(residue.get_id()[1])
        return np.array(ca_atoms), residue_ids
             



# ── Process all 3 proteins ───────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

for idx, pdb_id in enumerate(['1MBN','4HHB','1TIM']):
    print(f"\nProcessing {pdb_id}...")
    
    coords, res_ids = get_ca_atoms(pdb_id)
    print(f"  Cα atoms found: {len(coords)}")
    
    dist_matrix = pairwisedist(coords)
    contact_map = get_contact_map(dist_matrix, threshold=8.0)
    
    num_contacts = contact_map.sum() // 2  # divide by 2 (matrix is symmetric)
    print(f"  Total contacts below 8Å: {num_contacts}")
    print(f"  Max distance: {dist_matrix.max():.1f} Å")
    print(f"  Min distance (non-zero): {dist_matrix[dist_matrix>0].min():.1f} Å")
    
    # Plot distance matrix
    im1 = axes[0, idx].imshow(dist_matrix, cmap='viridis')
    axes[0, idx].set_title(f'{pdb_id} — Distance Matrix', fontsize=11)
    axes[0, idx].set_xlabel('Residue index')
    axes[0, idx].set_ylabel('Residue index')
    plt.colorbar(im1, ax=axes[0, idx], label='Distance (Å)')
    
    # Plot contact map
    axes[1, idx].imshow(contact_map, cmap='binary')
    axes[1, idx].set_title(f'{pdb_id} — Contact Map (<8Å)', fontsize=11)
    axes[1, idx].set_xlabel('Residue index')
    axes[1, idx].set_ylabel('Residue index')

plt.suptitle('Cα Distance Matrices and Contact Maps', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('distance_and_contact_maps.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nSaved figure to results/distance_and_contact_maps.png")



