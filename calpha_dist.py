import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Bio.PDB.PDBParser
from Bio.PDB import PDBList

###First exercise is to calculate the pairwise distance between Calpha atoms in the PDB file.
###For this first we need to load the PDB file, then extract the coordinates of c-alpha atoms and finally calculte the pairwise distance between them.
##This exercise explicitly use numpy arrays to store the coordinates and calculate the distances.

##Download the PDB file###
pdbl = PDBList()
pdbl.retrieve_pdb_file('4HHB', pdir='.', file_format='pdb')

###Getting coordinates###
parser = Bio.PDB.PDBParser()
structure = parser.get_structure('4hhb','pdb4hhb.ent')

ca_atoms = []
for model in structure:
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                ca_atoms.append(residue['CA'].get_coord())
ca_atoms = np.array(ca_atoms)
print(ca_atoms.shape)

###Calculate pairwise distance####
def pairwisedist(coords):
    distance = []
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
                diff = coords[i] - coords[j]
                dist = np.sqrt(np.sum(diff**2))
                distance.append(dist)
    return np.array(distance)
pairwise_distance = pairwisedist(ca_atoms)

###The issue with the above code it computes the pairwise distance while store it in a single row, it should give around n(n-1)/2 number of coloumns. rather than n*n matrix
###For that purpose you can use broadcasting of matrix, which basically turns the matrix into a row and then again into a coloumn, on computing the distance, the diaganol elements will be zero.
def pairwise_distance_matrix(coords_matrix):
    diff = coords_matrix[:, np.newaxis, :] - coords_matrix[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=-1))
    return dist
distances_matrix = pairwise_distance_matrix(ca_atoms)

###Plot a heat map which shows how the distance behaves over different residues###

plt.imshow(distances_matrix, cmap='viridis')
plt.colorbar()
plt.title('Pairwise Distance between C-alpha Atoms')
plt.xlabel('C-alpha Atom Index')
plt.ylabel('C-alpha Atom Index')
plt.show()
