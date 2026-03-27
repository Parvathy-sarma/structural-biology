import numpy as np
from Bio.PDB import PDBList, PDBParser
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split  # ML: used to split data into training and testing sets

# -----------------------------
# Step 1 — Download & load PDB
# -----------------------------
test_PDB = PDBList()
test_PDB.retrieve_pdb_file("1UBQ", pdir = ".", file_format= 'pdb')
parse = PDBParser()
structure = parse.get_structure(test_PDB, 'pdb1ubq.ent')
com = structure.center_of_mass()
com_protein = np.array(com)

# -----------------------------
# Step 2 — Extract C-alpha coordinates
# -----------------------------
ca_atoms = []
res_id = []

for model in structure:
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                ca_atoms.append(residue['CA'].get_coord())
                res_id.append(residue.get_id()[1])
ca_atoms = np.array(ca_atoms)
residue_id = np.array(res_id)
print(com_protein)

# -----------------------------
# Step 3 — Distance from center of mass
# -----------------------------
def dist_com(coords):
    diff = coords - com_protein
    dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
    return dist_matrix

distance = dist_com(ca_atoms)
print(distance.shape)

# -----------------------------
# Step 4 — Create binary labels
# -----------------------------
threshold = np.median(distance)
labels = (distance > threshold).astype(int)
# ML Explanation:
# Each residue is labeled as '0' (buried) or '1' (exposed) based on its distance from the center.
# This is what the ML model will try to predict.

# -----------------------------
# Step 5 — Pairwise distance & contact map (not ML, for structure visualization)
# -----------------------------
def pairwise_distance(coords):
    diff = coords[:,np.newaxis,:] - coords[np.newaxis,:,:]
    dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
    return dist_matrix
dist_com = pairwise_distance(ca_atoms)

def contactmaps(dist_matrix,threshold=8):
    return(dist_matrix < threshold).astype(int)
contact = contactmaps(dist_com)

# -----------------------------
# Step 6 — Prepare data for ML
# -----------------------------
X = distance.reshape(-1,1)  # ML: Features
y = labels                  # ML: Target labels
# Explanation:
# X = what the model sees (distance values)
# y = what the model should predict (buried vs exposed)

# -----------------------------
# Step 7 — Split dataset into training and testing
# -----------------------------
X_train, X_test, y_train, y_test = train_test_split(X,y, test_size =0.3, random_state=42)
# Explanation:
# Training set = used to teach the model
# Test set = used to see if the model can correctly classify new residues
# test_size=0.3 → 30% of residues are reserved for testing

# -----------------------------
# Step 8 — Train a simple ML model
# -----------------------------
from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
model.fit(X_train,y_train)
# Explanation:
# Logistic Regression is a simple algorithm that learns the rule mapping distance → buried/exposed.
# The model uses the training data to figure out the threshold separating classes.

# -----------------------------
# Step 9 — Make predictions on test set
# -----------------------------
from sklearn.metrics import accuracy_score, classification_report

y_pred = model.predict(X_test)
# Explanation:
# The model predicts labels for residues it hasn't seen (test set)
# This tests if it learned the pattern correctly

# -----------------------------
# Step 10 — Evaluate model performance
# -----------------------------
print("Accuracy:", accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))
# Explanation:
# Accuracy: fraction of correct predictions
# Classification report: gives precision, recall, f1-score
# Precision → when model predicts exposed, how often is it correct
# Recall → how many of the true exposed residues were found
# F1 → balance between precision and recall

# -----------------------------
# Step 11 — Visualize ML results
# -----------------------------
plt.figure()
plt.scatter(X, y, c=y)

x_vals = np.linspace(X.min(), X.max(), 100).reshape(-1, 1)
y_vals = model.predict(x_vals)

plt.plot(x_vals, y_vals)
plt.xlabel("Distance from COM")
plt.ylabel("Class (0=buried, 1=exposed)")
plt.title("ML Classification of Residue Exposure")
plt.show()
# Explanation:
# Scatter plot shows actual residues (buried/exposed)
# The line shows the decision boundary the model learned
# This is an intuitive way to see what the model is “thinking”

# -----------------------------
# Step 12 — Print dataset info
# -----------------------------
print("Number of residues:", len(distance))
print("Threshold used:", threshold)

