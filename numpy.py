import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio.PDB import PDBList

PDB = PDB.List().retrieve_pdb_file('4HHB')