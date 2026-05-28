import sys
import rdkit
from rdkit import Chem
print(f"RDKit Version: {rdkit.__version__}")
print(f"Python Version: {sys.version}")

try:
    from modules import io_module, structure_module, substructure_module, transformation_module, fingerprint_module, descriptor_module, reaction_module
    print("All custom modules imported successfully!")
except Exception as e:
    print(f"Import Error: {str(e)}")
    sys.exit(1)

# Sanity check 1: Parsing and drawing
mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
if mol is None:
    print("Failed to parse default SMILES.")
    sys.exit(1)
print("Sanity Check 1 passed: Default SMILES parsed.")

# Sanity check 2: 2D drawing
svg = structure_module.draw_molecule_svg(mol)
if "<svg" not in svg:
    print("Failed to generate molecular drawing SVG.")
    sys.exit(1)
print("Sanity Check 2 passed: SVG drawing generated.")

# Sanity check 3: Descriptors
desc = descriptor_module.calculate_key_descriptors(mol)
if "Molecular Weight" not in desc or desc["Molecular Weight"] < 200:
    print("Failed to calculate descriptors.")
    sys.exit(1)
print("Sanity Check 3 passed: Descriptors calculated successfully.")

# Sanity check 4: Similarity
ref = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
probe = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
sim = fingerprint_module.calculate_similarity(ref, probe, "Morgan (Circular)", "Tanimoto")
print(f"Sanity Check 4 passed: Tanimoto Similarity is {sim:.4f}")

# Sanity check 5: Reactions
rxn, err = reaction_module.create_reaction("[C:1](=[O:2])[O:3].[N:4]>>[C:1](=[O:2])[N:4]")
if err or rxn is None:
    print(f"Failed to parse reaction smarts: {err}")
    sys.exit(1)
print("Sanity Check 5 passed: Reaction SMARTS parsed.")

# Sanity check 6: Atoms and Bonds tables
df_atoms = structure_module.get_atoms_dataframe(mol)
df_bonds = structure_module.get_bonds_dataframe(mol)
if df_atoms.empty or df_bonds.empty:
    print("Failed to generate atoms/bonds dataframes.")
    sys.exit(1)
print("Sanity Check 6 passed: Atoms and Bonds detailed tables generated.")

print("All RDKit Studio Sanity Checks Passed Successfully!")
sys.exit(0)
