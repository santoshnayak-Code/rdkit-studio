from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import io
import gzip

def mol_from_input(input_str: str, input_format: str):
    """Parse a single molecule from SMILES, MolBlock, or InChI."""
    input_str = input_str.strip()
    if not input_str:
        return None
    
    try:
        if input_format == "SMILES":
            mol = Chem.MolFromSmiles(input_str)
        elif input_format == "MolBlock":
            mol = Chem.MolFromMolBlock(input_str)
        elif input_format == "InChI":
            mol = Chem.MolFromInchi(input_str)
        else:
            mol = None
        return mol
    except Exception:
        return None

def parse_sdf_data(file_content: bytes, is_gz: bool = False):
    """Parse molecules from SDF file content bytes."""
    mols = []
    names = []
    failures = 0
    
    try:
        if is_gz:
            decompressed = gzip.decompress(file_content)
            raw_data = io.BytesIO(decompressed)
        else:
            raw_data = io.BytesIO(file_content)
        
        # Use ForwardSDMolSupplier for file-like object parsing
        suppl = Chem.ForwardSDMolSupplier(raw_data)
        for i, mol in enumerate(suppl):
            if mol is not None:
                mols.append(mol)
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{i+1}"
                if not name.strip():
                    name = f"Mol_{i+1}"
                names.append(name)
            else:
                failures += 1
    except Exception as e:
        # Fallback to standard reading string-based if it fails
        pass
        
    return mols, names, failures

def parse_smiles_data(file_content: str):
    """Parse molecules from SMILES file string."""
    mols = []
    names = []
    failures = 0
    
    lines = file_content.strip().splitlines()
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        parts = line.split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else f"Mol_{i+1}"
        
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mol.SetProp("_Name", name)
                mols.append(mol)
                names.append(name)
            else:
                failures += 1
        except Exception:
            failures += 1
            
    return mols, names, failures

def mol_to_smiles(mol, kekule=False, isomeric=True):
    """Convert Mol object to SMILES."""
    if mol is None:
        return ""
    try:
        m = Chem.Mol(mol) # Copy the mol
        if kekule:
            Chem.Kekulize(m)
        return Chem.MolToSmiles(m, kekuleSmiles=kekule, isomericSmiles=isomeric)
    except Exception as e:
        return f"Error: {str(e)}"

def mol_to_molblock(mol, name=""):
    """Convert Mol object to MDL Mol Block."""
    if mol is None:
        return ""
    try:
        m = Chem.Mol(mol)
        if name:
            m.SetProp("_Name", name)
        # Ensure 2D coordinates exist for a clean depiction in block
        if m.GetNumConformers() == 0:
            AllChem.Compute2DCoords(m)
        return Chem.MolToMolBlock(m)
    except Exception as e:
        return f"Error: {str(e)}"

def mols_to_sdf_bytes(mols):
    """Convert list of mols to SDF bytes."""
    output = io.StringIO()
    writer = Chem.SDWriter(output)
    for m in mols:
        if m is not None:
            writer.write(m)
    writer.close()
    return output.getvalue().encode('utf-8')

def mols_to_smiles_bytes(mols):
    """Convert list of mols to SMILES file bytes."""
    output = io.StringIO()
    for i, m in enumerate(mols):
        if m is not None:
            smi = Chem.MolToSmiles(m)
            name = m.GetProp("_Name") if m.HasProp("_Name") else f"Mol_{i+1}"
            output.write(f"{smi}\t{name}\n")
    return output.getvalue().encode('utf-8')

def get_mols_metadata_df(mols, names):
    """Generate a pandas DataFrame of batch metadata."""
    data = []
    for m, name in zip(mols, names):
        if m is not None:
            data.append({
                "Name": name,
                "SMILES": Chem.MolToSmiles(m),
                "Atoms": m.GetNumAtoms(),
                "Bonds": m.GetNumBonds(),
                "Heavy Atoms": m.GetNumHeavyAtoms(),
                "Formula": Chem.rdMolDescriptors.CalcMolFormula(m)
            })
    return pd.DataFrame(data)
