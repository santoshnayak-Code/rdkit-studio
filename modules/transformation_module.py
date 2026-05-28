from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS
from rdkit.Chem import BRICS
from rdkit.Chem.Recap import RecapDecompose

def delete_substructure(mol, query_smarts: str):
    """Delete all occurrences of a substructure from a molecule."""
    if mol is None or not query_smarts.strip():
        return None, "Empty inputs."
        
    try:
        query = Chem.MolFromSmarts(query_smarts)
        if query is None:
            return None, "Invalid query SMARTS."
        
        prod = Chem.DeleteSubstructs(mol, query)
        # Sanitize to ensure valence and aromaticity are corrected
        Chem.SanitizeMol(prod)
        return prod, ""
    except Exception as e:
        return None, str(e)

def replace_substructure(mol, query_smarts: str, replacement_smiles: str):
    """Replace all occurrences of a substructure with a replacement group."""
    if mol is None or not query_smarts.strip() or not replacement_smiles.strip():
        return None, "Empty inputs."
        
    try:
        query = Chem.MolFromSmarts(query_smarts)
        if query is None:
            return None, "Invalid query SMARTS."
            
        repl = Chem.MolFromSmiles(replacement_smiles)
        if repl is None:
            return None, "Invalid replacement SMILES."
            
        # ReplaceSubstructs returns a tuple of molecules (all possible replacements)
        replaced_mols = Chem.ReplaceSubstructs(mol, query, repl)
        if replaced_mols:
            # Let's take the first result
            prod = replaced_mols[0]
            Chem.SanitizeMol(prod)
            return prod, ""
        return None, "No replacements could be made."
    except Exception as e:
        return None, str(e)

def replace_core(mol, core_smarts: str):
    """Isolate sidechains by cutting out a core substructure."""
    if mol is None or not core_smarts.strip():
        return None, "Empty inputs."
        
    try:
        core = Chem.MolFromSmarts(core_smarts)
        if core is None:
            return None, "Invalid core SMARTS."
            
        sidechains = Chem.ReplaceCore(mol, core)
        if sidechains is not None:
            # Return as SMILES or Mol
            Chem.SanitizeMol(sidechains)
            return sidechains, ""
        return None, "Core not found in molecule."
    except Exception as e:
        return None, str(e)

def get_murcko_scaffold(mol):
    """Calculate the Bemis-Murcko scaffold of a molecule."""
    if mol is None:
        return None, "Invalid molecule."
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return scaffold, ""
    except Exception as e:
        return None, str(e)

def decompose_brics(mol):
    """Decompose a molecule into fragments using BRICS rules."""
    if mol is None:
        return [], "Invalid molecule."
    try:
        # BRICSDecompose returns a set of SMILES strings of fragments
        frags = BRICS.BRICSDecompose(mol)
        # Convert SMILES to Mol objects
        frag_mols = []
        for f in frags:
            m = Chem.MolFromSmiles(f)
            if m is not None:
                frag_mols.append(m)
        return frag_mols, ""
    except Exception as e:
        return [], str(e)

def decompose_recap(mol):
    """Decompose a molecule into fragments using RECAP rules."""
    if mol is None:
        return [], "Invalid molecule."
    try:
        res = RecapDecompose(mol)
        # Get children nodes of decomposition tree
        leaves = res.GetLeaves()
        frag_mols = []
        for k in leaves.keys():
            m = Chem.MolFromSmiles(k)
            if m is not None:
                frag_mols.append(m)
        return frag_mols, ""
    except Exception as e:
        return [], str(e)

def find_maximum_common_substructure(mols, atom_compare="CompareElements", bond_compare="CompareOrder", 
                                      match_valences=True, ring_matches_ring=False):
    """Calculate the Maximum Common Substructure (MCS) of a list of molecules."""
    if not mols or len(mols) < 2:
        return None, "At least 2 molecules are required."
        
    try:
        # Map parameters
        atom_comp = rdFMCS.AtomCompare.CompareElements
        if atom_compare == "CompareAny":
            atom_comp = rdFMCS.AtomCompare.CompareAny
        elif atom_compare == "CompareAnyHeavyAtom":
            atom_comp = rdFMCS.AtomCompare.CompareAnyHeavyAtom
            
        bond_comp = rdFMCS.BondCompare.CompareOrder
        if bond_compare == "CompareAny":
            bond_comp = rdFMCS.BondCompare.CompareAny
        elif bond_compare == "CompareOrderExact":
            bond_comp = rdFMCS.BondCompare.CompareOrderExact
            
        mcs_res = rdFMCS.FindMCS(
            mols,
            atomCompare=atom_comp,
            bondCompare=bond_comp,
            matchValences=match_valences,
            ringMatchesRingOnly=ring_matches_ring
        )
        
        if mcs_res.numAtoms == 0:
            return None, "No common substructure found."
            
        # Parse result
        mcs_mol = Chem.MolFromSmarts(mcs_res.smartsString)
        return {
            "smarts": mcs_res.smartsString,
            "mol": mcs_mol,
            "num_atoms": mcs_res.numAtoms,
            "num_bonds": mcs_res.numBonds,
            "timeout": mcs_res.canceled
        }, ""
    except Exception as e:
        return None, str(e)
