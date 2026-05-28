from rdkit import Chem

def find_substructure_matches(mol, query_str: str, query_type: str = "SMARTS", use_chirality: bool = False):
    """Search for substructures using a SMARTS or SMILES query."""
    if mol is None or not query_str.strip():
        return False, [], None, ""
    
    try:
        if query_type == "SMARTS":
            query = Chem.MolFromSmarts(query_str)
            label = "SMARTS"
        else:
            query = Chem.MolFromSmiles(query_str)
            label = "SMILES"
            
        if query is None:
            return False, [], None, f"Invalid {label} syntax."
            
        matches = mol.GetSubstructMatches(query, useChirality=use_chirality)
        has_match = len(matches) > 0
        return has_match, list(matches), query, ""
    except Exception as e:
        return False, [], None, str(e)

def get_highlight_atoms_and_bonds(mol, query, matches):
    """Retrieve atom and bond indices to highlight for all matches."""
    if not matches or query is None or mol is None:
        return [], []
        
    highlight_atoms = set()
    highlight_bonds = set()
    
    for match in matches:
        # Add all matching atom indices
        for idx in match:
            highlight_atoms.add(idx)
            
        # Find matching bonds between the matching atoms
        for bond in query.GetBonds():
            aid1 = match[bond.GetBeginAtomIdx()]
            aid2 = match[bond.GetEndAtomIdx()]
            b = mol.GetBondBetweenAtoms(aid1, aid2)
            if b is not None:
                highlight_bonds.add(b.GetIdx())
                
    return list(highlight_atoms), list(highlight_bonds)

def get_atom_map_indices(query):
    """Retrieve atom map indices property from a query molecule."""
    if query is None:
        return {}
        
    ind_map = {}
    for atom in query.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num:
            # map number to internal atom index
            ind_map[map_num] = atom.GetIdx()
    return ind_map

class SidechainChecker(object):
    """Checks whether the sidechains connected to matched atoms satisfy specific constraints."""
    def __init__(self, query, prop_name="queryType"):
        self.ats_to_examine = [
            (x.GetIdx(), x.GetProp(prop_name)) 
            for x in query.GetAtoms() 
            if x.HasProp(prop_name)
        ]
        self.prop_name = prop_name
        self.matchers = {
            'alkyl': lambda at: not at.GetIsAromatic(),
            'all_carbon': lambda at: at.GetAtomicNum() == 6,
            'no_nitrogen': lambda at: at.GetAtomicNum() != 7
        }

    def __call__(self, mol, match_indices):
        seen = [0] * mol.GetNumAtoms()
        for idx in match_indices:
            seen[idx] = 1
            
        for q_idx, q_type in self.ats_to_examine:
            if q_type not in self.matchers:
                continue
            
            m_idx = match_indices[q_idx]
            stack = [mol.GetAtomWithIdx(m_idx)]
            
            while stack:
                atom = stack.pop(0)
                if not self.matchers[q_type](atom):
                    return False
                seen[atom.GetIdx()] = 1
                for nbr in atom.GetNeighbors():
                    if not seen[nbr.GetIdx()]:
                        stack.append(nbr)
        return True

def filter_matches_with_checker(mol, query, matches, constraint_type):
    """Filter list of matches using custom sidechain constraints."""
    if not matches or query is None:
        return []
        
    # Annotate query first atom to test constraint
    q_copy = Chem.Mol(query)
    # Give the first atom of the query the constraint type if there is a match
    # To demonstrate, we set queryType property on atom 0
    q_copy.GetAtomWithIdx(0).SetProp("queryType", constraint_type)
    
    checker = SidechainChecker(q_copy)
    filtered = []
    for match in matches:
        if checker(mol, match):
            filtered.append(match)
    return filtered
