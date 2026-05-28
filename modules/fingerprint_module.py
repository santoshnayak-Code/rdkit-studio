from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
import io
import matplotlib
matplotlib.use('Agg') # Non-interactive backend for server safety
import matplotlib.pyplot as plt

def calculate_fingerprint(mol, fp_type: str, radius: int = 2, n_bits: int = 2048):
    """Generate fingerprint bit vector or count vector based on selected type."""
    if mol is None:
        return None
        
    try:
        if fp_type == "RDKit (Topological)":
            return Chem.RDKFingerprint(mol, fpSize=n_bits)
        elif fp_type == "Morgan (Circular)":
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
        elif fp_type == "Atom Pairs":
            return AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
        elif fp_type == "Topological Torsions":
            return AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=n_bits)
        elif fp_type == "MACCS Keys":
            return MACCSkeys.GenMACCSKeys(mol)
        return None
    except Exception:
        return None

def calculate_similarity(mol1, mol2, fp_type: str, metric: str, radius: int = 2, n_bits: int = 2048, 
                         tversky_alpha: float = 0.5, tversky_beta: float = 0.5):
    """Compute similarity between two molecules using chosen fingerprint and metric."""
    if mol1 is None or mol2 is None:
        return 0.0
        
    try:
        fp1 = calculate_fingerprint(mol1, fp_type, radius, n_bits)
        fp2 = calculate_fingerprint(mol2, fp_type, radius, n_bits)
        
        if fp1 is None or fp2 is None:
            return 0.0
            
        if metric == "Tanimoto":
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        elif metric == "Dice":
            return DataStructs.DiceSimilarity(fp1, fp2)
        elif metric == "Cosine":
            return DataStructs.CosineSimilarity(fp1, fp2)
        elif metric == "Sokal":
            return DataStructs.SokalSimilarity(fp1, fp2)
        elif metric == "Kulczynski":
            return DataStructs.KulczynskiSimilarity(fp1, fp2)
        elif metric == "Tversky":
            return DataStructs.TverskySimilarity(fp1, fp2, tversky_alpha, tversky_beta)
        return 0.0
    except Exception:
        return 0.0

def get_morgan_bit_info(mol, radius=2, n_bits=2048):
    """Retrieve full list of active bits and their environments for Morgan fingerprint."""
    if mol is None:
        return {}, None
    info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits, bitInfo=info)
    return info, fp

def draw_morgan_bit_env(mol, bit_id, bit_info, width=300, height=300):
    """Render the exact sub-environment for a specific Morgan bit ID."""
    if mol is None or not bit_info or bit_id not in bit_info:
        return None
    try:
        # DrawMorganBit returns a PIL Image
        img = Draw.DrawMorganBit(mol, bit_id, bit_info, useSVG=True)
        return img
    except Exception:
        return None

def get_rdkit_bit_info(mol, n_bits=2048):
    """Retrieve active bits and environments for topological RDKit fingerprint."""
    if mol is None:
        return {}, None
    info = {}
    fp = Chem.RDKFingerprint(mol, fpSize=n_bits, bitInfo=info)
    return info, fp

def draw_rdkit_bit_env(mol, bit_id, bit_info):
    """Render the exact sub-environment for a specific RDKit bit ID."""
    if mol is None or not bit_info or bit_id not in bit_info:
        return None
    try:
        img = Draw.DrawRDKitBit(mol, bit_id, bit_info, useSVG=True)
        return img
    except Exception:
        return None

def generate_similarity_map_svg(ref_mol, probe_mol, fp_type="Morgan"):
    """Draw a similarity map of probe relative to reference, returning SVG string."""
    if ref_mol is None or probe_mol is None:
        return ""
        
    try:
        # Clear any active plt drawings
        plt.clf()
        fig, ax = plt.subplots(figsize=(6, 6))
        
        # Decide fingerprint function
        fp_fn = SimilarityMaps.GetMorganFingerprint
        if fp_type == "RDKit":
            fp_fn = SimilarityMaps.GetRDKFingerprint
        elif fp_type == "Atom Pairs":
            fp_fn = SimilarityMaps.GetAPFingerprint
        elif fp_type == "Topological Torsions":
            fp_fn = SimilarityMaps.GetTTFingerprint
            
        # Draw similarity map
        _, max_val = SimilarityMaps.GetSimilarityMapForFingerprint(
            ref_mol, probe_mol, fp_fn, ax=ax
        )
        
        # Save to SVG buffer
        buf = io.StringIO()
        fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True)
        plt.close(fig)
        return buf.getvalue()
    except Exception as e:
        return f"<svg><text y='20'>Error generating similarity map: {str(e)}</text></svg>"

def pick_diverse_mols(mols, num_to_pick, radius=2, n_bits=2048):
    """Select a diverse subset of molecules using Morgan fingerprints and MaxMinPicker."""
    if not mols or len(mols) <= num_to_pick or num_to_pick <= 0:
        return list(range(len(mols)))
        
    try:
        # Generate fingerprints for all molecules
        fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=radius, nBits=n_bits) for m in mols]
        
        # Define distance metric function (1 - Tanimoto)
        def dist_func(i, j):
            return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
            
        picker = MaxMinPicker()
        # lazy pick returns picked indices list
        picked = picker.LazyPick(dist_func, len(mols), num_to_pick)
        return list(picked)
    except Exception:
        return list(range(min(num_to_pick, len(mols))))
