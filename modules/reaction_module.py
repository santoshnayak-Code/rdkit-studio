from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
import base64
import io

def create_reaction(reaction_smarts: str):
    """Parse a chemical reaction from Reaction SMARTS."""
    if not reaction_smarts.strip():
        return None, "Empty Reaction SMARTS."
        
    try:
        rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return None, "Invalid Reaction SMARTS."
        return rxn, ""
    except Exception as e:
        return None, str(e)

def run_reaction(rxn, reactants_mols, protect_smarts: str = ""):
    """Run a chemical reaction on reactants, optionally protecting specific atom types."""
    if rxn is None or not reactants_mols:
        return [], "No reaction or reactants provided."
        
    try:
        # Create copies of reactants to avoid side-effects
        reactants = [Chem.Mol(r) for r in reactants_mols if r is not None]
        
        if len(reactants) != rxn.GetNumReactantTemplates():
            return [], f"Expected {rxn.GetNumReactantTemplates()} reactants, but got {len(reactants)}."
            
        # Apply protection if SMARTS query is provided
        if protect_smarts.strip():
            prot_query = Chem.MolFromSmarts(protect_smarts)
            if prot_query is not None:
                for r in reactants:
                    matches = r.GetSubstructMatches(prot_query)
                    for match in matches:
                        # Set RDKit's special protection property on matched atoms
                        for idx in match:
                            r.GetAtomWithIdx(idx).SetProp('_protected', '1')
        
        # Run reactants
        products_sets = rxn.RunReactants(tuple(reactants))
        
        # Format products
        formatted_products = []
        for p_set in products_sets:
            cleaned_set = []
            for p in p_set:
                try:
                    Chem.SanitizeMol(p)
                except Exception:
                    pass
                cleaned_set.append(p)
            formatted_products.append(cleaned_set)
            
        return formatted_products, ""
    except Exception as e:
        return [], str(e)

def draw_reaction_visualization(rxn, width=800, height=250):
    """Draw reaction scheme as HTML image tag containing base64 encoded PNG."""
    if rxn is None:
        return "<div>No reaction loaded.</div>"
        
    try:
        # Fallback to image drawing since it's most robust across platforms/installations
        img = rdChemReactions.ReactionToImage(rxn, subImgSize=(180, 150))
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        return f"""
        <div class="chemical-card" style="display: flex; justify-content: center; align-items: center; padding: 12px; margin-bottom: 20px;">
            <img src="data:image/png;base64,{img_str}" style="max-width: 100%; height: auto; border-radius: 8px;" />
        </div>
        """
    except Exception as e:
        return f"<div style='color: #ef4444;'>Failed to draw reaction: {str(e)}</div>"
