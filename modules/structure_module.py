from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import pandas as pd
import streamlit.components.v1 as components

def get_atoms_dataframe(mol):
    """Retrieve details of all atoms in a molecule."""
    if mol is None:
        return pd.DataFrame()
    
    atoms_data = []
    for atom in mol.GetAtoms():
        atoms_data.append({
            "Idx": atom.GetIdx(),
            "Symbol": atom.GetSymbol(),
            "Atomic Num": atom.GetAtomicNum(),
            "Formal Charge": atom.GetFormalCharge(),
            "Hybridization": str(atom.GetHybridization()),
            "Valence (Val/Exp/Imp)": f"{atom.GetTotalValence()}/{atom.GetExplicitValence()}/{atom.GetImplicitValence()}",
            "Aromatic?": "Yes" if atom.GetIsAromatic() else "No",
            "In Ring?": "Yes" if atom.IsInRing() else "No",
            "Ring Sizes (3-6)": ", ".join([str(s) for s in [3, 4, 5, 6] if atom.IsInRingSize(s)]) or "None"
        })
    return pd.DataFrame(atoms_data)

def get_bonds_dataframe(mol):
    """Retrieve details of all bonds in a molecule."""
    if mol is None:
        return pd.DataFrame()
    
    bonds_data = []
    for bond in mol.GetBonds():
        bonds_data.append({
            "Idx": bond.GetIdx(),
            "Begin Atom": f"{mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()} (Idx {bond.GetBeginAtomIdx()})",
            "End Atom": f"{mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()} (Idx {bond.GetEndAtomIdx()})",
            "Bond Type": str(bond.GetBondType()),
            "Aromatic?": "Yes" if bond.GetIsAromatic() else "No",
            "In Ring?": "Yes" if bond.IsInRing() else "No",
            "Ring Sizes (3-6)": ", ".join([str(s) for s in [3, 4, 5, 6] if bond.IsInRingSize(s)]) or "None"
        })
    return pd.DataFrame(bonds_data)

def toggle_hydrogens(mol, add_hs=True):
    """Add or remove Hydrogen atoms from the molecule."""
    if mol is None:
        return None
    try:
        if add_hs:
            return Chem.AddHs(mol)
        else:
            return Chem.RemoveHs(mol)
    except Exception:
        return mol

def draw_molecule_svg(mol, width=400, height=350, highlight_atoms=None, highlight_atom_colors=None, 
                      highlight_bonds=None, highlight_bond_colors=None, atom_notes=None, bond_notes=None,
                      add_stereo=False, add_indices=False):
    """Generate SVG depiction of a molecule with customized highlights and notes."""
    if mol is None:
        return ""
    
    try:
        # Work on a copy to avoid polluting the original molecule with notes
        m = Chem.Mol(mol)
        
        # In RDKit, MolToMolBlock or depictions require coordinate generation first
        if m.GetNumConformers() == 0:
            AllChem.Compute2DCoords(m)
        
        # Apply atom notes
        if atom_notes:
            for idx, note in atom_notes.items():
                if 0 <= idx < m.GetNumAtoms():
                    m.GetAtomWithIdx(idx).SetProp('atomNote', str(note))
                    
        # Apply bond notes
        if bond_notes:
            for idx, note in bond_notes.items():
                if 0 <= idx < m.GetNumBonds():
                    m.GetBondWithIdx(idx).SetProp('bondNote', str(note))
        
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        options = drawer.drawOptions()
        options.addStereoAnnotation = add_stereo
        options.addAtomIndices = add_indices
        options.clearBackground = False # Let transparent show style card background
        
        # Convert colors to RGB tuples (range 0 to 1) for RDKit
        rd_atom_colors = {}
        if highlight_atom_colors:
            for idx, color_rgb in highlight_atom_colors.items():
                rd_atom_colors[idx] = color_rgb
                
        rd_bond_colors = {}
        if highlight_bond_colors:
            for idx, color_rgb in highlight_bond_colors.items():
                rd_bond_colors[idx] = color_rgb
        
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, m, 
            highlightAtoms=highlight_atoms or [],
            highlightAtomColors=rd_atom_colors or None,
            highlightBonds=highlight_bonds or [],
            highlightBondColors=rd_bond_colors or None
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        return f"<svg><text y='20'>Error drawing molecule: {str(e)}</text></svg>"

def generate_3d_conformer(mol, random_seed=42):
    """Generate a optimized 3D conformer using the ETKDGv3 method and MMFF94."""
    if mol is None:
        return None
    try:
        m = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = random_seed
        
        embed_status = AllChem.EmbedMolecule(m, params)
        if embed_status != 0:
            # Fallback if ETKDG fails
            embed_status = AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=random_seed)
            
        if embed_status == 0:
            # Optimize structure
            try:
                AllChem.MMFFOptimizeMolecule(m)
            except Exception:
                pass
            return m
        return None
    except Exception:
        return None

def render_3dmol_viewer(mol_3d, height=450, style="stick"):
    """Render 3D molecular viewer using HTML5 3Dmol.js."""
    if mol_3d is None:
        return
    
    try:
        pdb_block = Chem.MolToPDBBlock(mol_3d)
        # Escape backticks and double quotes in PDB data to avoid JS breakages
        pdb_data_escaped = pdb_block.replace("`", "\\`").replace("\n", "\\n")
        
        style_js = "{stick: {colorscheme: 'Jmol', radius: 0.15}}"
        if style == "sphere":
            style_js = "{sphere: {colorscheme: 'Jmol', scale: 0.3}}"
        elif style == "line":
            style_js = "{line: {colorscheme: 'Jmol'}}"
        
        html_code = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
            <style>
                body {{
                    margin: 0;
                    padding: 0;
                    background-color: #0b0f17;
                    overflow: hidden;
                }}
                #viewer {{
                    width: 100vw;
                    height: 100vh;
                    position: relative;
                }}
            </style>
        </head>
        <body>
            <div id="viewer"></div>
            <script>
                var element = document.getElementById('viewer');
                var config = {{ backgroundColor: '#0f172a' }};
                var viewer = $3Dmol.createViewer(element, config);
                var modelData = `{pdb_data_escaped}`;
                viewer.addModel(modelData, "pdb");
                viewer.setStyle({{}}, {style_js});
                viewer.zoomTo();
                viewer.render();
            </script>
        </body>
        </html>
        """
        components.html(html_code, height=height)
    except Exception as e:
        st.error(f"Failed to render 3D Viewer: {str(e)}")

def mol_from_png_metadata(png_bytes):
    """Reconstruct a molecule from metadata embedded in a PNG image."""
    try:
        mol = Chem.MolFromPNGString(png_bytes)
        return mol
    except Exception:
        return None
