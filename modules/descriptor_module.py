from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def calculate_key_descriptors(mol):
    """Calculate key physical and chemical descriptors for a molecule."""
    if mol is None:
        return {}
        
    try:
        # Calculate descriptors using built-in methods
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        arom = Descriptors.NumAromaticRings(mol)
        f_csp3 = Descriptors.FractionCSP3(mol)
        
        # QED score
        qed_val = QED.qed(mol)
        
        # Formula
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        
        return {
            "Formula": formula,
            "Molecular Weight": round(mw, 2),
            "LogP (Octanol/Water)": round(logp, 2),
            "TPSA (Å²)": round(tpsa, 2),
            "H-Bond Donors": hbd,
            "H-Bond Acceptors": hba,
            "Rotatable Bonds": rotb,
            "Aromatic Rings": arom,
            "Fraction Csp3": round(f_csp3, 3),
            "QED Drug-likeness": round(qed_val, 3)
        }
    except Exception as e:
        return {"Error": str(e)}

def get_descriptors_dataframe(mol):
    """Retrieve descriptors as a pandas DataFrame."""
    desc = calculate_key_descriptors(mol)
    if not desc or "Error" in desc:
        return pd.DataFrame()
    return pd.DataFrame(list(desc.items()), columns=["Descriptor", "Value"])

def get_radar_chart_plotly(mol):
    """Generate an interactive Plotly radar chart of drug-likeness parameters."""
    desc = calculate_key_descriptors(mol)
    if not desc or "Error" in desc:
        return None
        
    # Standardize parameters from 0 to 1 for radar plot
    # Values normalized around typical Lipinski/Drug-like thresholds
    mw_norm = min(desc["Molecular Weight"] / 500.0, 1.2)
    logp_norm = min((desc["LogP (Octanol/Water)"] + 2.0) / 7.0, 1.2) # typical range -2 to 5
    tpsa_norm = min(desc["TPSA (Å²)"] / 140.0, 1.2)
    hbd_norm = min(desc["H-Bond Donors"] / 5.0, 1.2)
    hba_norm = min(desc["H-Bond Acceptors"] / 10.0, 1.2)
    rotb_norm = min(desc["Rotatable Bonds"] / 10.0, 1.2)
    qed_val = desc["QED Drug-likeness"]
    
    categories = [
        'Mol Wt (norm)', 'LogP (norm)', 'TPSA (norm)', 
        'HBD (norm)', 'HBA (norm)', 'Rotatable Bonds (norm)', 'QED Score'
    ]
    values = [mw_norm, logp_norm, tpsa_norm, hbd_norm, hba_norm, rotb_norm, qed_val]
    
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=values,
        theta=categories,
        fill='toself',
        fillcolor='rgba(6, 182, 212, 0.2)',
        line=dict(color='#06b6d4', width=2),
        name='Properties'
    ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1.2],
                color='#94a3b8',
                gridcolor='rgba(255,255,255,0.05)',
                tickvals=[0.2, 0.5, 0.8, 1.0]
            ),
            angularaxis=dict(
                color='#e2e8f0',
                gridcolor='rgba(255,255,255,0.05)'
            ),
            bgcolor='rgba(15, 23, 42, 0.1)'
        ),
        showlegend=False,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        margin=dict(l=40, r=40, t=20, b=20),
        height=320
    )
    return fig

def calculate_gasteiger_charges(mol):
    """Compute Gasteiger partial charges and return list of charges."""
    if mol is None:
        return []
        
    try:
        # Work on a copy
        m = Chem.Mol(mol)
        AllChem.ComputeGasteigerCharges(m)
        
        charges = []
        for atom in m.GetAtoms():
            if atom.HasProp('_GasteigerCharge'):
                val = float(atom.GetProp('_GasteigerCharge'))
                # Handle cases where Gasteiger charge is nan or inf
                if val != val or val == float('inf') or val == float('-inf'):
                    val = 0.0
                charges.append(val)
            else:
                charges.append(0.0)
        return charges
    except Exception:
        return [0.0] * mol.GetNumAtoms()

def draw_charge_distribution_svg(mol, width=400, height=350):
    """Draw a molecule color-coding atoms by partial Gasteiger charge."""
    if mol is None:
        return ""
        
    try:
        m = Chem.Mol(mol)
        charges = calculate_gasteiger_charges(m)
        
        # Color mapper: map Gasteiger charge (typically -0.5 to 0.5) to RGB
        # Negative charges: Red (high charge intensity -> pure red)
        # Positive charges: Blue (high charge intensity -> pure blue)
        # Near zero charges: Light gray/white
        
        highlight_atom_colors = {}
        highlight_atoms = []
        atom_notes = {}
        
        for idx, q in enumerate(charges):
            highlight_atoms.append(idx)
            atom_notes[idx] = f"{q:.2f}"
            
            # Clamp charge to range [-0.5, 0.5] for normalization
            q_norm = max(min(q, 0.5), -0.5)
            
            if q_norm < 0:
                # Map negative to Red-ish color
                factor = abs(q_norm) * 2.0 # range 0 to 1
                r = 1.0
                g = 1.0 - factor * 0.8
                b = 1.0 - factor * 0.8
            else:
                # Map positive to Blue-ish color
                factor = q_norm * 2.0 # range 0 to 1
                r = 1.0 - factor * 0.8
                g = 1.0 - factor * 0.8
                b = 1.0
                
            highlight_atom_colors[idx] = (r, g, b)
            
        # Draw SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        options = drawer.drawOptions()
        options.clearBackground = False
        
        # Set atom charge text annotations
        for idx, text in atom_notes.items():
            m.GetAtomWithIdx(idx).SetProp('atomNote', text)
            
        if m.GetNumConformers() == 0:
            AllChem.Compute2DCoords(m)
            
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, m,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_atom_colors
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        return f"<svg><text y='20'>Error generating charge view: {str(e)}</text></svg>"
