import streamlit as st
import pandas as pd
from rdkit import Chem
from modules.ui_helpers import init_page_styling, render_svg, card_header, metric_card
from modules import io_module, structure_module, substructure_module, transformation_module, fingerprint_module, descriptor_module, reaction_module

# Configure Streamlit page
st.set_page_config(
    page_title="RDKit Studio",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Load design system & CSS
init_page_styling()

# Initialize Session States with Defaults
if "active_mol" not in st.session_state:
    st.session_state.active_mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O") # Ibuprofen
    st.session_state.active_name = "Ibuprofen"

if "batch_mols" not in st.session_state:
    # Default batch dataset
    smiles_list = [
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen"),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine"),
        ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
        ("CC(=O)Nc1ccc(O)cc1", "Acetaminophen")
    ]
    st.session_state.batch_mols = [Chem.MolFromSmiles(s) for s, _ in smiles_list]
    st.session_state.batch_names = [n for _, n in smiles_list]
    for m, name in zip(st.session_state.batch_mols, st.session_state.batch_names):
        if m:
            m.SetProp("_Name", name)

# Sidebar - Molecular Input
st.sidebar.markdown("""
    <div style="text-align: center; margin-bottom: 20px;">
        <h2 style="margin: 0; background: linear-gradient(135deg, #38bdf8 0%, #06b6d4 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; font-weight: 700;">🧪 RDKit Studio</h2>
        <div style="color: #64748b; font-size: 0.85rem; letter-spacing: 0.05em; text-transform: uppercase;">Chemical Intelligence Suite</div>
    </div>
""", unsafe_allow_html=True)

st.sidebar.markdown("---")
st.sidebar.subheader("Active Molecule Input")

input_method = st.sidebar.radio("Select Input Method", ["SMILES", "MolBlock Text", "Upload File (.mol/.sdf/.smi)"])

new_mol = None
new_name = st.session_state.active_name

if input_method == "SMILES":
    smi_input = st.sidebar.text_input("Enter SMILES String", "CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    new_name = st.sidebar.text_input("Molecule Name", "Ibuprofen")
    if st.sidebar.button("Load SMILES"):
        mol = io_module.mol_from_input(smi_input, "SMILES")
        if mol:
            new_mol = mol
            st.sidebar.success("SMILES parsed successfully!")
        else:
            st.sidebar.error("Invalid SMILES string.")

elif input_method == "MolBlock Text":
    block_input = st.sidebar.text_area("Paste MDL Mol Block Content", height=200)
    new_name = st.sidebar.text_input("Molecule Name", "Custom Molecule")
    if st.sidebar.button("Load Mol Block"):
        mol = io_module.mol_from_input(block_input, "MolBlock")
        if mol:
            new_mol = mol
            st.sidebar.success("MolBlock parsed successfully!")
        else:
            st.sidebar.error("Invalid MolBlock format.")

elif input_method == "Upload File (.mol/.sdf/.smi)":
    uploaded_file = st.sidebar.file_uploader("Upload Chemical File", type=["mol", "sdf", "smi", "txt", "gz"])
    if uploaded_file is not None:
        file_bytes = uploaded_file.read()
        file_name = uploaded_file.name
        
        # Single Mol Upload
        if file_name.endswith(".mol"):
            mol_text = file_bytes.decode("utf-8")
            mol = io_module.mol_from_input(mol_text, "MolBlock")
            if mol:
                new_mol = mol
                new_name = file_name.replace(".mol", "")
                st.sidebar.success(f"Loaded single molecule: {new_name}")
            else:
                st.sidebar.error("Failed to parse Mol file.")
                
        # Batch Upload (SDF)
        elif file_name.endswith(".sdf") or file_name.endswith(".sdf.gz"):
            is_gz = file_name.endswith(".gz")
            mols, names, fails = io_module.parse_sdf_data(file_bytes, is_gz)
            if mols:
                st.session_state.batch_mols = mols
                st.session_state.batch_names = names
                new_mol = mols[0]
                new_name = names[0]
                st.sidebar.success(f"Loaded batch: {len(mols)} molecules ({fails} fails).")
            else:
                st.sidebar.error("Failed to parse SDF batch.")
                
        # Batch Upload (SMILES)
        elif file_name.endswith(".smi") or file_name.endswith(".txt"):
            smi_text = file_bytes.decode("utf-8")
            mols, names, fails = io_module.parse_smiles_data(smi_text)
            if mols:
                st.session_state.batch_mols = mols
                st.session_state.batch_names = names
                new_mol = mols[0]
                new_name = names[0]
                st.sidebar.success(f"Loaded batch: {len(mols)} molecules ({fails} fails).")
            else:
                st.sidebar.error("Failed to parse SMILES dataset.")

# Update active molecule if a new one was successfully parsed
if new_mol is not None:
    st.session_state.active_mol = new_mol
    st.session_state.active_name = new_name

# Active Molecule Sidebar Preview
if st.session_state.active_mol is not None:
    st.sidebar.markdown("---")
    st.sidebar.markdown(f"**Currently Active:** `{st.session_state.active_name}`")
    try:
        svg_prev = structure_module.draw_molecule_svg(st.session_state.active_mol, width=220, height=180)
        st.sidebar.markdown(
            f'<div style="background: rgba(15,23,42,0.4); border-radius: 8px; padding: 5px; display: flex; justify-content: center;">{svg_prev}</div>',
            unsafe_allow_html=True
        )
    except Exception:
        pass

# ----------------- MAIN INTERFACE -----------------
st.markdown(f'<h1 class="gradient-title">RDKit Chemical Studio</h1>', unsafe_allow_html=True)
st.markdown("An interactive research dashboard containing chemical structures, transformations, fingerprints, and reactions.")

# Main navigation tabs
tabs = st.tabs([
    "🔍 Active Molecule (I/O & 3D)",
    "📦 Batch Dataset (SDF/SMI)",
    "🎯 Substructure Search",
    "✂️ Edits & Scaffolds",
    "🧬 Fingerprints & Similarity",
    "📈 Descriptors & Charges",
    "⚗️ Reactions Simulator"
])

# ----------------- TAB 1: SINGLE MOLECULE VIEW & 3D -----------------
with tabs[0]:
    card_header("Active Molecule Viewer", "Inspect structure properties, 2D coordinates, and 3D conformers", "Single Mode")
    
    col_draw, col_info = st.columns([1, 1])
    
    with col_draw:
        st.markdown("**2D Structure Depiction**")
        
        # Rendering toggles
        c1, c2, c3 = st.columns(3)
        add_hs = c1.checkbox("Add Hydrogens", value=False, key="toggle_h_main")
        add_indices = c2.checkbox("Show Atom Indices", value=False)
        add_stereo = c3.checkbox("Show Stereochemistry", value=True)
        
        # Prepare molecule based on H-toggle
        m_display = structure_module.toggle_hydrogens(st.session_state.active_mol, add_hs)
        
        # Render SVG
        svg_code = structure_module.draw_molecule_svg(
            m_display, width=500, height=400, 
            add_indices=add_indices, add_stereo=add_stereo
        )
        st.markdown(f'<div style="background: rgba(15, 23, 42, 0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_code}</div>', unsafe_allow_html=True)
        
        # Download buttons
        s1, s2 = st.columns(2)
        s1.download_button(
            label="Download SVG Image",
            data=svg_code,
            file_name=f"{st.session_state.active_name}.svg",
            mime="image/svg+xml"
        )
        molblock_text = io_module.mol_to_molblock(st.session_state.active_mol, st.session_state.active_name)
        s2.download_button(
            label="Download MDL Mol File",
            data=molblock_text,
            file_name=f"{st.session_state.active_name}.mol",
            mime="chemical/x-mdl-molfile"
        )
        
    with col_info:
        st.markdown("**Interactive 3D Conformer Viewer**")
        style_3d = st.selectbox("3D Rendering Style", ["stick", "sphere", "line"])
        
        # Embed conformer
        mol_3d = structure_module.generate_3d_conformer(st.session_state.active_mol)
        if mol_3d:
            structure_module.render_3dmol_viewer(mol_3d, height=330, style=style_3d)
            pdb_block = Chem.MolToPDBBlock(mol_3d)
            st.download_button(
                label="Download 3D Conformer PDB",
                data=pdb_block,
                file_name=f"{st.session_state.active_name}_3D.pdb",
                mime="chemical/x-pdb"
            )
        else:
            st.warning("Failed to generate 3D conformer. Check that the molecule is structurally valid.")
            
    # PNG Metadata Reader Section
    st.markdown("---")
    st.subheader("Extract Molecule from PNG Metadata")
    st.write("RDKit automatically embeds molecular metadata in drawing PNGs. Upload a PNG molecule to extract its structure:")
    png_file = st.file_uploader("Upload RDKit PNG Image", type=["png"], key="png_meta_uploader")
    if png_file is not None:
        extracted_mol = structure_module.mol_from_png_metadata(png_file.read())
        if extracted_mol:
            st.success("Extracted successfully!")
            if st.button("Set as Active Molecule"):
                st.session_state.active_mol = extracted_mol
                st.session_state.active_name = png_file.name.replace(".png", "")
                st.rerun()
        else:
            st.error("No embedded RDKit molecular metadata found in this PNG.")

    # Atoms and Bonds DataFrames
    st.markdown("---")
    st.subheader("Detailed Structural Inspection")
    tab_atoms, tab_bonds = st.tabs(["Atoms Details Table", "Bonds Details Table"])
    with tab_atoms:
        df_atoms = structure_module.get_atoms_dataframe(st.session_state.active_mol)
        if not df_atoms.empty:
            st.dataframe(df_atoms, use_container_width=True)
            st.download_button(
                "Download Atoms CSV", 
                df_atoms.to_csv(index=False), 
                f"{st.session_state.active_name}_atoms.csv", 
                "text/csv"
            )
    with tab_bonds:
        df_bonds = structure_module.get_bonds_dataframe(st.session_state.active_mol)
        if not df_bonds.empty:
            st.dataframe(df_bonds, use_container_width=True)
            st.download_button(
                "Download Bonds CSV", 
                df_bonds.to_csv(index=False), 
                f"{st.session_state.active_name}_bonds.csv", 
                "text/csv"
            )

# ----------------- TAB 2: BATCH DATASET MANAGEMENT -----------------
with tabs[1]:
    card_header("Dataset Batch Manager", "Inspect currently loaded database sets and run diversity picking", "Batch Mode")
    
    if st.session_state.batch_mols:
        # Show stats
        num_mols = len(st.session_state.batch_mols)
        st.write(f"Currently loaded batch contains **{num_mols}** molecules.")
        
        # Display dataset table
        metadata_df = io_module.get_mols_metadata_df(st.session_state.batch_mols, st.session_state.batch_names)
        st.dataframe(metadata_df, use_container_width=True)
        
        col_down1, col_down2, col_down3 = st.columns(3)
        
        # Downloads
        sdf_bytes = io_module.mols_to_sdf_bytes(st.session_state.batch_mols)
        col_down1.download_button("Download Dataset as SDF File", sdf_bytes, "dataset.sdf", "chemical/x-mdl-sdfile")
        
        smi_bytes = io_module.mols_to_smiles_bytes(st.session_state.batch_mols)
        col_down2.download_button("Download Dataset as SMILES File", smi_bytes, "dataset.smi", "text/plain")
        
        col_down3.download_button("Download Metadata Table as CSV", metadata_df.to_csv(index=False), "dataset_metadata.csv", "text/csv")
        
        # Diverse Picking Section
        st.markdown("---")
        st.subheader("Diverse Molecule Picking (MaxMinPicker)")
        st.write("Uses circular Morgan fingerprints to identify a subset of molecules with the maximum diversity:")
        
        pick_count = st.number_input("Number of diverse molecules to select", min_value=1, max_value=num_mols, value=min(3, num_mols))
        if st.button("Run Diversity Picker"):
            picked_indices = fingerprint_module.pick_diverse_mols(st.session_state.batch_mols, pick_count)
            picked_mols = [st.session_state.batch_mols[i] for i in picked_indices]
            picked_names = [st.session_state.batch_names[i] for i in picked_indices]
            
            st.success(f"Picked {len(picked_indices)} diverse molecules:")
            
            # Show grid of picked molecules
            cols = st.columns(len(picked_indices))
            for i, (m, name) in enumerate(zip(picked_mols, picked_names)):
                cols[i].write(f"**{name}** (Index {picked_indices[i]})")
                svg_pick = structure_module.draw_molecule_svg(m, width=200, height=180)
                cols[i].markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 8px; padding: 5px; display: flex; justify-content: center;">{svg_pick}</div>', unsafe_allow_html=True)
                
            # Allow saving this subset
            if st.button("Set Selection as Active Dataset"):
                st.session_state.batch_mols = picked_mols
                st.session_state.batch_names = picked_names
                st.rerun()
    else:
        st.info("No batch dataset loaded. Upload an SDF or SMILES file in the sidebar to begin.")

# ----------------- TAB 3: SUBSTRUCTURE SEARCHING -----------------
with tabs[2]:
    card_header("Substructure Search Engine", "Query target molecules using SMILES or SMARTS", "Searching")
    
    col_q, col_opt = st.columns([1.5, 1])
    with col_q:
        q_str = st.text_input("Enter Query (SMILES or SMARTS)", "ccO") # Phenol-like oxygen on aromatic ring
        q_type = st.radio("Query Language", ["SMARTS", "SMILES"], horizontal=True)
    with col_opt:
        use_chirality = st.checkbox("Match Stereochemistry (useChirality)", value=False)
        constraint_toggle = st.checkbox("Apply Sidechain Constraints", value=False)
        constraint_type = st.selectbox("Sidechain Constraint Type", ["alkyl", "all_carbon", "no_nitrogen"], disabled=not constraint_toggle)
        
    st.markdown("---")
    
    col_active, col_dataset = st.columns([1, 1.2])
    
    with col_active:
        st.subheader("Match Active Molecule")
        has_match, matches, query, error_msg = substructure_module.find_substructure_matches(
            st.session_state.active_mol, q_str, q_type, use_chirality
        )
        
        if error_msg:
            st.error(error_msg)
        elif has_match:
            # Apply checker if checked
            if constraint_toggle:
                matches = substructure_module.filter_matches_with_checker(
                    st.session_state.active_mol, query, matches, constraint_type
                )
                has_match = len(matches) > 0
                
            if has_match:
                st.success(f"Matched! Found **{len(matches)}** substructure occurrence(s).")
                # Highlight atoms and bonds
                hatoms, hbonds = substructure_module.get_highlight_atoms_and_bonds(st.session_state.active_mol, query, matches)
                
                # Option to select highlight color
                hl_color = st.color_picker("Highlight Color", "#ef4444")
                # Convert hex to RGB tuple (0.0 to 1.0)
                hex_color = hl_color.lstrip('#')
                rgb_color = tuple(int(hex_color[i:i+2], 16)/255.0 for i in (0, 2, 4))
                
                color_map = {idx: rgb_color for idx in hatoms}
                bond_map = {idx: rgb_color for idx in hbonds}
                
                # Show highlighted molecule
                svg_match = structure_module.draw_molecule_svg(
                    st.session_state.active_mol, width=450, height=350,
                    highlight_atoms=hatoms, highlight_atom_colors=color_map,
                    highlight_bonds=hbonds, highlight_bond_colors=bond_map
                )
                st.markdown(f'<div style="background: rgba(15, 23, 42, 0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_match}</div>', unsafe_allow_html=True)
                
                # Atom mapping indices
                maps = substructure_module.get_atom_map_indices(query)
                if maps:
                    st.write("Atom Map indices in SMARTS:")
                    st.json(maps)
            else:
                st.warning("Query matched structure, but failed custom sidechain constraints.")
        else:
            st.warning("No substructure match found for active molecule.")
            
    with col_dataset:
        st.subheader("Search Batch Dataset")
        if st.session_state.batch_mols and query is not None:
            matched_batch = []
            matched_names = []
            
            for m, name in zip(st.session_state.batch_mols, st.session_state.batch_names):
                b_match, b_indices, _, _ = substructure_module.find_substructure_matches(m, q_str, q_type, use_chirality)
                if b_match:
                    if constraint_toggle:
                        b_indices = substructure_module.filter_matches_with_checker(m, query, b_indices, constraint_type)
                    if b_indices:
                        matched_batch.append(m)
                        matched_names.append(name)
                        
            st.write(f"Found **{len(matched_batch)}** matches in dataset.")
            
            if matched_batch:
                # Render matched molecules in grid
                grid_cols = st.columns(2)
                for idx, (m, name) in enumerate(zip(matched_batch[:6], matched_names[:6])): # limit to first 6
                    col = grid_cols[idx % 2]
                    col.write(f"**{name}**")
                    
                    b_match, b_indices, _, _ = substructure_module.find_substructure_matches(m, q_str, q_type, use_chirality)
                    hatoms, hbonds = substructure_module.get_highlight_atoms_and_bonds(m, query, b_indices)
                    
                    svg_m = structure_module.draw_molecule_svg(m, width=220, height=180, highlight_atoms=hatoms, highlight_bonds=hbonds)
                    col.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 8px; padding: 5px; display: flex; justify-content: center;">{svg_m}</div>', unsafe_allow_html=True)
                
                if len(matched_batch) > 6:
                    st.write(f"... and {len(matched_batch) - 6} more molecules.")
                    
                # Download button for matched SDF
                sdf_matched = io_module.mols_to_sdf_bytes(matched_batch)
                st.download_button("Download Matched Molecules (SDF)", sdf_matched, "matched_results.sdf", "chemical/x-mdl-sdfile")
        else:
            st.info("Load a query and batch dataset to search across multiple molecules.")

# ----------------- TAB 4: TRANSFORMATIONS & DECOMPOSITION -----------------
with tabs[3]:
    card_header("Chemical Transformations & Scaffolds", "Modify chemical structures or extract molecular fragments", "Edit & Decomposition")
    
    st.subheader("1. Substructure Edits")
    edit_mode = st.selectbox("Select Edit Action", ["Delete Substructure", "Replace Substructure", "Replace Core (Isolate Sidechains)"])
    
    col_edit1, col_edit2 = st.columns([1, 1])
    
    with col_edit1:
        st.write("**Original Structure**")
        orig_svg = structure_module.draw_molecule_svg(st.session_state.active_mol, width=400, height=300)
        st.markdown(f'<div style="background: rgba(15, 23, 42, 0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{orig_svg}</div>', unsafe_allow_html=True)
        
    with col_edit2:
        st.write("**Transformed Structure**")
        prod_mol = None
        error_msg = ""
        
        if edit_mode == "Delete Substructure":
            del_query = st.text_input("Enter Query SMARTS to Delete", "C(=O)O") # Delete carboxyl acid groups by default
            if del_query:
                prod_mol, error_msg = transformation_module.delete_substructure(st.session_state.active_mol, del_query)
                
        elif edit_mode == "Replace Substructure":
            rep_query = st.text_input("Query SMARTS to Replace", "C(=O)O")
            rep_with = st.text_input("Replacement SMILES Group", "N") # replace acid with amine by default (make amide/amine)
            if rep_query and rep_with:
                prod_mol, error_msg = transformation_module.replace_substructure(st.session_state.active_mol, rep_query, rep_with)
                
        elif edit_mode == "Replace Core (Isolate Sidechains)":
            core_query = st.text_input("Core SMARTS to Cut Out", "c1ccccc1") # cut benzene ring core by default
            if core_query:
                prod_mol, error_msg = transformation_module.replace_core(st.session_state.active_mol, core_query)
                
        if error_msg:
            st.error(error_msg)
        elif prod_mol is not None:
            prod_svg = structure_module.draw_molecule_svg(prod_mol, width=400, height=300)
            st.markdown(f'<div style="background: rgba(15, 23, 42, 0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{prod_svg}</div>', unsafe_allow_html=True)
            
            # Action to save output
            if st.button("Set Transformed as Active Molecule"):
                st.session_state.active_mol = prod_mol
                st.session_state.active_name = f"{st.session_state.active_name}_edit"
                st.rerun()
        else:
            st.warning("Awaiting query/replacement inputs.")
            
    # Fragmentation Engines Section
    st.markdown("---")
    st.subheader("2. Fragmentation Scaffolds (Bemis-Murcko, BRICS, RECAP)")
    
    frag_type = st.radio("Decomposition Type", ["Bemis-Murcko Scaffold", "BRICS Decomposition", "RECAP Deconstruction"], horizontal=True)
    
    if st.button("Run Decomposition"):
        if frag_type == "Bemis-Murcko Scaffold":
            scaff, err = transformation_module.get_murcko_scaffold(st.session_state.active_mol)
            if err:
                st.error(err)
            elif scaff is not None:
                st.success("Calculated Bemis-Murcko Scaffold:")
                col_s1, col_s2 = st.columns(2)
                col_s1.write(f"Scaffold SMILES: `{Chem.MolToSmiles(scaff)}`")
                svg_scaff = structure_module.draw_molecule_svg(scaff, width=350, height=280)
                col_s2.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_scaff}</div>', unsafe_allow_html=True)
                
                # Download
                scaff_block = io_module.mol_to_molblock(scaff)
                st.download_button("Download Scaffold (MOL)", scaff_block, "scaffold.mol", "chemical/x-mdl-molfile")
                
        else:
            if frag_type == "BRICS Decomposition":
                frags, err = transformation_module.decompose_brics(st.session_state.active_mol)
                label = "BRICS"
            else:
                frags, err = transformation_module.decompose_recap(st.session_state.active_mol)
                label = "RECAP"
                
            if err:
                st.error(err)
            elif frags:
                st.success(f"Generated {len(frags)} fragments via {label}:")
                frag_cols = st.columns(min(len(frags), 4))
                for i, f in enumerate(frags):
                    col = frag_cols[i % 4]
                    col.write(f"Fragment {i+1}: `{Chem.MolToSmiles(f)}`")
                    svg_f = structure_module.draw_molecule_svg(f, width=200, height=180)
                    col.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 8px; padding: 5px; display: flex; justify-content: center;">{svg_f}</div>', unsafe_allow_html=True)
                
                # Download
                sdf_frags = io_module.mols_to_sdf_bytes(frags)
                st.download_button(f"Download {label} Fragments (SDF)", sdf_frags, "fragments.sdf", "chemical/x-mdl-sdfile")
            else:
                st.warning("Molecule could not be fragmented using these rules.")
                
    # Maximum Common Substructure Section (MCS)
    st.markdown("---")
    st.subheader("3. Maximum Common Substructure (MCS)")
    st.write("Run Maximum Common Substructure on the batch dataset:")
    
    col_mcs_opt1, col_mcs_opt2 = st.columns(2)
    atom_compare = col_mcs_opt1.selectbox("Atom Match Rule", ["CompareElements", "CompareAny", "CompareAnyHeavyAtom"])
    bond_compare = col_mcs_opt2.selectbox("Bond Match Rule", ["CompareOrder", "CompareAny", "CompareOrderExact"])
    
    if st.button("Compute MCS"):
        if st.session_state.batch_mols and len(st.session_state.batch_mols) >= 2:
            mcs_res, err = transformation_module.find_maximum_common_substructure(
                st.session_state.batch_mols, atom_compare, bond_compare
            )
            if err:
                st.error(err)
            elif mcs_res:
                st.success(f"MCS Found! SMARTS: `{mcs_res['smarts']}`")
                col_m1, col_m2 = st.columns([1, 1.5])
                with col_m1:
                    st.write(f"**MCS Structure ({mcs_res['num_atoms']} atoms, {mcs_res['num_bonds']} bonds)**")
                    svg_mcs = structure_module.draw_molecule_svg(mcs_res["mol"], width=300, height=250)
                    st.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_mcs}</div>', unsafe_allow_html=True)
                    
                with col_m2:
                    st.write("**Highlighted in Dataset**")
                    # Draw first 3 dataset molecules with MCS highlighted
                    mcs_cols = st.columns(3)
                    for i, (m, name) in enumerate(zip(st.session_state.batch_mols[:3], st.session_state.batch_names[:3])):
                        # Get matching atom indices
                        has_match, b_indices, _, _ = substructure_module.find_substructure_matches(m, mcs_res["smarts"], "SMARTS")
                        hatoms, hbonds = substructure_module.get_highlight_atoms_and_bonds(m, mcs_res["mol"], b_indices)
                        
                        mcs_cols[i].write(f"**{name}**")
                        svg_h_mcs = structure_module.draw_molecule_svg(m, width=200, height=180, highlight_atoms=hatoms, highlight_bonds=hbonds)
                        mcs_cols[i].markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 8px; padding: 5px; display: flex; justify-content: center;">{svg_h_mcs}</div>', unsafe_allow_html=True)
        else:
            st.warning("Please upload a batch dataset containing at least 2 molecules in the sidebar.")

# ----------------- TAB 5: FINGERPRINTS & MOLECULAR SIMILARITY -----------------
with tabs[4]:
    card_header("Fingerprinting & Similarity Maps", "Generate binary descriptors, inspect bit details and visual similarity gradients", "Fingerprints")
    
    fp_select = st.selectbox("Select Fingerprint Type", [
        "Morgan (Circular)", "RDKit (Topological)", "Atom Pairs", "Topological Torsions", "MACCS Keys"
    ])
    
    st.markdown("---")
    
    col_sim_map, col_bit_env = st.columns([1.2, 1])
    
    with col_sim_map:
        st.subheader("1. Molecular Similarity Calculator & Maps")
        comp_smi = st.text_input("Enter Secondary Molecule SMILES to Compare", "CC(=O)Oc1ccccc1C(=O)O") # Default Aspirin
        
        sim_metric = st.selectbox("Similarity Metric", ["Tanimoto", "Dice", "Cosine", "Sokal", "Kulczynski", "Tversky"])
        
        # Tversky weights
        t_alpha, t_beta = 0.5, 0.5
        if sim_metric == "Tversky":
            t_col1, t_col2 = st.columns(2)
            t_alpha = t_col1.slider("Tversky Alpha (weight reference)", 0.0, 1.0, 0.5)
            t_beta = t_col2.slider("Tversky Beta (weight probe)", 0.0, 1.0, 0.5)
            
        comp_mol = Chem.MolFromSmiles(comp_smi)
        
        if comp_mol is not None:
            sim_score = fingerprint_module.calculate_similarity(
                st.session_state.active_mol, comp_mol, fp_select, sim_metric,
                radius=2, n_bits=2048, tversky_alpha=t_alpha, tversky_beta=t_beta
            )
            
            metric_card(f"{sim_metric} Similarity Score", f"{sim_score:.4f}")
            
            # Similarity Map Drawing
            st.write("**Visual Similarity Heatmap (probe vs active reference)**")
            st.write("Displays atom contribution gradients (Green matches reference, Orange differs):")
            
            # Map type (use Morgan or topological RDKit)
            map_fp = "Morgan"
            if "RDKit" in fp_select:
                map_fp = "RDKit"
            elif "Atom Pairs" in fp_select:
                map_fp = "Atom Pairs"
            elif "Torsions" in fp_select:
                map_fp = "Topological Torsions"
                
            svg_sim_map = fingerprint_module.generate_similarity_map_svg(
                st.session_state.active_mol, comp_mol, fp_type=map_fp
            )
            st.markdown(f'<div style="background: rgba(15, 23, 42, 0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_sim_map}</div>', unsafe_allow_html=True)
            
            st.download_button("Download Similarity Map SVG", svg_sim_map, "similarity_map.svg", "image/svg+xml")
        else:
            st.error("Invalid comparison molecule SMILES.")
            
    with col_bit_env:
        st.subheader("2. Fingerprint Bit Environment Visualizer")
        st.write("Inspect which exact molecular environments are captured by individual fingerprint bits:")
        
        if "Morgan" in fp_select:
            radius_sel = st.slider("Morgan Radius", 1, 3, 2)
            bit_info, fp_vect = fingerprint_module.get_morgan_bit_info(st.session_state.active_mol, radius=radius_sel)
            
            active_bits = list(bit_info.keys())
            st.write(f"Total active bits: **{len(active_bits)}**")
            
            bit_id_sel = st.selectbox("Select Active Bit to Visualize", sorted(active_bits))
            
            if bit_id_sel:
                # Render Morgan bit
                svg_bit = fingerprint_module.draw_morgan_bit_env(st.session_state.active_mol, bit_id_sel, bit_info)
                if svg_bit:
                    st.write(f"Sub-environment for Bit **{bit_id_sel}**:")
                    st.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_bit}</div>', unsafe_allow_html=True)
                    st.download_button("Download Bit SVG", svg_bit, f"bit_{bit_id_sel}.svg", "image/svg+xml")
                else:
                    st.warning("Could not render this bit sub-environment.")
        elif "RDKit" in fp_select:
            bit_info, fp_vect = fingerprint_module.get_rdkit_bit_info(st.session_state.active_mol)
            active_bits = list(bit_info.keys())
            st.write(f"Total active bits: **{len(active_bits)}**")
            bit_id_sel = st.selectbox("Select Active Bit to Visualize", sorted(active_bits))
            
            if bit_id_sel:
                svg_bit = fingerprint_module.draw_rdkit_bit_env(st.session_state.active_mol, bit_id_sel, bit_info)
                if svg_bit:
                    st.write(f"Sub-environment for Bit **{bit_id_sel}**:")
                    st.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_bit}</div>', unsafe_allow_html=True)
                    st.download_button("Download Bit SVG", svg_bit, f"bit_{bit_id_sel}.svg", "image/svg+xml")
                else:
                    st.warning("Could not render this bit sub-environment.")
        else:
            st.info("Sub-environment bit visualization is supported for Morgan (Circular) and RDKit (Topological) fingerprints.")

# ----------------- TAB 6: DESCRIPTORS & CHARGE MAPS -----------------
with tabs[5]:
    card_header("Molecular Descriptors & Partial Charges", "Calculate bulk physical properties and map Gasteiger partial charges", "Descriptors & Charges")
    
    col_desc, col_chg = st.columns([1, 1.2])
    
    with col_desc:
        st.subheader("1. Physical & Chemical Descriptors")
        df_desc = descriptor_module.get_descriptors_dataframe(st.session_state.active_mol)
        if not df_desc.empty:
            # Display properties cleanly
            st.dataframe(df_desc, use_container_width=True, height=360)
            
            st.download_button(
                "Download Descriptors CSV", 
                df_desc.to_csv(index=False), 
                f"{st.session_state.active_name}_descriptors.csv", 
                "text/csv"
            )
            
            # Plotly Radar Chart of drug likeness
            st.markdown("---")
            st.write("**Lipinski / Drug-likeness Radar Profile**")
            fig_radar = descriptor_module.get_radar_chart_plotly(st.session_state.active_mol)
            if fig_radar:
                st.plotly_chart(fig_radar, use_container_width=True)
                
    with col_chg:
        st.subheader("2. Gasteiger Partial Charge Map")
        st.write("Color-coded drawing representing atom partial charges (Red = Negative, Blue = Positive):")
        
        # Render Gasteiger charge SVG
        svg_charge = descriptor_module.draw_charge_distribution_svg(st.session_state.active_mol, width=500, height=420)
        st.markdown(f'<div style="background: rgba(15, 23, 42, 0.4); border-radius: 12px; padding: 10px; display: flex; justify-content: center;">{svg_charge}</div>', unsafe_allow_html=True)
        
        st.download_button("Download Charge Map SVG", svg_charge, "charge_map.svg", "image/svg+xml")

# ----------------- TAB 7: REACTIONS SIMULATOR -----------------
with tabs[6]:
    card_header("Chemical Reactions Simulator", "Define synthetic pathways using Reaction SMARTS, protect groups, and generate products", "Synthesis & Simulation")
    
    # Reaction presets
    presets = {
        "Amide Bond Formation": "[C:1](=[O:2])[O:3].[N:4]>>[C:1](=[O:2])[N:4]",
        "Click Triazole Coupling": "[C:1]#[C:2].[N:3]=[N+:4]=[N-:5]>>[c:1]1[c:2][n:3][n:4][n:5]1",
        "Esterification": "[C:1](=[O:2])[O:3].[O:4]>>[C:1](=[O:2])[O:4]",
        "Suzuki Coupling (Aromatic C-C)": "[c:1][B:2]([O:3])[O:4].[c:5][Cl,Br,I:6]>>[c:1][c:5]"
    }
    
    st.write("Select a preset reaction or enter a custom Reaction SMARTS:")
    preset_choice = st.selectbox("Reaction Presets", list(presets.keys()) + ["Custom Reaction SMARTS"])
    
    if preset_choice == "Custom Reaction SMARTS":
        rxn_smarts = st.text_input("Reaction SMARTS", "[C:1](=[O:2])[O:3].[N:4]>>[C:1](=[O:2])[N:4]")
    else:
        rxn_smarts = presets[preset_choice]
        st.info(f"Loaded SMARTS: `{rxn_smarts}`")
        
    rxn, rxn_err = reaction_module.create_reaction(rxn_smarts)
    
    if rxn_err:
        st.error(rxn_err)
    elif rxn:
        # Draw reaction
        reaction_svg = reaction_module.draw_reaction_visualization(rxn)
        st.markdown(reaction_svg, unsafe_allow_html=True)
        
        # Reactants input
        st.markdown("---")
        st.subheader("Reactant Inputs")
        num_reactants = rxn.GetNumReactantTemplates()
        st.write(f"This reaction requires **{num_reactants}** reactant(s).")
        
        reactants_inputs = []
        # Defaults based on preset
        defaults = ["CC(=O)O", "CCN", "C"]
        if "Click" in preset_choice:
            defaults = ["C#CC", "CCN=[N+]=[N-]", "C"]
        elif "Suzuki" in preset_choice:
            defaults = ["OB(O)c1ccccc1", "Ic1ccccc1", "C"]
            
        for i in range(num_reactants):
            def_val = defaults[i] if i < len(defaults) else "C"
            reactants_inputs.append(st.text_input(f"Reactant {i+1} SMILES", value=def_val, key=f"reactant_{i}"))
            
        # Parse inputs
        reactants_mols = []
        ready = True
        for i, ri in enumerate(reactants_inputs):
            m = Chem.MolFromSmiles(ri)
            if m is None:
                st.error(f"Reactant {i+1} is invalid.")
                ready = False
            else:
                reactants_mols.append(m)
                
        # Advanced Atom Protection Section
        st.markdown("---")
        st.subheader("Advanced Reaction Controls (Atom Protection)")
        st.write("Specify a SMARTS pattern of functional groups to protect (respects RDKit `_protected` flags):")
        protect_smarts = st.text_input("SMARTS to Protect (e.g. '[N;$(NC=[O,S])]' protects amide nitrogens from coupling)", "")
        
        if ready and st.button("Run Reaction Simulator"):
            prod_sets, err = reaction_module.run_reaction(rxn, reactants_mols, protect_smarts)
            
            if err:
                st.error(err)
            elif prod_sets:
                st.success(f"Reaction succeeded! Generated **{len(prod_sets)}** possible product configuration(s).")
                
                for s_idx, p_set in enumerate(prod_sets):
                    st.markdown(f"**Product Set {s_idx + 1}:**")
                    p_cols = st.columns(min(len(p_set), 4))
                    for p_idx, p in enumerate(p_set):
                        col = p_cols[p_idx % 4]
                        smi_p = Chem.MolToSmiles(p)
                        col.write(f"Product {p_idx+1}: `{smi_p}`")
                        svg_p = structure_module.draw_molecule_svg(p, width=200, height=180)
                        col.markdown(f'<div style="background: rgba(15,23,42,0.4); border-radius: 8px; padding: 5px; display: flex; justify-content: center;">{svg_p}</div>', unsafe_allow_html=True)
                        
                # Download button for all generated products
                all_prods = []
                for p_set in prod_sets:
                    for p in p_set:
                        all_prods.append(p)
                sdf_prods = io_module.mols_to_sdf_bytes(all_prods)
                st.download_button("Download All Generated Products (SDF)", sdf_prods, "products.sdf", "chemical/x-mdl-sdfile")
            else:
                st.warning("No reaction took place. Make sure the reactants match the reaction templates.")
