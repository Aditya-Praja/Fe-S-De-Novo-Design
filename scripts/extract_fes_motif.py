from pathlib import Path
import argparse
import json

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

def atom_dist(atom1, atom2):
    return float(np.linalg.norm(atom1.coord - atom2.coord))
    

def find_2fe2s_clusters(structures):
    clusters = []
    
    for residue in structures.get_residues():
        if residue.resname.strip().upper() != "FES":
            continue
        
        atoms = list(residue.get_atoms())
        fe_atoms = [
            atom for atom in atoms
            if atom.element.strip().upper() == "FE"
        ]
        
        s_atoms = [
            atom for atom in atoms
            if atom.element.strip().upper() == "S"
        ]
        
        if len(fe_atoms) == 2 and len(s_atoms) == 2:
            clusters.append((residue, fe_atoms, s_atoms))
            
    return clusters
    

def find_coordinating_cys(structure, fe_atoms, cutoff=2.8):
    coordinating = []
    
    for residue in structure.get_residues():
        if residue.resname.strip().upper() != "CYS":
            continue
        
        if "SG" not in residue:
            continue
        
        sg = residue["SG"]
        dists = [atom_dist(sg, fe) for fe in fe_atoms]
        
        if min(dists) <= cutoff:
            coordinating.append({
                "residue": residue,
                "sg_atom": sg,
                "fe_distances": dists,
                "nearest_fe_index": int(np.argmin(dists))
            })

    return coordinating
    
    
class MotifSelect(Select):

    def __init__(self, selected_residues):
        self.selected_ids = {id(r) for r in selected_residues}

    def accept_model(self, model):
        return model.id == 0

    def accept_residue(self, residue):
        return id(residue) in self.selected_ids


def save_selected_pdb(structure, residues, output_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), select=gMotifSelect(residues))


def residue_label(residue):
    chain = residue.get_parent()
    return {
        "chain": chain.id,
        "name": residue.resname,
        "num": residue.id[1],
        "insertion": residue.id[2].strip()
    }
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, help="Input PDB file", required=True)
    parser.add_argument("--output", type=Path, help="Output PDB file", required=True)
    parser.add_argument("--cutoff", type=float, default=2.8)
    parser.add_argument("--cluster", type=int, help="Cluster index", default=0)
    args = parser.parse_args()
    
    if args.cutoff <= 0:
        parser.error("--cutoff must be positive")
        
    args.output.mkdir(parents=True, exist_ok=True)
    
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("reference", str(args.input))
    structure = structure[0]
    
    clusters = find_2fe2s_clusters(structure)
    print("Found clusters:", len(clusters))
    
    if not clusters:
        raise ValueError("No clusters found")
    
    if not 0 <= args.cluster < len(clusters):
        parser.error("Invalid cluster index")
        
    cluster, fe_atoms, s_atoms = clusters[args.cluster]
    
    ligands = find_coordinating_cys(structure, fe_atoms, cutoff=args.cutoff)
    
    print("Found ligands:", len(ligands))
    
    save_selected_pdb(
        structure, 
        [cluster],
        args.output / "cluster.pdb"
    )
    
    print("Saving cluster pdb")
    
    motif_residues = [cluster] + [item["residue"] for item in ligands]
    save_selected_pdb(
        structure, 
        motif_residues, 
        args.output / "motif.pdb"
    )
    
    print("Saving motif pdb")
    
    metadata = {
        "source_pdb": args.input.name,
        "model_index": 0,
        "cluster": residue_label(cluster),
        "cluster_formula": "Fe2S2",
        "cutoff_angstrom": args.cutoff,
        "cluster_atoms": [
            {
                "name": atom.name,
                "element": atom.element,
                "coordinates": atom.coord.tolist(),
            }
            for atom in fe_atoms + s_atoms
        ],
        "coordinating_residues": [
            {
                **residue_label(item["residue"]),
                "sg_coordinates": item["sg_atom"].coord.tolist(),
                "distances_to_fe": item["fe_distances"],
                "nearest_fe_index": item["nearest_fe_index"],
            }
            for item in ligands
        ],
    }
    
    with open(args.output / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
        
    if len(ligands) != 4:
        print("Warning: not 4 ligands")
        
    print(f"Saved extracted motif to {args.output}")
    
if __name__ == "__main__":
    main()
    
        
    