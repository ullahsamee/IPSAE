"""
Core functionality for ipSAE calculations.
"""

import sys
import os
import math
import json
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Set
from pathlib import Path

np.set_printoptions(threshold=np.inf)


class IpsaeCalculator:
    """Main class for calculating ipSAE scores."""
    
    def __init__(self, pae_cutoff: float = 10.0, dist_cutoff: float = 10.0):
        """
        Initialize the calculator.
        
        Args:
            pae_cutoff: PAE cutoff value
            dist_cutoff: Distance cutoff value
        """
        self.pae_cutoff = pae_cutoff
        self.dist_cutoff = dist_cutoff
        
    def calculate(self, pae_file: Union[str, Path], structure_file: Union[str, Path], 
                  output_dir: Optional[Union[str, Path]] = None) -> Dict:
        """
        Calculate ipSAE scores for the given files.
        
        Args:
            pae_file: Path to PAE file (JSON or NPZ)
            structure_file: Path to structure file (PDB or CIF)
            output_dir: Output directory (defaults to structure file directory)
            
        Returns:
            Dictionary containing calculated scores
        """
        # Convert to Path objects
        pae_file = Path(pae_file)
        structure_file = Path(structure_file)
        
        if output_dir is None:
            output_dir = structure_file.parent
        else:
            output_dir = Path(output_dir)
            
        # Determine file types
        file_type = self._determine_file_type(pae_file, structure_file)
        
        # Load and process data
        residues, cb_residues, chains, token_mask = self._load_structure(structure_file, file_type)
        pae_matrix, plddt, cb_plddt, global_scores = self._load_pae_data(pae_file, file_type, token_mask, len(residues), residues)
        
        # Calculate scores
        results = self._calculate_scores(residues, cb_residues, chains, pae_matrix, plddt, cb_plddt, file_type, global_scores)
        
        # Save outputs
        self._save_outputs(results, output_dir, structure_file.stem)
        
        return results
    
    def _determine_file_type(self, pae_file: Path, structure_file: Path) -> str:
        """Determine the file type (AF2, AF3, or Boltz1)."""
        if structure_file.suffix == ".pdb":
            return "af2"
        elif structure_file.suffix == ".cif" and pae_file.suffix == ".json":
            return "af3"
        elif structure_file.suffix == ".cif" and pae_file.suffix == ".npz":
            return "boltz1"
        else:
            raise ValueError(f"Unsupported file combination: {pae_file.suffix} and {structure_file.suffix}")
    
    def _load_structure(self, structure_file: Path, file_type: str) -> Tuple[List, List, np.ndarray, List]:
        """Load structure file and extract residue information."""
        residues = []
        cb_residues = []
        chains = []
        token_mask = []
        atomsitefield_num = 0
        atomsitefield_dict = {}
        
        residue_set = {"ALA", "ARG", "ASN", "ASP", "CYS",
                      "GLN", "GLU", "GLY", "HIS", "ILE",
                      "LEU", "LYS", "MET", "PHE", "PRO",
                      "SER", "THR", "TRP", "TYR", "VAL",
                      "DA", "DC", "DT", "DG", "A", "C", "U", "G"}

        with open(structure_file, 'r') as f:
            for line in f:
                if line.startswith("_atom_site."):
                    line = line.strip()
                    (atomsite, fieldname) = line.split(".")
                    atomsitefield_dict[fieldname] = atomsitefield_num
                    atomsitefield_num += 1
                    
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if file_type in ["af3", "boltz1"]:  # CIF format
                        atom = self._parse_cif_atom_line(line, atomsitefield_dict)
                    else:  # PDB format
                        atom = self._parse_pdb_atom_line(line)
                        
                    if atom is None:  # ligand atom
                        token_mask.append(0)
                        continue

                    if atom['atom_name'] == "CA" or "C1" in atom['atom_name']:
                        token_mask.append(1)
                        residues.append({
                            'atom_num': atom['atom_num'],
                            'coor': np.array([atom['x'], atom['y'], atom['z']]),
                            'res': atom['residue_name'],
                            'chainid': atom['chain_id'],
                            'resnum': atom['residue_seq_num'],
                            'residue': f"{atom['residue_name']:3}   {atom['chain_id']:3} {atom['residue_seq_num']:4}"
                        })
                        chains.append(atom['chain_id'])

                    if atom['atom_name'] == "CB" or "C3" in atom['atom_name'] or (atom['residue_name']=="GLY" and atom['atom_name']=="CA"):
                        cb_residues.append({
                            'atom_num': atom['atom_num'],
                            'coor': np.array([atom['x'], atom['y'], atom['z']]),
                            'res': atom['residue_name'],
                            'chainid': atom['chain_id'],
                            'resnum': atom['residue_seq_num'],
                            'residue': f"{atom['residue_name']:3}   {atom['chain_id']:3} {atom['residue_seq_num']:4}"
                        })

                    # add nucleic acids and non-CA atoms in PTM residues to tokens
                    if atom['atom_name'] != "CA" and "C1" not in atom['atom_name'] and atom['residue_name'] not in residue_set:
                        token_mask.append(0)
        
        return residues, cb_residues, np.array(chains), token_mask
    
    def _parse_pdb_atom_line(self, line: str) -> Dict:
        """Parse PDB ATOM/HETATM line."""
        atom_num = line[6:11].strip()
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_seq_num = line[22:26].strip()
        x = line[30:38].strip()
        y = line[38:46].strip()
        z = line[46:54].strip()

        # Convert string numbers to integers or floats as appropriate
        atom_num = int(atom_num)
        residue_seq_num = int(residue_seq_num)
        x = float(x)
        y = float(y)
        z = float(z)

        return {
            'atom_num': atom_num,
            'atom_name': atom_name,
            'residue_name': residue_name,
            'chain_id': chain_id,
            'residue_seq_num': residue_seq_num,
            'x': x,
            'y': y,
            'z': z
        }
    
    def _parse_cif_atom_line(self, line: str, fielddict: Dict) -> Optional[Dict]:
        """Parse mmCIF ATOM/HETATM line."""
        linelist = line.split()
        atom_num = linelist[fielddict['id']]
        atom_name = linelist[fielddict['label_atom_id']]
        residue_name = linelist[fielddict['label_comp_id']]
        chain_id = linelist[fielddict['label_asym_id']]
        residue_seq_num = linelist[fielddict['label_seq_id']]
        x = linelist[fielddict['Cartn_x']]
        y = linelist[fielddict['Cartn_y']]
        z = linelist[fielddict['Cartn_z']]

        if residue_seq_num == ".":
            return None  # ligand atom

        # Convert string numbers to integers or floats as appropriate
        atom_num = int(atom_num)
        residue_seq_num = int(residue_seq_num)
        x = float(x)
        y = float(y)
        z = float(z)

        return {
            'atom_num': atom_num,
            'atom_name': atom_name,
            'residue_name': residue_name,
            'chain_id': chain_id,
            'residue_seq_num': residue_seq_num,
            'x': x,
            'y': y,
            'z': z
        }
    
    def _load_pae_data(self, pae_file: Path, file_type: str, token_mask: List, numres: int, residues: List) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
        """Load PAE data and extract plddt and pae_matrix."""
        if file_type == "af2":
            return self._load_af2_data(pae_file, numres)
        elif file_type == "af3":
            return self._load_af3_data(pae_file, token_mask, numres, residues)
        elif file_type == "boltz1":
            return self._load_boltz1_data(pae_file, token_mask, numres)
        else:
            raise ValueError(f"Unknown file type: {file_type}")
    
    def _load_af2_data(self, pae_file: Path, numres: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
        """Load AlphaFold2 PAE data."""
        if not pae_file.exists():
            raise FileNotFoundError(f"AF2 PAE file does not exist: {pae_file}")
            
        if pae_file.suffix == '.pkl':
            data = np.load(pae_file, allow_pickle=True)
        else:
            with open(pae_file, 'r') as file:
                data = json.load(file)
        
        global_scores = {}
        if 'iptm' in data:
            global_scores['iptm'] = float(data['iptm'])
        else:
            global_scores['iptm'] = -1.0
            
        if 'ptm' in data:
            global_scores['ptm'] = float(data['ptm'])
        else:
            global_scores['ptm'] = -1.0
        
        if 'plddt' in data:
            plddt = np.array(data['plddt'])
            cb_plddt = np.array(data['plddt'])
        else:
            plddt = np.zeros(numres)
            cb_plddt = np.zeros(numres)
            
        if 'pae' in data:
            pae_matrix = np.array(data['pae'])
        elif 'predicted_aligned_error' in data:
            pae_matrix = np.array(data['predicted_aligned_error'])
        else:
            raise ValueError("No PAE data found in AF2 file")
            
        return pae_matrix, plddt, cb_plddt, global_scores
    
    def _load_af3_data(self, pae_file: Path, token_mask: List, numres: int, residues: List) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
        """Load AlphaFold3 PAE data."""
        if not pae_file.exists():
            raise FileNotFoundError(f"AF3 PAE file does not exist: {pae_file}")
            
        with open(pae_file, 'r') as file:
            data = json.load(file)
        
        # Extract atom plddts and convert to residue plddts
        atom_plddts = np.array(data['atom_plddts'])
        CA_atom_num = np.array([res['atom_num']-1 for res in residues])
        CB_atom_num = CA_atom_num  # Simplified for now
        
        plddt = atom_plddts[CA_atom_num]
        cb_plddt = atom_plddts[CB_atom_num]
        
        if 'pae' in data:
            pae_matrix_af3 = np.array(data['pae'])
        else:
            raise ValueError("No PAE data found in AF3 file")
        
        # Extract residue-level PAE matrix
        token_array = np.array(token_mask)
        pae_matrix = pae_matrix_af3[np.ix_(token_array.astype(bool), token_array.astype(bool))]
        
        # Load global scores from summary file
        global_scores = self._load_af3_global_scores(pae_file)
        
        return pae_matrix, plddt, cb_plddt, global_scores
    
    def _load_af3_global_scores(self, pae_file: Path) -> Dict:
        """Load AF3 global scores from summary file."""
        global_scores = {}
        summary_file_path = None
        
        if "confidences" in str(pae_file):
            summary_file_path = str(pae_file).replace("confidences", "summary_confidences")
        elif "full_data" in str(pae_file):
            summary_file_path = str(pae_file).replace("full_data", "summary_confidences")

        if summary_file_path and Path(summary_file_path).exists():
            with open(summary_file_path, 'r') as file:
                data_summary = json.load(file)
            global_scores['chain_pair_iptm'] = data_summary.get('chain_pair_iptm', {})
        else:
            global_scores['chain_pair_iptm'] = {}
            
        return global_scores
    
    def _load_boltz1_data(self, pae_file: Path, token_mask: List, numres: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
        """Load Boltz1 PAE data."""
        if not pae_file.exists():
            raise FileNotFoundError(f"Boltz1 PAE file does not exist: {pae_file}")
        
        # Load PAE data
        data_pae = np.load(pae_file)
        pae_matrix_boltz1 = np.array(data_pae['pae'])
        token_array = np.array(token_mask)
        pae_matrix = pae_matrix_boltz1[np.ix_(token_array.astype(bool), token_array.astype(bool))]
        
        # Load pLDDT data
        plddt_file_path = str(pae_file).replace("pae", "plddt")
        if Path(plddt_file_path).exists():
            data_plddt = np.load(plddt_file_path)
            plddt_boltz1 = np.array(100.0 * data_plddt['plddt'])
            plddt = plddt_boltz1[token_array.astype(bool)]
            cb_plddt = plddt_boltz1[token_array.astype(bool)]
        else:
            plddt = np.zeros(np.sum(token_array))
            cb_plddt = np.zeros(np.sum(token_array))
        
        # Load global scores
        global_scores = self._load_boltz1_global_scores(pae_file)
        
        return pae_matrix, plddt, cb_plddt, global_scores
    
    def _load_boltz1_global_scores(self, pae_file: Path) -> Dict:
        """Load Boltz1 global scores from summary file."""
        global_scores = {}
        summary_file_path = str(pae_file).replace("pae", "confidence").replace(".npz", ".json")
        
        if Path(summary_file_path).exists():
            with open(summary_file_path, 'r') as file:
                data_summary = json.load(file)
            global_scores['pair_chains_iptm'] = data_summary.get('pair_chains_iptm', {})
        else:
            global_scores['pair_chains_iptm'] = {}
            
        return global_scores
    
    def _calculate_scores(self, residues: List, cb_residues: List, chains: np.ndarray, 
                         pae_matrix: np.ndarray, plddt: np.ndarray, cb_plddt: np.ndarray, 
                         file_type: str, global_scores: Dict) -> Dict:
        """Calculate all scores including ipSAE, pDockQ, pDockQ2, and LIS."""
        numres = len(residues)
        unique_chains = np.unique(chains)
        residue_types = np.array([res['res'] for res in residues])
        
        # Calculate CB coordinates and distances
        coordinates = np.array([res['coor'] for res in cb_residues])
        distances = np.sqrt(((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :])**2).sum(axis=2))
        
        # Classify chains (protein vs nucleic acid)
        chain_dict = self._classify_chains(chains, residue_types)
        chain_pair_type = self._init_chainpairdict_zeros(unique_chains)
        
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue
                if chain_dict[chain1] == 'nucleic_acid' or chain_dict[chain2] == 'nucleic_acid':
                    chain_pair_type[chain1][chain2] = 'nucleic_acid'
                else:
                    chain_pair_type[chain1][chain2] = 'protein'
        
        # Calculate various scores
        pdockq_scores = self._calculate_pdockq_scores(unique_chains, chains, distances, cb_plddt, pae_matrix, numres)
        lis_scores = self._calculate_lis_scores(unique_chains, chains, pae_matrix)
        ipsae_scores = self._calculate_ipsae_scores(unique_chains, chains, pae_matrix, residues, chain_pair_type, numres, distances)
        
        # Prepare results
        results = {
            'numres': numres,
            'unique_chains': unique_chains,
            'chains': chains,
            'distances': distances,
            'pae_matrix': pae_matrix,
            'plddt': plddt,
            'cb_plddt': cb_plddt,
            'file_type': file_type,
            'residues': residues,
            'pdockq_scores': pdockq_scores,
            'lis_scores': lis_scores,
            'ipsae_scores': ipsae_scores,
            'global_scores': global_scores,
            'chain_dict': chain_dict
        }
        
        return results
    
    def _classify_chains(self, chains: np.ndarray, residue_types: np.ndarray) -> Dict:
        """Classify chains as protein or nucleic acid."""
        nuc_residue_set = {"DA", "DC", "DT", "DG", "A", "C", "U", "G"}
        chain_types = {}
        
        unique_chains = np.unique(chains)
        for chain in unique_chains:
            indices = np.where(chains == chain)[0]
            chain_residues = residue_types[indices]
            nuc_count = sum(residue in nuc_residue_set for residue in chain_residues)
            chain_types[chain] = 'nucleic_acid' if nuc_count > 0 else 'protein'
        
        return chain_types
    
    def _init_chainpairdict_zeros(self, chainlist):
        """Initialize nested dictionary with zeros."""
        return {chain1: {chain2: 0 for chain2 in chainlist if chain1 != chain2} for chain1 in chainlist}
    
    def _init_chainpairdict_set(self, chainlist):
        """Initialize nested dictionary with empty sets."""
        return {chain1: {chain2: set() for chain2 in chainlist if chain1 != chain2} for chain1 in chainlist}
    
    def _calculate_pdockq_scores(self, unique_chains, chains, distances, cb_plddt, pae_matrix, numres):
        """Calculate pDockQ and pDockQ2 scores."""
        pDockQ = self._init_chainpairdict_zeros(unique_chains)
        pDockQ2 = self._init_chainpairdict_zeros(unique_chains)
        pDockQ_unique_residues = self._init_chainpairdict_set(unique_chains)
        pDockQ_cutoff = 8.0

        # Calculate pDockQ
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue
                npairs = 0
                for i in range(numres):
                    if chains[i] != chain1:
                        continue
                    valid_pairs = (chains == chain2) & (distances[i] <= pDockQ_cutoff)
                    npairs += np.sum(valid_pairs)
                    if valid_pairs.any():
                        pDockQ_unique_residues[chain1][chain2].add(i)
                        chain2residues = np.where(valid_pairs)[0]
                        for residue in chain2residues:
                            pDockQ_unique_residues[chain1][chain2].add(residue)

                if npairs > 0:
                    nres = len(list(pDockQ_unique_residues[chain1][chain2]))
                    mean_plddt = cb_plddt[list(pDockQ_unique_residues[chain1][chain2])].mean()
                    x = mean_plddt * math.log10(npairs)
                    pDockQ[chain1][chain2] = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
                else:
                    pDockQ[chain1][chain2] = 0.0

        # Calculate pDockQ2
        ptm_func_vec = np.vectorize(self.ptm_func)
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue
                npairs = 0
                sum_ptm = 0.0
                for i in range(numres):
                    if chains[i] != chain1:
                        continue
                    valid_pairs = (chains == chain2) & (distances[i] <= pDockQ_cutoff)
                    if valid_pairs.any():
                        npairs += np.sum(valid_pairs)
                        pae_list = pae_matrix[i][valid_pairs]
                        pae_list_ptm = ptm_func_vec(pae_list, 10.0)
                        sum_ptm += pae_list_ptm.sum()

                if npairs > 0:
                    nres = len(list(pDockQ_unique_residues[chain1][chain2]))
                    mean_plddt = cb_plddt[list(pDockQ_unique_residues[chain1][chain2])].mean()
                    mean_ptm = sum_ptm / npairs
                    x = mean_plddt * mean_ptm
                    pDockQ2[chain1][chain2] = 1.31 / (1 + math.exp(-0.075 * (x - 84.733))) + 0.005
                else:
                    pDockQ2[chain1][chain2] = 0.0

        return {'pDockQ': pDockQ, 'pDockQ2': pDockQ2, 'unique_residues': pDockQ_unique_residues}
    
    def _calculate_lis_scores(self, unique_chains, chains, pae_matrix):
        """Calculate LIS scores."""
        LIS = self._init_chainpairdict_zeros(unique_chains)

        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue

                mask = (chains[:, None] == chain1) & (chains[None, :] == chain2)
                selected_pae = pae_matrix[mask]

                if selected_pae.size > 0:
                    valid_pae = selected_pae[selected_pae <= 12]
                    if valid_pae.size > 0:
                        scores = (12 - valid_pae) / 12
                        LIS[chain1][chain2] = np.mean(scores)
                    else:
                        LIS[chain1][chain2] = 0.0
                else:
                    LIS[chain1][chain2] = 0.0

        return LIS
    
    def _calculate_ipsae_scores(self, unique_chains, chains, pae_matrix, residues, chain_pair_type, numres, distances):
        """Calculate comprehensive ipSAE scores."""
        # Initialize data structures
        iptm_d0chn_byres = {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
        ipsae_d0chn_byres = {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
        ipsae_d0dom_byres = {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
        ipsae_d0res_byres = {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

        iptm_d0chn_asym = self._init_chainpairdict_zeros(unique_chains)
        ipsae_d0res_asym = self._init_chainpairdict_zeros(unique_chains)
        ipsae_d0chn_asym = self._init_chainpairdict_zeros(unique_chains)
        ipsae_d0dom_asym = self._init_chainpairdict_zeros(unique_chains)
        
        ipsae_d0res_max = self._init_chainpairdict_zeros(unique_chains)
        ipsae_d0chn_max = self._init_chainpairdict_zeros(unique_chains)
        ipsae_d0dom_max = self._init_chainpairdict_zeros(unique_chains)

        unique_residues_chain1 = self._init_chainpairdict_set(unique_chains)
        unique_residues_chain2 = self._init_chainpairdict_set(unique_chains)
        dist_unique_residues_chain1 = self._init_chainpairdict_set(unique_chains)
        dist_unique_residues_chain2 = self._init_chainpairdict_set(unique_chains)

        n0chn = self._init_chainpairdict_zeros(unique_chains)
        n0dom = self._init_chainpairdict_zeros(unique_chains)
        d0chn = self._init_chainpairdict_zeros(unique_chains)
        d0dom = self._init_chainpairdict_zeros(unique_chains)
        n0res = self._init_chainpairdict_zeros(unique_chains)
        d0res = self._init_chainpairdict_zeros(unique_chains)
        
        ptm_func_vec = np.vectorize(self.ptm_func)

        # Calculate scores for each chain pair
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 == chain2:
                    continue

                n0chn[chain1][chain2] = np.sum(chains == chain1) + np.sum(chains == chain2)
                d0chn[chain1][chain2] = self.calc_d0(n0chn[chain1][chain2], chain_pair_type[chain1][chain2])
                ptm_matrix_d0chn = ptm_func_vec(pae_matrix, d0chn[chain1][chain2])

                valid_pairs_iptm = (chains == chain2)
                valid_pairs_matrix = (chains == chain2) & (pae_matrix < self.pae_cutoff)

                for i in range(numres):
                    if chains[i] != chain1:
                        continue

                    valid_pairs_ipsae = valid_pairs_matrix[i]
                    iptm_d0chn_byres[chain1][chain2][i] = ptm_matrix_d0chn[i, valid_pairs_iptm].mean() if valid_pairs_iptm.any() else 0.0
                    ipsae_d0chn_byres[chain1][chain2][i] = ptm_matrix_d0chn[i, valid_pairs_ipsae].mean() if valid_pairs_ipsae.any() else 0.0

                    if valid_pairs_ipsae.any():
                        iresnum = residues[i]['resnum']
                        unique_residues_chain1[chain1][chain2].add(iresnum)
                        for j in np.where(valid_pairs_ipsae)[0]:
                            jresnum = residues[j]['resnum']
                            unique_residues_chain2[chain1][chain2].add(jresnum)
                    
                    # Track unique residues contributing to iptm in interface
                    valid_pairs = (chains == chain2) & (pae_matrix[i] < self.pae_cutoff) & (distances[i] < self.dist_cutoff)
                    if valid_pairs.any():
                        iresnum=residues[i]['resnum']
                        dist_unique_residues_chain1[chain1][chain2].add(iresnum)
                        for j in np.where(valid_pairs)[0]:
                            jresnum=residues[j]['resnum']
                            dist_unique_residues_chain2[chain1][chain2].add(jresnum)


                # Calculate domain-based scores
                residues_1 = len(unique_residues_chain1[chain1][chain2])
                residues_2 = len(unique_residues_chain2[chain1][chain2])
                n0dom[chain1][chain2] = residues_1 + residues_2
                d0dom[chain1][chain2] = self.calc_d0(n0dom[chain1][chain2], chain_pair_type[chain1][chain2])

                ptm_matrix_d0dom = ptm_func_vec(pae_matrix, d0dom[chain1][chain2])
                valid_pairs_matrix = (chains == chain2) & (pae_matrix < self.pae_cutoff)

                n0res_byres_all = np.sum(valid_pairs_matrix, axis=1)
                d0res_byres_all = self.calc_d0_array(n0res_byres_all, chain_pair_type[chain1][chain2])

                for i in range(numres):
                    if chains[i] != chain1:
                        continue
                    valid_pairs = valid_pairs_matrix[i]
                    ipsae_d0dom_byres[chain1][chain2][i] = ptm_matrix_d0dom[i, valid_pairs].mean() if valid_pairs.any() else 0.0

                    ptm_row_d0res = ptm_func_vec(pae_matrix[i], d0res_byres_all[i])
                    ipsae_d0res_byres[chain1][chain2][i] = ptm_row_d0res[valid_pairs].mean() if valid_pairs.any() else 0.0

                # Calculate asymmetric scores (max values for each direction)
                interchain_values = ipsae_d0res_byres[chain1][chain2]
                max_index = np.argmax(interchain_values)
                ipsae_d0res_asym[chain1][chain2] = interchain_values[max_index]
                n0res[chain1][chain2] = n0res_byres_all[max_index]
                d0res[chain1][chain2] = d0res_byres_all[max_index]

                interchain_values = ipsae_d0chn_byres[chain1][chain2]
                max_index = np.argmax(interchain_values)
                ipsae_d0chn_asym[chain1][chain2] = interchain_values[max_index]

                interchain_values = iptm_d0chn_byres[chain1][chain2]
                max_index = np.argmax(interchain_values)
                iptm_d0chn_asym[chain1][chain2] = interchain_values[max_index]

                interchain_values = ipsae_d0dom_byres[chain1][chain2]
                max_index = np.argmax(interchain_values)
                ipsae_d0dom_asym[chain1][chain2] = interchain_values[max_index]

        # Calculate maximum scores for each chain pair (symmetric)
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 >= chain2:
                    continue

                # Calculate max values between A->B and B->A
                ipsae_d0res_max[chain1][chain2] = max(ipsae_d0res_asym[chain1][chain2], ipsae_d0res_asym[chain2][chain1])
                ipsae_d0res_max[chain2][chain1] = ipsae_d0res_max[chain1][chain2]

                ipsae_d0chn_max[chain1][chain2] = max(ipsae_d0chn_asym[chain1][chain2], ipsae_d0chn_asym[chain2][chain1])
                ipsae_d0chn_max[chain2][chain1] = ipsae_d0chn_max[chain1][chain2]

                ipsae_d0dom_max[chain1][chain2] = max(ipsae_d0dom_asym[chain1][chain2], ipsae_d0dom_asym[chain2][chain1])
                ipsae_d0dom_max[chain2][chain1] = ipsae_d0dom_max[chain1][chain2]

        return {
            'ipsae_d0res_asym': ipsae_d0res_asym,
            'ipsae_d0chn_asym': ipsae_d0chn_asym,
            'ipsae_d0dom_asym': ipsae_d0dom_asym,
            'iptm_d0chn_asym': iptm_d0chn_asym,
            'ipsae_d0res_max': ipsae_d0res_max,
            'ipsae_d0chn_max': ipsae_d0chn_max,
            'ipsae_d0dom_max': ipsae_d0dom_max,
            'n0chn': n0chn,
            'n0dom': n0dom,
            'd0chn': d0chn,
            'd0dom': d0dom,
            'n0res': n0res,
            'd0res': d0res,
            'unique_residues_chain1': unique_residues_chain1,
            'unique_residues_chain2': unique_residues_chain2,
            'dist_unique_residues_chain1': dist_unique_residues_chain1,
            'dist_unique_residues_chain2': dist_unique_residues_chain2,
            'ipsae_d0res_byres': ipsae_d0res_byres,
            'ipsae_d0chn_byres': ipsae_d0chn_byres,
            'ipsae_d0dom_byres': ipsae_d0dom_byres,
            'iptm_d0chn_byres': iptm_d0chn_byres,
            'n0res_byres': n0res_byres_all,
            'd0res_byres': d0res_byres_all
        }
    
    def _save_outputs(self, results: Dict, output_dir: Path, stem: str):
        """Save comprehensive output files matching original script format."""
        pae_string = str(int(self.pae_cutoff))
        if self.pae_cutoff < 10:
            pae_string = "0" + pae_string
        dist_string = str(int(self.dist_cutoff))
        if self.dist_cutoff < 10:
            dist_string = "0" + dist_string

        # Main output file
        output_file = output_dir / f"{stem}_{pae_string}_{dist_string}.txt"
        byres_file = output_dir / f"{stem}_{pae_string}_{dist_string}_byres.txt"
        pml_file = output_dir / f"{stem}_{pae_string}_{dist_string}.pml"

        self._write_main_output(output_file, results, stem)
        self._write_byres_output(byres_file, results)
        self._write_pymol_script(pml_file, results)
        
        print(f"Results saved to: {output_file}")
        print(f"By-residue results saved to: {byres_file}")
        print(f"PyMOL script saved to: {pml_file}")
    
    def _write_main_output(self, output_file: Path, results: Dict, stem: str):
        """Write the main output file with chain pair scores."""
        unique_chains = results['unique_chains']
        pdockq_scores = results['pdockq_scores']
        lis_scores = results['lis_scores']
        ipsae_scores = results['ipsae_scores']
        global_scores = results['global_scores']
        file_type = results['file_type']

        # Create chain pairs
        chainpairs = set()
        for chain1 in unique_chains:
            for chain2 in unique_chains:
                if chain1 >= chain2:
                    continue
                chainpairs.add(chain1 + "-" + chain2)

        with open(output_file, 'w') as f:
            f.write("\nChn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")

            for pair in sorted(chainpairs):
                chain_a, chain_b = pair.split("-")
                
                for pair_tuple in [(chain_a, chain_b), (chain_b, chain_a)]:
                    chain1, chain2 = pair_tuple

                    # Get global ipTM value
                    if file_type == "af2":
                        iptm_af = global_scores.get('iptm', -1.0)
                    elif file_type == "af3":
                        chain_pair_iptm = global_scores.get('chain_pair_iptm', [])
                        nchain1 = ord(chain1) - ord('A')
                        nchain2 = ord(chain2) - ord('A')
                        if nchain1 < len(chain_pair_iptm) and nchain2 < len(chain_pair_iptm[nchain1]):
                            iptm_af = chain_pair_iptm[nchain1][nchain2]
                        else:
                            iptm_af = 0.0
                    elif file_type == "boltz1":
                        pair_chains_iptm = global_scores.get('pair_chains_iptm', {})
                        nchain1 = ord(chain1) - ord('A')
                        nchain2 = ord(chain2) - ord('A')
                        iptm_af = pair_chains_iptm.get(str(nchain1), {}).get(str(nchain2), 0.0)
                    else:
                        iptm_af = -1.0

                    residues_1 = len(ipsae_scores['unique_residues_chain1'][chain1][chain2])
                    residues_2 = len(ipsae_scores['unique_residues_chain2'][chain1][chain2])
                    dist_residues_1 = len(ipsae_scores['dist_unique_residues_chain1'][chain1][chain2])
                    dist_residues_2 = len(ipsae_scores['dist_unique_residues_chain2'][chain1][chain2])

                    # Asymmetric scores
                    outstring = f'{chain1}    {chain2}     {int(self.pae_cutoff):3d}  {int(self.dist_cutoff):3d}  {"asym":5} ' + (
                        f'{ipsae_scores["ipsae_d0res_asym"][chain1][chain2]:8.6f}    '
                        f'{ipsae_scores["ipsae_d0chn_asym"][chain1][chain2]:8.6f}    '
                        f'{ipsae_scores["ipsae_d0dom_asym"][chain1][chain2]:8.6f}    '
                        f'{iptm_af:5.3f}    '
                        f'{ipsae_scores["iptm_d0chn_asym"][chain1][chain2]:8.6f}    '
                        f'{pdockq_scores["pDockQ"][chain1][chain2]:8.4f}   '
                        f'{pdockq_scores["pDockQ2"][chain1][chain2]:8.4f}   '
                        f'{lis_scores[chain1][chain2]:8.4f}   '
                        f'{int(ipsae_scores["n0res"][chain1][chain2]):5d}  '
                        f'{int(ipsae_scores["n0chn"][chain1][chain2]):5d}  '
                        f'{int(ipsae_scores["n0dom"][chain1][chain2]):5d}  '
                        f'{ipsae_scores["d0res"][chain1][chain2]:6.2f}  '
                        f'{ipsae_scores["d0chn"][chain1][chain2]:6.2f}  '
                        f'{ipsae_scores["d0dom"][chain1][chain2]:6.2f}  '
                        f'{residues_1:5d}   '
                        f'{residues_2:5d}   '
                        f'{dist_residues_1:5d}   '
                        f'{dist_residues_2:5d}   '
                        f'{stem}\n')
                    f.write(outstring)

                    # Max scores (only for chain1 > chain2 to avoid duplication)
                    if chain1 > chain2:
                        max_iptm = iptm_af 
                        if file_type == "boltz1":
                           max_iptm=max(global_scores.get('pair_chains_iptm', {}).get(str(ord(chain1)-ord('A')), {}).get(str(ord(chain2)-ord('A')), 0.0), global_scores.get('pair_chains_iptm', {}).get(str(ord(chain2)-ord('A')), {}).get(str(ord(chain1)-ord('A')), 0.0))

                        max_pdockq2 = max(pdockq_scores["pDockQ2"][chain1][chain2], pdockq_scores["pDockQ2"][chain2][chain1])
                        max_lis = (lis_scores[chain1][chain2] + lis_scores[chain2][chain1]) / 2.0
                        residues_1_max = max(len(ipsae_scores['unique_residues_chain2'][chain1][chain2]), len(ipsae_scores['unique_residues_chain1'][chain2][chain1]))
                        residues_2_max = max(len(ipsae_scores['unique_residues_chain1'][chain1][chain2]), len(ipsae_scores['unique_residues_chain2'][chain2][chain1]))
                        dist_residues_1_max = max(len(ipsae_scores['dist_unique_residues_chain2'][chain1][chain2]), len(ipsae_scores['dist_unique_residues_chain1'][chain2][chain1]))
                        dist_residues_2_max = max(len(ipsae_scores['dist_unique_residues_chain1'][chain1][chain2]), len(ipsae_scores['dist_unique_residues_chain2'][chain2][chain1]))


                        outstring = f'{chain2}    {chain1}     {int(self.pae_cutoff):3d}  {int(self.dist_cutoff):3d}  {"max":5} ' + (
                            f'{ipsae_scores["ipsae_d0res_max"][chain1][chain2]:8.6f}    '
                            f'{ipsae_scores["ipsae_d0chn_max"][chain1][chain2]:8.6f}    '
                            f'{ipsae_scores["ipsae_d0dom_max"][chain1][chain2]:8.6f}    '
                            f'{max_iptm:5.3f}    '
                            f'{max(ipsae_scores["iptm_d0chn_asym"][chain1][chain2], ipsae_scores["iptm_d0chn_asym"][chain2][chain1]):8.6f}    '
                            f'{pdockq_scores["pDockQ"][chain1][chain2]:8.4f}   '
                            f'{max_pdockq2:8.4f}   '
                            f'{max_lis:8.4f}   '
                            f'{int(max(ipsae_scores["n0res"][chain1][chain2], ipsae_scores["n0res"][chain2][chain1])):5d}  '
                            f'{int(ipsae_scores["n0chn"][chain1][chain2]):5d}  '
                            f'{int(max(ipsae_scores["n0dom"][chain1][chain2], ipsae_scores["n0dom"][chain2][chain1])):5d}  '
                            f'{max(ipsae_scores["d0res"][chain1][chain2], ipsae_scores["d0res"][chain2][chain1]):6.2f}  '
                            f'{ipsae_scores["d0chn"][chain1][chain2]:6.2f}  '
                            f'{max(ipsae_scores["d0dom"][chain1][chain2], ipsae_scores["d0dom"][chain2][chain1]):6.2f}  '
                            f'{residues_1_max:5d}   '
                            f'{residues_2_max:5d}   '
                            f'{dist_residues_1_max:5d}   '
                            f'{dist_residues_2_max:5d}   '
                            f'{stem}\n')
                        f.write(outstring)

                f.write("\n")
    
    def _write_byres_output(self, byres_file: Path, results: Dict):
        """Write by-residue output file."""
        residues = results['residues']
        chains = results['chains']
        plddt = results['plddt']
        unique_chains = results['unique_chains']
        ipsae_scores = results['ipsae_scores']

        with open(byres_file, 'w') as f:
            f.write("i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT      n0chn  n0dom  n0res    d0chn     d0dom     d0res   ipTM_pae  ipSAE_d0chn ipSAE_d0dom    ipSAE \n")
            
            for i, residue in enumerate(residues):
                for chain1 in unique_chains:
                    for chain2 in unique_chains:
                        if chain1 == chain2:
                            continue
                        if chains[i] != chain1:
                            continue

                        outstring = f'{i+1:<4d}    ' + (
                            f'{chain1:4}      '
                            f'{chain2:4}      '
                            f'{residue["resnum"]:4d}           '
                            f'{residue["res"]:3}        '
                            f'{plddt[i]:8.2f}         '
                            f'{int(ipsae_scores["n0chn"][chain1][chain2]):5d}  '
                            f'{int(ipsae_scores["n0dom"][chain1][chain2]):5d}  '
                            f'{int(ipsae_scores["n0res_byres"][i]):5d}  '
                            f'{ipsae_scores["d0chn"][chain1][chain2]:8.3f}  '
                            f'{ipsae_scores["d0dom"][chain1][chain2]:8.3f}  '
                            f'{ipsae_scores["d0res_byres"][i]:8.3f}   '
                            f'{ipsae_scores["iptm_d0chn_byres"][chain1][chain2][i]:8.4f}    '
                            f'{ipsae_scores["ipsae_d0chn_byres"][chain1][chain2][i]:8.4f}    '
                            f'{ipsae_scores["ipsae_d0dom_byres"][chain1][chain2][i]:8.4f}    '
                            f'{ipsae_scores["ipsae_d0res_byres"][chain1][chain2][i]:8.4f}\n'
                        )
                        f.write(outstring)
    
    def _write_pymol_script(self, pml_file: Path, results: Dict):
        """Write PyMOL visualization script."""
        unique_chains = results['unique_chains']
        ipsae_scores = results['ipsae_scores']
        
        chaincolor = {
            'A': 'magenta', 'B': 'marine', 'C': 'lime', 'D': 'orange',
            'E': 'yellow', 'F': 'cyan', 'G': 'lightorange', 'H': 'pink',
            'I': 'deepteal', 'J': 'forest', 'K': 'lightblue', 'L': 'slate',
            'M': 'violet', 'N': 'arsenic', 'O': 'iodine', 'P': 'silver',
            'Q': 'red', 'R': 'sulfur', 'S': 'purple', 'T': 'olive',
            'U': 'palegreen', 'V': 'green', 'W': 'blue', 'X': 'palecyan',
            'Y': 'limon', 'Z': 'chocolate'
        }

        with open(pml_file, 'w') as f:
            f.write("# PyMOL script for ipSAE visualization\n")
            f.write("# Generated by ipSAE package\n\n")
            
            for chain1 in unique_chains:
                for chain2 in unique_chains:
                    if chain1 >= chain2:
                        continue
                    
                    color1 = chaincolor.get(chain1, 'magenta')
                    color2 = chaincolor.get(chain2, 'marine')
                    
                    chain1_residues = ipsae_scores['unique_residues_chain1'][chain1][chain2]
                    chain2_residues = ipsae_scores['unique_residues_chain2'][chain1][chain2]
                    
                    if chain1_residues and chain2_residues:
                        chain1_range = self._contiguous_ranges(chain1_residues)
                        chain2_range = self._contiguous_ranges(chain2_residues)
                        
                        chain_pair = f'color_{chain1}_{chain2}'
                        chain1_selection = f'chain  {chain1} and resi {chain1_range}'
                        chain2_selection = f'chain  {chain2} and resi {chain2_range}'
                        
                        f.write(f'alias {chain_pair}, color gray80, all; color {color1}, {chain1_selection}; color {color2}, {chain2_selection}\n\n')
    
    def _contiguous_ranges(self, numbers):
        """Convert set of numbers to contiguous ranges for PyMOL."""
        if not numbers:
            return ""
        
        sorted_numbers = sorted(numbers)
        ranges = []
        start = sorted_numbers[0]
        end = start

        for number in sorted_numbers[1:]:
            if number == end + 1:
                end = number
            else:
                if start == end:
                    ranges.append(f"{start}")
                else:
                    ranges.append(f"{start}-{end}")
                start = end = number

        # Append the last range
        if start == end:
            ranges.append(f"{start}")
        else:
            ranges.append(f"{start}-{end}")

        return '+'.join(ranges)
    
    @staticmethod
    def ptm_func(x, d0):
        """PTM function calculation."""
        return 1.0 / (1 + (x / d0) ** 2.0)
    
    @staticmethod
    def calc_d0(L, pair_type='protein'):
        """Calculate d0 value."""
        L = float(L)
        if L < 27:
            L = 27
        min_value = 1.0
        if pair_type == 'nucleic_acid':
            min_value = 2.0
        d0 = 1.24 * (L - 15) ** (1.0/3.0) - 1.8
        return max(min_value, d0)
    
    @staticmethod
    def calc_d0_array(L, pair_type='protein'):
        """Calculate d0 values for arrays."""
        L = np.array(L, dtype=float)
        L = np.maximum(27, L)
        min_value = 1.0
        if pair_type == 'nucleic_acid':
            min_value = 2.0
        return np.maximum(min_value, 1.24 * (L - 15) ** (1.0/3.0) - 1.8)


def calculate_ipsae(pae_file: Union[str, Path], structure_file: Union[str, Path], 
                   pae_cutoff: float = 10.0, dist_cutoff: float = 10.0,
                   output_dir: Optional[Union[str, Path]] = None) -> Dict:
    """
    Convenience function to calculate ipSAE scores.
    
    Args:
        pae_file: Path to PAE file
        structure_file: Path to structure file
        pae_cutoff: PAE cutoff value
        dist_cutoff: Distance cutoff value
        output_dir: Output directory
        
    Returns:
        Dictionary containing calculated scores
    """
    calculator = IpsaeCalculator(pae_cutoff, dist_cutoff)
    return calculator.calculate(pae_file, structure_file, output_dir)
