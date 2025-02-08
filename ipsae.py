# Usage:

#  python ipsae.py <path_to_json_file> <path_to_af2_pdb_file>  <pae_cutoff> <dist_cutoff>
#  python ipsae.py <path_to_json_file> <path_to_af3_cif_file>  <pae_cutoff> <dist_cutoff>

import sys
import numpy as np
import json
np.set_printoptions(threshold=np.inf)  # for printing out full numpy arrays for debugging

# Input and output files and parameters
json_file_path =   sys.argv[1]
pdb_path =         sys.argv[2]
pae_cutoff =       float(sys.argv[3])
dist_cutoff =      float(sys.argv[4])
pae_string =       str(int(pae_cutoff))
if pae_cutoff<10:  pae_string="0"+pae_string
dist_string =      str(int(dist_cutoff))
if dist_cutoff<10: dist_string="0"+dist_string

if ".pdb" in pdb_path:
    path_stem =     f'{pdb_path.replace(".pdb","")}_{pae_string}_{dist_string}'
    af2 = True
    af3 = False
elif ".cif" in pdb_path:
    path_stem =     f'{pdb_path.replace(".cif","")}_{pae_string}_{dist_string}'
    af2 = False
    af3 = True
else:
    print("Wrong PDB file type ", pdb_path)
    exit()
    
file_path =        path_stem + ".txt"
file2_path =       path_stem + "_byres.txt"
pml_path =         path_stem + ".pml"
OUT =              open(file_path,'w')
PML =              open(pml_path,'w')
OUT2 =             open(file2_path,'w')

# Ensure correct usage
if len(sys.argv) < 5:
    print("Usage for AF2:")
    print("   python ipsae.py <path_to_json_file> <path_to_pdb_file> <pae_cutoff> <dist_cutoff>")
    print("   Example: python ipsae.py RAF1_KSR1_scores_rank_001_alphafold2_multimer_v3_model_4_seed_003.json RAF1_KSR1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_003.pdb 10 10")
    print("")
    print("Usage for AF3:")
    print("python ipsae.py <path_to_json_file> <path_to_pdb_file> <pae_cutoff> <dist_cutoff>")
    print("   Example: python ipsae.py fold_aurka_tpx2_full_data_0.json  fold_aurka_tpx2_model_0.cif 10 10")
    sys.exit(1)

# Define the ptm function
def ptm_func(x,d0):
    return 1.0/(1+(x/d0)**2.0)  
ptm_func_vec=np.vectorize(ptm_func)  # vector version

# Define the d0 functions for numbers and arrays; minimum value = 1.0; from Yang and Skolnick, PROTEINS: Structure, Function, and Bioinformatics 57:702–710 (2004)
def calc_d0(L):
    L=float(L)
    if L<27: return 1.0
    return 1.24*(L-15)**(1.0/3.0) - 1.8

def calc_d0_array(L):
    # Convert L to a NumPy array if it isn't already one (enables flexibility in input types)
    L = np.array(L, dtype=float)
    # Ensure all values of L are at least 19.0
    L = np.maximum(L, 26.523)
    # Calculate d0 using the vectorized operation
    return 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8

# Define the parse_atom_line function for PDB lines (by column)
# parsed_line = parse_atom_line(line)
# line = "ATOM    123  CA  ALA A  15     11.111  22.222  33.333  1.00 20.00           C"
def parse_pdb_atom_line(line):
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

def parse_cif_atom_line(line):
    # ligands do not have residue numbers but modified residues do. Return "None" for ligand.
    # for parsing AF3 mmCIF files
    #ATOM   1294 N  N     . ARG A 1 159 ? 5.141   -14.096 10.526  1.00 95.62 159 A 1
    #ATOM   1295 C  CA    . ARG A 1 159 ? 4.186   -13.376 11.366  1.00 96.27 159 A 1
    #ATOM   1296 C  C     . ARG A 1 159 ? 2.976   -14.235 11.697  1.00 96.42 159 A 1
    #ATOM   1297 O  O     . ARG A 1 159 ? 2.654   -15.174 10.969  1.00 95.46 159 A 1
    #ATOM   1298 C  CB    . ARG A 1 159 ? 3.798   -12.052 10.695  1.00 95.57 159 A 1
    #HETATM 1305 N  N     . TPO A 1 160 ? 2.328   -13.853 12.742  1.00 96.42 160 A 1
    #HETATM 1306 C  CA    . TPO A 1 160 ? 1.081   -14.560 13.218  1.00 96.78 160 A 1
    #HETATM 1307 C  C     . TPO A 1 160 ? -2.115  -11.668 12.263  1.00 96.19 160 A 1
    #HETATM 1308 O  O     . TPO A 1 160 ? -1.790  -11.556 11.113  1.00 95.75 160 A 1
    # ...
    #HETATM 2639 MG MG    . MG  D 4 .   ? -7.262  2.709   4.825   1.00 91.47 1   D 1 
    #HETATM 2640 MG MG    . MG  E 5 .   ? -4.994  2.251   8.755   1.00 85.96 1   E 1 
    #HETATM 2608 P  PG    . ATP C 3 .   ? -6.858  4.182   10.275  1.00 84.94 1   C 1 
    #HETATM 2609 O  O1G   . ATP C 3 .   ? -6.178  5.238   11.074  1.00 75.56 1   C 1 
    #HETATM 2610 O  O2G   . ATP C 3 .   ? -5.889  3.166   9.748   1.00 75.15 1   C 1 

    linelist =        line.split()
    atom_num =        linelist[1]
    atom_name =       linelist[3]
    residue_name =    linelist[5]
    chain_id =        linelist[6]
    residue_seq_num = linelist[8]
    x =               linelist[10]
    y =               linelist[11]
    z =               linelist[12]

    if residue_seq_num == ".": return None   # ligand atom

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

# function for pymol scripts
def contiguous_ranges(numbers):
    if not numbers:  # Check if the set is empty
        return
    
    sorted_numbers = sorted(numbers)  # Sort the numbers
    start = sorted_numbers[0]
    end = start
    ranges = []  # List to store ranges

    def format_range(start, end):
        if start == end:
            return f"{start}"
        else:
            return f"{start}-{end}"

    for number in sorted_numbers[1:]:
        if number == end + 1:
            end = number
        else:
            ranges.append(format_range(start, end))
            start = end = number
    
    # Append the last range after the loop
    ranges.append(format_range(start, end))

    # Join all ranges with a plus sign and print the result
    string='+'.join(ranges)
    return(string)


# Load residues from AlphaFold PDB or mmCIF file into lists; each residue is a dictionary
# Save CA coordinates
# Read PDB file to get distances, chainids, and residue numbers
residues = []
chains = []
with open(pdb_path, 'r') as PDB:
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if af2:
                atom=parse_pdb_atom_line(line)
            else:
                atom=parse_cif_atom_line(line)
            if atom is None: continue
            if atom['atom_name'] == "CA":
                residues.append({
                    'atom_num': atom['atom_num'],
                    'coor': np.array([atom['x'], atom['y'], atom['z']]),
                    'res': atom['residue_name'],
                    'chainid': atom['chain_id'],
                    'resnum': atom['residue_seq_num'],
                    'residue': f"{atom['residue_name']:3}   {atom['chain_id']:3} {atom['residue_seq_num']:4}"
                    
                })
                chains.append(atom['chain_id'])

# Convert to np arrays
numres = len(residues)
CA_atom_num=  np.array([res['atom_num']-1 for res in residues])  # for AF3 atom indexing from 0
coordinates = np.array([res['coor']       for res in residues])
chains = np.array(chains)
unique_chains = np.unique(chains)

# Calculate distance matrix using NumPy broadcasting
distances = np.sqrt(((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :])**2).sum(axis=2))

# Load json data and extract plddt and pae_matrix (and ptm_matrix if available)
with open(json_file_path, 'r') as file:
    data = json.load(file)

#for item in data:
#    print("json",item, data[item])
#    print("\n")

#fold_aurka_0_tpx2_0_full_data_0.json
#fold_aurka_0_tpx2_0_summary_confidences_0.json
#fold_aurka_0_tpx2_0_model_0.cif


ptm_matrix_af2=np.zeros((numres,numres))
if af2:
    iptm_af2 = float(data['iptm'])
    ptm_af2  = float(data['ptm'])
    plddt =      np.array(data['plddt'])
    pae_matrix = np.array(data['pae'])
    if 'ptm_matrix' in data:
        ptm_matrix_af2 = np.array(data['ptm_matrix'])
        have_ptm_matrix_af2=1
    else:
        ptm_matrix_af2=np.zeros((numres,numres))
        have_ptm_matrix_af2=0

if af3 and ("full_data" in json_file_path):

    atom_plddts=np.array(data['atom_plddts'])
    plddt=atom_plddts[CA_atom_num]  # pull out residue plddts from Calpha atoms
    pae_matrix_af3 = np.array(data['pae'])

    # Get pairwise residue PAE matrix by identifying one token per protein residue.
    # Modified residues have separate tokens for each atom, so need to pull out Calpha atom as token
    # Skip ligands
    
    # Create a composite key for residue and chain
    token_chain_ids = np.array(data['token_chain_ids'])
    token_res_ids =   np.array(data['token_res_ids'])
    composite_key =   np.core.defchararray.add(token_res_ids.astype(str), token_chain_ids)
    #    print("composite_key",composite_key)
    
    # Detect changes including the first element
    changes = np.concatenate(([True], composite_key[1:] != composite_key[:-1]))

    # Initialize the mask array as zeroes
    mask = np.zeros_like(token_res_ids)

    # Set mask to 1 at each new composite key
    mask[changes] = 1

    # Loop to handle the second occurrence within the same chain
    #    i=0; print("mask",i,composite_key[i],mask[i])
    for i in range(1, len(composite_key)):
        if token_chain_ids[i] not in unique_chains:  # set ligands to 0
            mask[i]=0
            continue
        if composite_key[i] == composite_key[i - 1]:
            mask[i] = 0  # Set all repeats to 0
            if composite_key[i] == composite_key[i - 1] and composite_key[i] != composite_key[i - 2]:
                mask[i] = 1  # Set the second occurrence to 1 and the first occurrence to 0
                mask[i-1] = 0
                #                print("mask2",i,composite_key[i-1],composite_key[i-1], composite_key[i],mask[i-2],mask[i-1],mask[i])


    for i in range(0,len(composite_key)):
        print("mask",i,composite_key[i],mask[i])

    # Set pae_matrix for AF3 from subset of full PAE matrix from json file
    pae_matrix = pae_matrix_af3[np.ix_(mask.astype(bool), mask.astype(bool))]

    print("Filtered PAE Matrix:")
    print(pae_matrix)        

    # Get iptm matrix from AF3 summary_confidences file
    json_summary_file_path=json_file_path.replace("full_data","summary_confidences")
    with open(json_summary_file_path,'r') as file:
        data_summary=json.load(file)
    #    for item in data_summary:
    #        print("json_sum",item,data_summary[item])
    #        print("\n")
        
    af3_chain_pair_iptm=   {chain1: {chain2: 0     for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
    af3_chain_pair_iptm_data=data_summary['chain_pair_iptm']
    for chain1 in unique_chains:
        nchain1=  ord(chain1) - ord('A')
        for chain2 in unique_chains:
            if chain1 == chain2: continue
            nchain2=ord(chain2) - ord('A')
            af3_chain_pair_iptm[chain1][chain2]=af3_chain_pair_iptm_data[nchain1][nchain2]
            print("iptm",chain1,chain2,af3_chain_pair_iptm[chain1][chain2])
    
# Compute chain-pair-specific interchain PTM and PAE, count valid pairs, and count unique residues
# First, create dictionaries of appropriate size: top keys are chain1 and chain2 where chain1 != chain2
# Nomenclature:
# iptm_af2_d0num = calculate iptm from ptm matrix from af2 itself (if present in json file); d0num = numres in full AF2 (multichain) complex
# iptm_pae_d0chn = calculate iptm from ptm matrix calculated from PAEs (from json file); d0chn = numres in chain pair = len(chain1) + len(chain2)
# ipsa_af2_d0num = calculate ipsae from ptm matrix from af2 itself (if available)
# ipsa_pae_d0chn = calculate ipsae from ptm matrix from PAEs
# ipsa_pae_d0dom = calculate ipsae from ptm matrix from PAEs; d0 from number of residues in chain1 and chain2 that have interchain PAE<cutoff
# ipsa_pae_d0res = calculate ipsae from ptm matrix from PAEs; d0 from number of residues in chain2 that have interchain PAE<cutoff given residue in chain1
# 
# for each chain_pair iptm/ipsae, there is (for example)
# ipsa_pae_d0res_byres = by-residue array;
# ipsa_pae_d0res_asym  = asymmetric pair value (A->B is different from B->A)
# ipsa_pae_d0res_max   = maximum of A->B and B->A value
# ipsa_pae_d0res_asymres = identify of residue that provides each asym maximum
# ipsa_pae_d0res_maxres =  identify of residue that provides each maximum over both chains
#
# n0num = number of residues in whole complex provided by AF2 model
# n0chn = number of residues in chain pair = len(chain1) + len(chain2)
# n0dom = number of residues in chain pair that have good PAE values (<cutoff)
# n0res = number of residues in chain2 that have good PAE residues for each residue of chain1
iptm_af2_d0num_byres =        {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
iptm_pae_d0chn_byres =        {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_af2_d0num_byres =        {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0chn_byres =        {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0dom_byres =        {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0res_byres =        {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

iptm_af2_d0num_asym  =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
iptm_pae_d0chn_asym  =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_af2_d0num_asym  =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0chn_asym  =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0dom_asym  =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0res_asym  =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

iptm_af2_d0num_max   =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
iptm_pae_d0chn_max   =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_af2_d0num_max   =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0chn_max   =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0dom_max   =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0res_max   =        {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

iptm_af2_d0num_asymres  =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
iptm_pae_d0chn_asymres  =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_af2_d0num_asymres  =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0chn_asymres  =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0dom_asymres  =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0res_asymres  =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

iptm_af2_d0num_maxres   =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
iptm_pae_d0chn_maxres   =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_af2_d0num_maxres   =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0chn_maxres   =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0dom_maxres   =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
ipsa_pae_d0res_maxres   =     {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

n0num                 =       numres
n0chn                 =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
n0dom                 =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
n0dom_max             =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
n0res                 =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
n0res_max             =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
n0res_byres           =       {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

d0num                 =       calc_d0(n0num)
d0chn                 =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
d0dom                 =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
d0dom_max             =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
d0res                 =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
d0res_max             =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
d0res_byres           =       {chain1: {chain2: np.zeros(numres) for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}

valid_pair_counts      =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
dist_valid_pair_counts =       {chain1: {chain2: 0                for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
unique_residues_chain1 =      {chain1: {chain2: set()            for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
unique_residues_chain2 =      {chain1: {chain2: set()            for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
dist_unique_residues_chain1 = {chain1: {chain2: set()            for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}
dist_unique_residues_chain2 = {chain1: {chain2: set()            for chain2 in unique_chains if chain1 != chain2} for chain1 in unique_chains}


# calculate regular ipTM with and without PAE cutoff and ptm_matrix calculated by AF2

for chain1 in unique_chains:
    for chain2 in unique_chains:
        if chain1 == chain2:
            continue

        n0chn[chain1][chain2]=np.sum( chains==chain1) + np.sum(chains==chain2) # total number of residues in chain1 and chain2
        d0chn[chain1][chain2]=calc_d0(n0chn[chain1][chain2])
        ptm_matrix_d0chn=np.zeros((numres,numres))
        ptm_matrix_d0chn=ptm_func_vec(pae_matrix,d0chn[chain1][chain2])

        valid_pairs_iptm = (chains == chain2) 
        chain2_mask = (chains == chain2)
        valid_pairs_matrix = (pae_matrix < pae_cutoff) & chain2_mask

        for i in range(numres):
            if chains[i] != chain1:
                continue

            valid_pairs_ipsa = valid_pairs_matrix[i]  # row for residue i of chain1
            iptm_af2_d0num_byres[chain1][chain2][i] =   ptm_matrix_af2[i, valid_pairs_iptm].mean() if valid_pairs_iptm.any() else 0.0
            iptm_pae_d0chn_byres[chain1][chain2][i] = ptm_matrix_d0chn[i, valid_pairs_iptm].mean() if valid_pairs_iptm.any() else 0.0
            ipsa_af2_d0num_byres[chain1][chain2][i] =   ptm_matrix_af2[i, valid_pairs_ipsa].mean() if valid_pairs_ipsa.any() else 0.0
            ipsa_pae_d0chn_byres[chain1][chain2][i] = ptm_matrix_d0chn[i, valid_pairs_ipsa].mean() if valid_pairs_ipsa.any() else 0.0

            # Track unique residues contributing to the IPSA for chain1,chain2
            valid_pair_counts[chain1][chain2] += np.sum(valid_pairs_ipsa)
            if valid_pairs_ipsa.any():
                iresnum=residues[i]['resnum']
                unique_residues_chain1[chain1][chain2].add(iresnum)
                for j in np.where(valid_pairs_ipsa)[0]:
                    jresnum=residues[j]['resnum']
                    unique_residues_chain2[chain1][chain2].add(jresnum)
                    
            # Track unique residues contributing to iptm in interface
            valid_pairs = (chains == chain2) & (pae_matrix[i] < pae_cutoff) & (distances[i] < dist_cutoff)
            dist_valid_pair_counts[chain1][chain2] += np.sum(valid_pairs)

            # Track unique residues contributing to the IPTM
            if valid_pairs.any():
                iresnum=residues[i]['resnum']
                dist_unique_residues_chain1[chain1][chain2].add(iresnum)
                for j in np.where(valid_pairs)[0]:
                    jresnum=residues[j]['resnum']
                    dist_unique_residues_chain2[chain1][chain2].add(jresnum)

# for comparison of ptm matrix from af2 and ptm matrix from PAE values
#for i in range(numres):
#    for j in range(numres):
#        if chains[i]>=chains[j]: continue
#        print(i+1,j+1,chains[i],chains[j],ptm_matrix_af2[i,j],ptm_matrix_d0chn[i,j])
        
# calculate for ptm from pae and for d0=numres and d0=total(pae<cutoff)
OUT2.write("chn1 chn2 i numres plddt  n0chn n0dom   n0res   d0num      d0chn      d0dom   d0res tm_af2_d0num tm_pae_d0chn af2_d0num pae_d0chn pae_d0dom pae_d0res\n")
for chain1 in unique_chains:
    for chain2 in unique_chains:
        if chain1 == chain2:
            continue
        residues_1 = len(unique_residues_chain1[chain1][chain2])
        residues_2 = len(unique_residues_chain2[chain1][chain2])
        n0dom[chain1][chain2]=residues_1+residues_2
        d0dom[chain1][chain2]=calc_d0(n0dom[chain1][chain2])

        ptm_matrix_d0dom=np.zeros((numres,numres))
        ptm_matrix_d0dom=ptm_func_vec(pae_matrix,d0dom[chain1][chain2])

        chain2_mask = (chains == chain2)
        valid_pairs_matrix = (pae_matrix < pae_cutoff) & chain2_mask
        # Assuming valid_pairs_matrix is already defined
        n0res_byres_all = np.sum(valid_pairs_matrix, axis=1)
        d0res_byres_all = calc_d0_array(n0res_byres_all)

        n0res_byres[chain1][chain2] = n0res_byres_all
        d0res_byres[chain1][chain2] = d0res_byres_all

        
        for i in range(numres):
            if chains[i] != chain1:
                continue
            valid_pairs = valid_pairs_matrix[i]
            ipsa_pae_d0dom_byres[chain1][chain2][i] = ptm_matrix_d0dom[i, valid_pairs].mean() if valid_pairs.any() else 0.0

            #            if n0res_byres[chain1][chain2][i]==0: continue

            ptm_row_d0res=np.zeros((numres))
            ptm_row_d0res=ptm_func_vec(pae_matrix[i], d0res_byres[chain1][chain2][i])
            ipsa_pae_d0res_byres[chain1][chain2][i] = ptm_row_d0res[valid_pairs].mean() if valid_pairs.any() else 0.0
            outstring = f'{chain1}  {chain2} '+ (
                f'{i:5d}  '
                f'{int(numres):5d}  '
                f'{plddt[i]:8.3f}   '
                f'{int(n0chn[chain1][chain2]):5d}  '
                f'{int(n0dom[chain1][chain2]):5d}  '
                f'{int(n0res_byres[chain1][chain2][i]):5d}  '
                f'{d0num:8.3f}  '
                f'{d0chn[chain1][chain2]:8.3f}  '
                f'{d0dom[chain1][chain2]:8.3f}  '
                f'{d0res_byres[chain1][chain2][i]:8.3f}  '
                f'{iptm_af2_d0num_byres[chain1][chain2][i]:8.4f}  '
                f'{iptm_pae_d0chn_byres[chain1][chain2][i]:8.4f}  '
                f'{ipsa_af2_d0num_byres[chain1][chain2][i]:8.4f}  '
                f'{ipsa_pae_d0chn_byres[chain1][chain2][i]:8.4f}  '
                f'{ipsa_pae_d0dom_byres[chain1][chain2][i]:8.4f}  '
                f'{ipsa_pae_d0res_byres[chain1][chain2][i]:8.4f}\n'
            )
            OUT2.write(outstring)
            #            print(outstring)
            
# Compute interchain PTM for each chain pair and regular PTM from AF2 and d0=numres
for chain1 in unique_chains:
    for chain2 in unique_chains:
        if chain1 == chain2:
            continue

        interchain_values = iptm_af2_d0num_byres[chain1][chain2]
        max_index = np.argmax(interchain_values)
        iptm_af2_d0num_asym[chain1][chain2] = interchain_values[max_index]
        iptm_af2_d0num_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

        interchain_values = iptm_pae_d0chn_byres[chain1][chain2]
        max_index = np.argmax(interchain_values)
        iptm_pae_d0chn_asym[chain1][chain2] = interchain_values[max_index]
        iptm_pae_d0chn_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

        interchain_values = ipsa_af2_d0num_byres[chain1][chain2]
        max_index = np.argmax(interchain_values)
        ipsa_af2_d0num_asym[chain1][chain2] = interchain_values[max_index]
        ipsa_af2_d0num_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

        interchain_values = ipsa_pae_d0chn_byres[chain1][chain2]
        max_index = np.argmax(interchain_values)
        ipsa_pae_d0chn_asym[chain1][chain2] = interchain_values[max_index]
        ipsa_pae_d0chn_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

        interchain_values = ipsa_pae_d0dom_byres[chain1][chain2]
        max_index = np.argmax(interchain_values)
        ipsa_pae_d0dom_asym[chain1][chain2] = interchain_values[max_index]
        ipsa_pae_d0dom_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"

        interchain_values = ipsa_pae_d0res_byres[chain1][chain2]
        max_index = np.argmax(interchain_values)
        ipsa_pae_d0res_asym[chain1][chain2] = interchain_values[max_index]
        ipsa_pae_d0res_asymres[chain1][chain2] = residues[max_index]['residue'] if max_index is not None else "None"
        n0res[chain1][chain2]=n0res_byres[chain1][chain2][max_index]
        d0res[chain1][chain2]=d0res_byres[chain1][chain2][max_index]

        if chain1 > chain2:
            if iptm_af2_d0num_asym[chain1][chain2] >= iptm_af2_d0num_asym[chain2][chain1]:
                iptm_af2_d0num_max[chain1][chain2]=iptm_af2_d0num_asym[chain1][chain2]
                iptm_af2_d0num_maxres[chain1][chain2]=iptm_af2_d0num_asymres[chain1][chain2]
                iptm_af2_d0num_max[chain2][chain1]=iptm_af2_d0num_asym[chain1][chain2]
                iptm_af2_d0num_maxres[chain2][chain1]=iptm_af2_d0num_asymres[chain1][chain2]
            else:
                iptm_af2_d0num_max[chain1][chain2]=iptm_af2_d0num_asym[chain2][chain1]
                iptm_af2_d0num_maxres[chain1][chain2]=iptm_af2_d0num_asymres[chain2][chain1]
                iptm_af2_d0num_max[chain2][chain1]=iptm_af2_d0num_asym[chain2][chain1]
                iptm_af2_d0num_maxres[chain2][chain1]=iptm_af2_d0num_asymres[chain2][chain1]
                
            if iptm_pae_d0chn_asym[chain1][chain2] >= iptm_pae_d0chn_asym[chain2][chain1]:
                iptm_pae_d0chn_max[chain1][chain2]=iptm_pae_d0chn_asym[chain1][chain2]
                iptm_pae_d0chn_maxres[chain1][chain2]=iptm_pae_d0chn_asymres[chain1][chain2]
                iptm_pae_d0chn_max[chain2][chain1]=iptm_pae_d0chn_asym[chain1][chain2]
                iptm_pae_d0chn_maxres[chain2][chain1]=iptm_pae_d0chn_asymres[chain1][chain2]
            else:
                iptm_pae_d0chn_max[chain1][chain2]=iptm_pae_d0chn_asym[chain2][chain1]
                iptm_pae_d0chn_maxres[chain1][chain2]=iptm_pae_d0chn_asymres[chain2][chain1]
                iptm_pae_d0chn_max[chain2][chain1]=iptm_pae_d0chn_asym[chain2][chain1]
                iptm_pae_d0chn_maxres[chain2][chain1]=iptm_pae_d0chn_asymres[chain2][chain1]
                
            if ipsa_af2_d0num_asym[chain1][chain2] >= ipsa_af2_d0num_asym[chain2][chain1]:
                ipsa_af2_d0num_max[chain1][chain2]=ipsa_af2_d0num_asym[chain1][chain2]
                ipsa_af2_d0num_maxres[chain1][chain2]=ipsa_af2_d0num_asymres[chain1][chain2]
                ipsa_af2_d0num_max[chain2][chain1]=ipsa_af2_d0num_asym[chain1][chain2]
                ipsa_af2_d0num_maxres[chain2][chain1]=ipsa_af2_d0num_asymres[chain1][chain2]
            else:
                ipsa_af2_d0num_max[chain1][chain2]=ipsa_af2_d0num_asym[chain2][chain1]
                ipsa_af2_d0num_maxres[chain1][chain2]=ipsa_af2_d0num_asymres[chain2][chain1]
                ipsa_af2_d0num_max[chain2][chain1]=ipsa_af2_d0num_asym[chain2][chain1]
                ipsa_af2_d0num_maxres[chain2][chain1]=ipsa_af2_d0num_asymres[chain2][chain1]
                
            if ipsa_pae_d0chn_asym[chain1][chain2] >= ipsa_pae_d0chn_asym[chain2][chain1]:
                ipsa_pae_d0chn_max[chain1][chain2]=ipsa_pae_d0chn_asym[chain1][chain2]
                ipsa_pae_d0chn_maxres[chain1][chain2]=ipsa_pae_d0chn_asymres[chain1][chain2]
                ipsa_pae_d0chn_max[chain2][chain1]=ipsa_pae_d0chn_asym[chain1][chain2]
                ipsa_pae_d0chn_maxres[chain2][chain1]=ipsa_pae_d0chn_asymres[chain1][chain2]
            else:
                ipsa_pae_d0chn_max[chain1][chain2]=ipsa_pae_d0chn_asym[chain2][chain1]
                ipsa_pae_d0chn_maxres[chain1][chain2]=ipsa_pae_d0chn_asymres[chain2][chain1]
                ipsa_pae_d0chn_max[chain2][chain1]=ipsa_pae_d0chn_asym[chain2][chain1]
                ipsa_pae_d0chn_maxres[chain2][chain1]=ipsa_pae_d0chn_asymres[chain2][chain1]
                
            if ipsa_pae_d0dom_asym[chain1][chain2] >= ipsa_pae_d0dom_asym[chain2][chain1]:
                ipsa_pae_d0dom_max[chain1][chain2]=ipsa_pae_d0dom_asym[chain1][chain2]
                ipsa_pae_d0dom_maxres[chain1][chain2]=ipsa_pae_d0dom_asymres[chain1][chain2]
                ipsa_pae_d0dom_max[chain2][chain1]=ipsa_pae_d0dom_asym[chain1][chain2]
                ipsa_pae_d0dom_maxres[chain2][chain1]=ipsa_pae_d0dom_asymres[chain1][chain2]
                n0dom_max[chain1][chain2]=n0dom[chain1][chain2]
                n0dom_max[chain2][chain1]=n0dom[chain1][chain2]
                d0dom_max[chain1][chain2]=d0dom[chain1][chain2]
                d0dom_max[chain2][chain1]=d0dom[chain1][chain2]
            else:
                ipsa_pae_d0dom_max[chain1][chain2]=ipsa_pae_d0dom_asym[chain2][chain1]
                ipsa_pae_d0dom_maxres[chain1][chain2]=ipsa_pae_d0dom_asymres[chain2][chain1]
                ipsa_pae_d0dom_max[chain2][chain1]=ipsa_pae_d0dom_asym[chain2][chain1]
                ipsa_pae_d0dom_maxres[chain2][chain1]=ipsa_pae_d0dom_asymres[chain2][chain1]
                n0dom_max[chain1][chain2]=n0dom[chain2][chain1]
                n0dom_max[chain2][chain1]=n0dom[chain2][chain1]
                d0dom_max[chain1][chain2]=d0dom[chain2][chain1]
                d0dom_max[chain2][chain1]=d0dom[chain2][chain1]
                
            if ipsa_pae_d0res_asym[chain1][chain2] >= ipsa_pae_d0res_asym[chain2][chain1]:
                ipsa_pae_d0res_max[chain1][chain2]=ipsa_pae_d0res_asym[chain1][chain2]
                ipsa_pae_d0res_maxres[chain1][chain2]=ipsa_pae_d0res_asymres[chain1][chain2]
                ipsa_pae_d0res_max[chain2][chain1]=ipsa_pae_d0res_asym[chain1][chain2]
                ipsa_pae_d0res_maxres[chain2][chain1]=ipsa_pae_d0res_asymres[chain1][chain2]
                n0res_max[chain1][chain2]=n0res[chain1][chain2]
                n0res_max[chain2][chain1]=n0res[chain1][chain2]
                d0res_max[chain1][chain2]=d0res[chain1][chain2]
                d0res_max[chain2][chain1]=d0res[chain1][chain2]
            else:
                ipsa_pae_d0res_max[chain1][chain2]=ipsa_pae_d0res_asym[chain2][chain1]
                ipsa_pae_d0res_maxres[chain1][chain2]=ipsa_pae_d0res_asymres[chain2][chain1]
                ipsa_pae_d0res_max[chain2][chain1]=ipsa_pae_d0res_asym[chain2][chain1]
                ipsa_pae_d0res_maxres[chain2][chain1]=ipsa_pae_d0res_asymres[chain2][chain1]
                n0res_max[chain1][chain2]=n0res[chain2][chain1]
                n0res_max[chain2][chain1]=n0res[chain2][chain1]
                d0res_max[chain1][chain2]=d0res[chain2][chain1]
                d0res_max[chain2][chain1]=d0res[chain2][chain1]
                
        
chaincolor={'A':'magenta', 'B':'marine', 'C':'lime', 'D':'orange', 'E':'yellow', 'F':'cyan'}


#   Chn1 Chn2 PAE tm_af2_d0num  tm_pae_d0chn  af2_d0num  pae_d0chn  pae_d0dom  pae_d0res  numres n0chn  n0dom  n0res  d0num   d0chn   d0dom   d0res   nres_a  nres_b  dist_a  dist_b pdb
#   A    B    10   0.8776        0.8736        0.8776     0.8736     0.8460     0.7780     236    236    182    114    5.70    5.70    5.03    3.94     68     114      42      58   7QII_YSCX_YEREN_YSCY_YEREN_rank024_iptm0.88_ptm0.73_nrx_3_004_10_10
#   B    A    10   0.5391        0.5182        0.8461     0.8275     0.8008     0.6459     236    236    185     73    5.70    5.70    5.07    3.00    111      74      58      41   7QII_YSCX_YEREN_YSCY_YEREN_rank024_iptm0.88_ptm0.73_nrx_3_004_10_10

OUT.write("\nChn1 Chn2 PAE Size    tm_af2_d0num  tm_pae_d0chn  af2_d0num  pae_d0chn  pae_d0dom  pae_d0res  numres n0chn  n0dom  n0res  d0num   d0chn   d0dom   d0res   nres_a  nres_b  dist_a  dist_b  pdb\n")
for chain1 in unique_chains:
    color1=chaincolor[chain1]
    for chain2 in unique_chains:
        color2=chaincolor[chain2]
        if chain1 != chain2:
            residues_1 = len(unique_residues_chain1[chain1][chain2])
            residues_2 = len(unique_residues_chain2[chain1][chain2])
            dist_residues_1 = len(dist_unique_residues_chain1[chain1][chain2])
            dist_residues_2 = len(dist_unique_residues_chain2[chain1][chain2])
            pairs = valid_pair_counts[chain1][chain2]
            dist_pairs = dist_valid_pair_counts[chain1][chain2]

            if ipsa_pae_d0res_asym[chain1][chain2] >= ipsa_pae_d0res_asym[chain2][chain1]: size="big"
            else: size="small"

            outstring=f'{chain1}    {chain2}   {int(pae_cutoff):3d}  {size:5}   ' + (
            f'{iptm_af2_d0num_asym[chain1][chain2]:8.6f}      '
            f'{iptm_pae_d0chn_asym[chain1][chain2]:8.6f}      '
            f'{ipsa_af2_d0num_asym[chain1][chain2]:8.6f}   '
            f'{ipsa_pae_d0chn_asym[chain1][chain2]:8.6f}   '
            f'{ipsa_pae_d0dom_asym[chain1][chain2]:8.6f}   '
            f'{ipsa_pae_d0res_asym[chain1][chain2]:8.6f}   '
            f'{int(numres):5d}  '
            f'{int(n0chn[chain1][chain2]):5d}  '
            f'{int(n0dom[chain1][chain2]):5d}  '
            f'{int(n0res[chain1][chain2]):5d}  '
            f'{d0num:6.2f}  '
            f'{d0chn[chain1][chain2]:6.2f}  '
            f'{d0dom[chain1][chain2]:6.2f}  '
            f'{d0res[chain1][chain2]:6.2f}  '
            f'{residues_1:5d}   '
            f'{residues_2:5d}   '
            f'{dist_residues_1:5d}   '
            f'{dist_residues_2:5d}   '
            f'{path_stem}\n')
            OUT.write(outstring)
            PML.write("#" + outstring)
            size="max"
            if chain1 > chain2:
                residues_1 = max(len(unique_residues_chain2[chain1][chain2]), len(unique_residues_chain1[chain2][chain1]))
                residues_2 = max(len(unique_residues_chain1[chain1][chain2]), len(unique_residues_chain2[chain2][chain1]))
                dist_residues_1 = max(len(dist_unique_residues_chain2[chain1][chain2]), len(dist_unique_residues_chain1[chain2][chain1]))
                dist_residues_2 = max(len(dist_unique_residues_chain1[chain1][chain2]), len(dist_unique_residues_chain2[chain2][chain1]))
                
                outstring=f'{chain2}    {chain1}   {int(pae_cutoff):3d}  {size:5}   ' + (
                    f'{iptm_af2_d0num_max[chain1][chain2]:8.6f}      '
                    f'{iptm_pae_d0chn_max[chain1][chain2]:8.6f}      '
                    f'{ipsa_af2_d0num_max[chain1][chain2]:8.6f}   '
                    f'{ipsa_pae_d0chn_max[chain1][chain2]:8.6f}   '
                    f'{ipsa_pae_d0dom_max[chain1][chain2]:8.6f}   '
                    f'{ipsa_pae_d0res_max[chain1][chain2]:8.6f}   '
                    f'{int(numres):5d}  '
                    f'{int(n0chn[chain1][chain2]):5d}  '
                    f'{int(n0dom_max[chain1][chain2]):5d}  '
                    f'{int(n0res_max[chain1][chain2]):5d}  '
                    f'{d0num:6.2f}  '
                    f'{d0chn[chain1][chain2]:6.2f}  '
                    f'{d0dom_max[chain1][chain2]:6.2f}  '
                    f'{d0res_max[chain1][chain2]:6.2f}  '
                    f'{residues_1:5d}   '
                    f'{residues_2:5d}   '
                    f'{dist_residues_1:5d}   '
                    f'{dist_residues_2:5d}   '
                    f'{path_stem}\n')
                OUT.write(outstring)
                PML.write("#" + outstring)
                
            chain_pair= f'color_{chain1}_{chain2}'
            chain1_residues = f'chain  {chain1} and resi {contiguous_ranges(unique_residues_chain1[chain1][chain2])}'
            chain2_residues = f'chain  {chain2} and resi {contiguous_ranges(unique_residues_chain2[chain1][chain2])}'
            PML.write(f'alias {chain_pair}, color {color1}, {chain1_residues}; color {color2}, {chain2_residues}\n\n')
            
