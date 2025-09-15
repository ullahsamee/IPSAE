# ipSAE: for binders screening
![ipSAE Animation](ipsae.gif)

## Paper: What’s wrong with AlphaFold’s 𝑖𝑝𝑇𝑀 score and how to fix it https://doi.org/10.1101/2025.02.10.637595
## Installation

```bash
pip install ipsae
```

## Usage

Run `ipSAE` tool from the command line terminal in Linux. It takes the Alphafold2, Alphold3 and Boltz2 multimers PAE file in .json and .pdb or .cif file format as input, and distance cutoff values.

### Command-line Help

```
$ ipsae --help
```
```
ipSAE
Scoring function for interprotein interactions
-------------------------------------------------
Version: 1.0
Author: Samee Ullah
Email: sameeullah@bs.qau.edu.pk
LinkedIn: https://www.linkedin.com/in/ullah-samee/
-------------------------------------------------
Based on the original work of Roland Dunbrack.

usage: ipsae [-h] pae_file structure_file [cutoffs ...]

Calculate ipSAE scores for protein-protein interactions

positional arguments:
  pae_file        PAE file (JSON or NPZ)
  structure_file  Structure file (PDB or CIF)
  cutoffs         Optional: PAE and distance cutoff values. Defaults are 15
                  for AF2 and 10 for AF3/Boltz1.

options:
  -h, --help      show this help message and exit

Examples:
For AlphaFold2:
  ipsae alphafold2_multimer.json alphafold2_model.pdb 15 15

For AlphaFold3:
  ipsae full_data_0.json model_0.cif 10 10

For Boltz:
  ipsae pae_model_0.npz Boltz_model.cif 10 10
```

### Examples

#### AlphaFold2

```bash
ipsae scores.json model.pdb 15 15
```
default distance cutoffs of 15.

#### AlphaFold3

```bash
ipsae full_data_0.json model_0.cif 10 10
```
default distance cutoffs of 10.

#### Boltz1

```bash
ipsae pae_model_0.npz model_0.cif 10 10
```
default distance cutoffs of 10.

## Output Files

The tool generates three output files:

1.  **`<name>_<pae>_<dist>.txt`**: The main output file with chain-chain scores.
2.  **`<name>_<pae>_<dist>_byres.txt`**: By-residue scores.
3.  **`<name>_<pae>_<dist>.pml`**: A PyMOL script for visualizing the ipSAE scores.

### Main Output File Format

The main output file contains a table with the following columns:

| Chn1 | Chn2 | PAE | Dist | Type | ipSAE | ipSAE_d0chn | ipSAE_d0dom | ipTM_af | ipTM_d0chn | pDockQ | pDockQ2 | LIS | n0res | n0chn | n0dom | d0res | d0chn | d0dom | nres1 | nres2 | dist1 | dist2 | Model |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

## Credits

This tool is based on the original work of Roland Dunbrack.
