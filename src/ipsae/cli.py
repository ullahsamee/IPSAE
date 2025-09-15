"""
Command-line interface for ipSAE.
"""

import argparse
import sys
from pathlib import Path
from .core import calculate_ipsae


def main():
    """Main CLI entry point."""
    
    header = """
ipSAE
Scoring function for interprotein interactions
-------------------------------------------------
Version: 1.0
Author: Samee Ullah
Email: sameeullah@bs.qau.edu.pk
LinkedIn: https://www.linkedin.com/in/ullah-samee/
-------------------------------------------------
Based on the original work of Roland Dunbrack.
"""

    examples = """
Examples:
For AlphaFold2:
  ipsae alphafold2_multimer.json alphafold2_model.pdb 15 15

For AlphaFold3:
  ipsae full_data_0.json model_0.cif 10 10

For Boltz:
  ipsae pae_model_0.npz Boltz_model.cif 10 10
"""

    if '--help' in sys.argv or '-h' in sys.argv:
        print(header)
        print(examples)
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="Calculate ipSAE scores for protein-protein interactions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=examples,
        add_help=False # disable the default help
    )
    
    parser.add_argument("pae_file", help="PAE file (JSON or NPZ)")
    parser.add_argument("structure_file", help="Structure file (PDB or CIF)")
    parser.add_argument("cutoffs", type=float, nargs='*', help="Optional: PAE and distance cutoff values. Defaults are 15 for AF2 and 10 for AF3/Boltz1.")
    
    args = parser.parse_args()
    
    pae_cutoff = None
    dist_cutoff = None

    if len(args.cutoffs) == 2:
        pae_cutoff = args.cutoffs[0]
        dist_cutoff = args.cutoffs[1]
    elif len(args.cutoffs) == 0:
        # Determine file type to set default cutoffs
        if ".pdb" in args.structure_file:
            pae_cutoff = 15.0
            dist_cutoff = 15.0
        elif ".cif" in args.structure_file:
            pae_cutoff = 10.0
            dist_cutoff = 10.0
        else:
            # Default to 10 if structure file type is ambiguous
            pae_cutoff = 10.0
            dist_cutoff = 10.0
    else:
        parser.error("Please provide two cutoff values (PAE and distance) or none to use defaults.")

    try:
        results = calculate_ipsae(
            args.pae_file,
            args.structure_file, 
            pae_cutoff,
            dist_cutoff
        )
        print("Calculation completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
