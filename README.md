# IPSAE
Scoring function for interprotein interactions in AlphaFold2 and AlphaFold3

# Usage:                                                                                                                                                                                                     
AlphaFold2:
     python ipsae.py <path_to_json_file> <path_to_af2_pdb_file>  <pae_cutoff> <dist_cutoff>   
     python ipsae.py RAF1_KSR1_scores_rank_001_alphafold2_multimer_v3_model_4_seed_003.json RAF1_KSR1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_003.pdb 10 10

AlphaFold3:
python ipsae.py <path_to_json_file> <path_to_af3_cif_file>  <pae_cutoff> <dist_cutoff>                                                                                                                    
python ipsae.py fold_aurka_tpx2_full_data_0.json  fold_aurka_tpx2_model_0.cif 10 10
