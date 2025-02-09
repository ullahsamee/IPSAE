# IPSAE
Scoring function for interprotein interactions in AlphaFold2 and AlphaFold3

# Usage:                                                                                                                                                                                                     
AlphaFold2:

     python ipsae.py <path_to_json_file> <path_to_af2_pdb_file> <pae_cutoff> <dist_cutoff>   
     python ipsae.py RAF1_KSR1_scores_rank_001_alphafold2_multimer_v3_model_4_seed_003.json RAF1_KSR1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_003.pdb 10 10

AlphaFold3:

     python ipsae.py <path_to_json_file> <path_to_af3_cif_file> <pae_cutoff> <dist_cutoff>                    
     python ipsae.py fold_aurka_tpx2_full_data_0.json fold_aurka_tpx2_model_0.cif 10 10

# Output

AlphaFold2 complex of RAF1, KSR1, and MEK1

    Chn1 Chn2 PAE Dist   Type   ipTM_af       ipSAE  ipTM_d0chn ipSAE_d0chn ipSAE_d0dom    n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2   pdb
    A    B     15   15   asym     0.430    0.377610    0.662062    0.874327    0.763488      289   1571    550    6.25   12.57    8.27    259     291      92      89   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    B    A     15   15   asym     0.430    0.477229    0.664869    0.879583    0.767981      285   1571    541    6.21   12.57    8.21    255     286      89      92   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    A    B     15   15   max      0.430    0.477229    0.664869    0.879583    0.767981      285   1571    541    6.21   12.57    8.21    286     291      92      89   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    
    A    C     15   15   asym     0.430    0.748725    0.644246    0.797952    0.734552      362   1041    635    6.91   10.71    8.77    273     362     107      92   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    C    A     15   15   asym     0.430    0.424375    0.607214    0.807980    0.735115      285   1041    597    6.21   10.71    8.55    312     285      84     106   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    A    C     15   15   max      0.430    0.748725    0.644246    0.807980    0.735115      362   1041    597    6.91   10.71    8.55    285     362     107      92   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    
    B    C     15   15   asym     0.430    0.542510    0.326022    0.581319    0.418832      350   1316    566    6.81   11.74    8.37    213     353       0       0   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    C    B     15   15   asym     0.430    0.282392    0.316282    0.612731    0.450441      287   1316    566    6.23   11.74    8.37    278     288       0       0   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000
    B    C     15   15   max      0.430    0.542510    0.326022    0.612731    0.450441      350   1316    566    6.81   11.74    8.37    288     353       0       0   RAF1_KSR1_MEK1_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000

