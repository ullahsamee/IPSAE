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

    Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn    n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model
    A    B     15   15   asym  0.563173    0.824097    0.678527    0.460    0.362276      292   1571    541    6.28   12.57    8.21    249     292      90      83   fold_raf1_ksr1_mek1_model_0
    B    A     15   15   asym  0.546718    0.818835    0.671409    0.460    0.442774      278   1571    539    6.14   12.57    8.20    259     280      83      92   fold_raf1_ksr1_mek1_model_0
    A    B     15   15   max   0.563173    0.824097    0.678527    0.460    0.442774      292   1571    541    6.28   12.57    8.21    280     292      92      83   fold_raf1_ksr1_mek1_model_0

    A    C     15   15   asym  0.261359    0.484143    0.338764    0.510    0.428632      319   1041    491    6.54   10.71    7.88    171     320       0       0   fold_raf1_ksr1_mek1_model_0
    C    A     15   15   asym  0.232004    0.486093    0.350874    0.510    0.270827      260   1041    513    5.96   10.71    8.03    250     263       0       0   fold_raf1_ksr1_mek1_model_0
    A    C     15   15   max   0.261359    0.486093    0.350874    0.510    0.428632      319   1041    513    6.54   10.71    8.03    263     320       0       0   fold_raf1_ksr1_mek1_model_0

    B    C     15   15   asym  0.636110    0.829537    0.732184    0.770    0.751586      344   1316    607    6.76   11.74    8.61    262     345      69      72   fold_raf1_ksr1_mek1_model_0
    C    B     15   15   asym  0.559980    0.802212    0.692571    0.770    0.343057      291   1316    598    6.27   11.74    8.56    307     291      71      69   fold_raf1_ksr1_mek1_model_0
    B    C     15   15   max   0.636110    0.829537    0.732184    0.770    0.751586      344   1316    607    6.76   11.74    8.61    291     345      69      72   fold_raf1_ksr1_mek1_model_0




**Chn1**=first chain

**Chn2**=second chain

**PAE**=PAE cutoff

**Dist**=Distance cutoff for CA-CA contacts

**Type**="asym" or "max"; asym means asymmetric ipTM/ipSAE values; max is maximum of asym values

**ipTM_af**=AlphaFold ipTM values. For AF2, this is for whole complex. For AF3, this is symmetric pairwise value.   

**ipSAE**=ipSAE value for given PAE cutoff and d0 determined by number of residues in 2nd chain with PAE<cutoff 

**ipTM_d0chn**=ipTM calculated from PAE matrix and d0 = sum of chain lengths 

**ipSAE_d0chn**=ipSAE calculated with PAE cutoff and d0 = sum of chain lengths

**ipSAE_d0dom**=ipSAE calculated with PAE cutoff and d0 = total number of residues in both chains with any interchain PAE<cutoff

**n0res**=number of residues for d0 in ipSAE calculation

**n0chn**=number of residues in d0 in ipSAE_d0chn calculation

**n0dom**=number of residues in d0 in ipSAE_d0dom calculation

**d0res**=d0 for ipSAE

**d0chn**=d0 for ipSAE_d0chn

**d0dom**=d0 for ipSAE_d0dom

**nres1**=number of residues in chain1 with PAE<cutoff with residues in chain2

**nres2**=number of residues in chain2 with PAE<cutoff with residues in chain1

**dist1**=number of residues in chain 1 with PAE<cutoff and dist<cutoff from chain2

**dist2**=number of residues in chain 2 with PAE<cutoff and dist<cutoff from chain1

**pdb**=AlphaFold filename


