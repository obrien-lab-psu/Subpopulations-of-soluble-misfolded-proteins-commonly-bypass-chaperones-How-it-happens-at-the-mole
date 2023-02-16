Quanitfying the significance of the consistency between simualted post-translational structural distibution and LiP-MS experimental data

(1) Use the workflow at the following link [https://github.com/obrien-lab/cg_simtk_protein_folding] to take a set of coarse grained protein dynatmics simulations and generate a set of metastable macrostates. 

(2) Then use the workflow here [https://github.com/obrien-lab-psu/entanglement_analysis_public] to identify the set of unique misfolded conformations that have changes in self-entanglement in each metastable state

(3) Get distribution of theoretical proteinase K cut sites
Usage of codes/get_theoretical_LiP_MS_peptides_v2.1.py
        python codes/get_theoretical_LiP_MS_peptides_v2.1.py [1] [2] [3] [4] [5] [6] [7] [8]
        [1] path to all-atom reference pdb
        [2] path to join prob pkl file
        [3] outfile_basename
        [4] num iterations
        [5] nproc
        [6] path to PKcutsite observed probability
        [7] path to observed AA prob across proteome
        [8] out_path

        Examples
        python codes/get_theoretical_LiP_MS_peptides_v2.0.py 2fym/2fym_chain_a_renum_rebuilt_mini_v3.pdb whole_proteome_joint_prob.pkl 2fym/p2fym 100000000 10 whole_proteome_obs_PKcutsite_prob.txt AA_prob_across_pdb_v1.1.txt
        python codes/get_theoretical_LiP_MS_peptides_v2.0.py 1p7l/1p7l_chain_a_renum_rebuilt_mini_v3.pdb whole_proteome_joint_prob.pkl 1p7l/p1p7l 100000000 10 whole_proteome_obs_PKcutsite_prob.txt AA_prob_across_pdb_v1.1.txt


(4) Backmap representative entangled structures and get per residue SASA changes relative to the native state

    (4a) backmap each frame of the DCD to the all-atom structure
            python codes/backmap_dcd_v1.1.py
            [1] = dcd
            [2] = psf
            [3] = working dir
            [4] = outfile basename
            [5] = AA PDB

    (4b) collect the backmapped PDBs into a single dcd
            python codes/collect_backmapped_frames_to_dcd_v1.1.py
            [1] path to input backmapped pdbs files
            [2] path to outfile


    (4c) calculate the solvant accessible surface area with the Shrake Rupley method
            Usage of codes/get_SASA_change_v1.2.py
            python codes/get_SASA_change_v1.2.py [1] [2] [3] [4] [5]
            [1] path to all-atom reference pdb
            [2] outfile basename
            [3] path to backmapped dcd
            [4] path to backmapped dcd topology
            [5] out_path

    Examples
    python codes/get_SASA_change_v1.2.py inpfiles/1zmr_model_clean.pdb 1zmr inpfiles/sampled_traj_AA.dcd inpfiles/1zmr_model_clean.pdb ./


(5) permutation test to determine the probability of our model being this consistent given random chance
    python codes/overlap_stat_test_v4.0.py [1] [2] [3] [4] [5] [6] [7] [8] [9] [10]

    [1] path to ent  state .pkl file
    [2] path to LiPMS_peptide file
    [3] path to outfile
    [4] residue buffer for residues considered around PK cutsite when doing overlap analysis
    [5] path to theoretical peptide file containing list of theoretical PK cut sites
    [6] path to sasa change file
    [7] max_res in protein model
    [8] reps for permutation 
    [9] number processors


    Examples
    python codes/overlap_stat_test_v4.0.py 1p7l/p1p7l_u_ent_states.pkl 1p7l/p1p7l_LiPMS_peptides_Cbuff.csv 1p7l/p1p7l_overlap_v4.0Cbuff_5rb_final 5 1p7l/p1p7l_any_pk_site_data.txt 1p7l/p1p7l_sasa_change.pkl 384 100000 20
    python codes/overlap_stat_test_v4.0.py 2fym/p2fym_u_ent_states.pkl 2fym/p2fym_LiPMS_peptides_Cbuff.csv 2fym/p2fym_overlap_v4.0Cbuff_5rb_final 5 2fym/p2fym_any_pk_site_data.txt 2fym/p2fym_sasa_change.pkl 432 100000 20
