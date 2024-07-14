#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Select switch designs from ProteinMPNN MSD
"""

# Import libraries
import json
import numpy as np
import os

# Based on Dayhoff classification
RES_TO_INT = {'A': 1, # Small 
              'G': 2, # Small, very flexible
              'P': 3, # Special backbone geometry
              'C': 4, # Disulfide forming
              'S': 5, 'T': 5, # Small, OH containing
              'D': 6, 'E': 6, 'N': 6, 'Q': 6, # Acids and amides
              'F': 7, 'Y': 7, 'W': 7, # Aromatic
              'H': 7, # pH sensitive at physiological pH
              'K': 8, 'R': 8, # Large bases
              'M': 9, 'L': 9, 'V': 9, 'I': 9, # Hydrophobic 
              '-': 21, '.': 21, '~': 21, 'X': 21 # Unknown or gap
              }

def find_ss_similar(summarized_data_path, fasta_dir, state_1_rmsd_threshold, state_1_var_threshold, state_1_plddt_threshold, 
                    state_2_rmsd_threshold, state_2_var_threshold, state_2_plddt_threshold, state_id = 'LUCS'):
    
    # Find sequences adopting a single state that are similar in seqid to switch designs

    with open(summarized_data_path, 'r') as rf:
        summarized_data = json.load(rf)
    
    low_variance_high_plddt_state_1_designs = {}
    low_variance_high_plddt_state_2_designs = {}
         
    for design_id in summarized_data:

        variance = np.var([summarized_data[design_id]['%s_rmsd_bb_remodeled_%d'%(state_id, i)] for i in range (1,6)])
        best_plddt = max([summarized_data[design_id]['pLDDT_bb_remodeled_%d'%(i)] for i in range (1,6)])
        state_1_RMSDs = [summarized_data[design_id]['%s_rmsd_bb_remodeled_%d'%('WT', i)] for i in range (1,6)]
        state_2_RMSDs = [summarized_data[design_id]['%s_rmsd_bb_remodeled_%d'%('LUCS', i)] for i in range (1,6)]
       
        if any (rmsd < state_1_rmsd_threshold for rmsd in state_1_RMSDs) :
            if variance < state_1_var_threshold and best_plddt > state_1_plddt_threshold:
                fasta = os.path.join(fasta_dir, '%s.fasta'%design_id)
                with open(fasta, 'r') as rf:
                    sequence = list(rf.readlines()[-1].rstrip())
                low_variance_high_plddt_state_1_designs[design_id] = sequence
                
        if any (rmsd < state_2_rmsd_threshold for rmsd in state_2_RMSDs):
            if variance < state_2_var_threshold and best_plddt > state_2_plddt_threshold:
                fasta = os.path.join(fasta_dir, '%s.fasta'%design_id)
                with open(fasta, 'r') as rf:
                    sequence = list(rf.readlines()[-1].rstrip())
                low_variance_high_plddt_state_2_designs[design_id] = sequence
        
              
    for ss_id_1, ss_design_1 in low_variance_high_plddt_state_1_designs.items():
        for ss_id_2, ss_design_2 in low_variance_high_plddt_state_2_designs.items():
                count = 0
                for position in range(len(ss_design_1)):
                    if [ss_design_1[position-1]] != [ss_design_2[position-1]]:
                        count += 1.0
                if count < 6:
                    print('Single state ID WT:', ss_id_1,
                          'Single state ID LUCS:', ss_id_2,
                          'Differences:', count)
    
def find_single_state_designs(summarized_data_path, rmsd_threshold = 2.0, variation_threshold = 0.1 ):
    
    '''
    Find designs where all models < rmsd threshold for a given state and return confidence
    '''

    with open(summarized_data_path, 'r') as rf:
        summarized_data = json.load(rf)
    
    single_state_designs = []
    
    for design_id in summarized_data:
        plddt_of_lucs_helix = [summarized_data[design_id]['pLDDT_bb_remodeled_%d'%(i)] for i in range (1,6)]
        best_by_plddt = plddt_of_lucs_helix.index(max(plddt_of_lucs_helix))+1
        rmsd_of_best_by_plddt_to_state = summarized_data[design_id]['rmsd_bb_remodeled_%d'%(best_by_plddt)]
        var = np.var([summarized_data[design_id]['rmsd_bb_remodeled_%d'%(i)] for i in range (1,6)])
        
        if rmsd_of_best_by_plddt_to_state < rmsd_threshold and var < variation_threshold:
            single_state_designs.append({'design_id': design_id,
                                    'rmsd': rmsd_of_best_by_plddt_to_state,
                                    'variance': var,
                                    'plddt': max(plddt_of_lucs_helix),
                                    'best model': best_by_plddt})

    return single_state_designs

if __name__ == '__main__':
    
    # The af directory should have subdirectories containing the *.pdb AF2 predictions of each design
    af_dir = '$PATH_TO_AF2_PREDICTIONS'
    # The fasta directory should contain all *.fasta files for each switch design run through AF2
    fasta_dr = '$PATH_TO_FASTA_FILES'
    summarized_data_path = os.path.join(af_dir, 'summarized_structure_pred_eval.json')
    
