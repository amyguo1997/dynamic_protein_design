#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate RMSD of reshaped helix of AF2 models to state 1 and state 2
"""

import numpy as np
import json
import argparse
import os

from pyrosetta import *
from pyrosetta import rosetta


def RMSD(coords_1, coords_2):
    
    '''
    Return RMSD given two sets of xyz coordinates
    '''

    differences = [coords_1[i] - coords_2[i] for i in range(0, len(coords_1))]
    rmsd = np.sqrt(sum(np.dot(d,d) for d in differences)/len(differences))
    
    return rmsd

def align(pose, reference_pose, helical_residues):
    
    '''
    Superimpose helical test_pose residues onto state_2_pose
    '''
    
    superimpose_residues = rosetta.utility.vector1_unsigned_long() # need this type
    for residue in helical_residues:
        superimpose_residues.append(residue)
        
    rosetta.protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(pose, reference_pose, superimpose_residues)

def get_helical_residues(state_2_pose):
    
    '''
    Return helical residue positions in target pose
    '''
    
    # Get helical residues from state_2_pose remodeled region
    helical_residues = []
    
    for residue in range(1,state_2_pose.total_residue()+1):
        if state_2_pose.secstruct(residue) == 'H':
            helical_residues.append(residue)
    
    return helical_residues

def get_atom_xyz(pose):
    
    '''
    Return xyz position for CA atoms
    '''
    
    xyz_coords = [pose.residue(i).xyz('CA') for i in range(1, pose.total_residue()+1)]

    return xyz_coords

def get_bb_remodeled_residues(design_info_path):
    
    with open(design_info_path, 'r') as rf:
        data = json.load(rf)
    bb_remodeled_residues = data['bb_remodeled_residues']
    
    return bb_remodeled_residues

if __name__ == '__main__':
    
    pyrosetta.init()
    sfxn = get_fa_scorefunction()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-af_dir', type=str, dest='alphafold_directory', required=True, help='Directory of AF2 predictions')
    parser.add_argument('-state_2_pdb', type=str, dest='state_2_pdb_path', required=True, help='PDB file path of state 2 structure')
    parser.add_argument('-state_1_pdb', type=str, dest='state_1_pdb_path', required=True, help='PDB file path of state 1 structure')
    parser.add_argument('-idx_range', nargs='+', dest='idx_range', type=int, help = 'start and end range of residues to calculate overall RMSD (to remove disordered termini from calculation) e.g. 1 90')
    parser.add_argument('-design_info', type=str, dest='design_info_path', required=True, help='Design info file containing bb_remodeled_residues')

    args = parser.parse_args()
    
    # Get relevant residue ranges for RMSD calculations
    idx_range = range(args.idx_range[0], args.idx_range[1]+1)
    bb_remodeled_residues = get_bb_remodeled_residues(args.design_info_path)

    # Get ss assignments
    DSSP = rosetta.protocols.moves.DsspMover()
    state_2_pose = pose_from_pdb(args.state_2_pdb_path)
    state_1_pose = pose_from_pdb(args.state_1_pdb_path)
    DSSP.apply(state_2_pose)  
    DSSP.apply(state_1_pose)
    
    # Get helical residues (union of WT/target helical residues)
    helical_residues_state_2 = get_helical_residues(state_2_pose)
    helical_residues_state_1 = get_helical_residues(state_1_pose)
    helical_residues = list(set(helical_residues_state_2 + helical_residues_state_1))
    
    # Get CA coordinates
    state_2_coords = get_atom_xyz(state_2_pose)
    state_2_coords_overall = [state_2_coords[i-1] for i in idx_range]
    state_2_coords_bb_remodeled = [state_2_coords[i-1] for i in bb_remodeled_residues if i in helical_residues]
    state_1_coords = get_atom_xyz(state_1_pose)
    state_1_coords_overall = [state_1_coords[i-1] for i in idx_range]
    state_1_coords_bb_remodeled = [state_1_coords[i-1] for i in bb_remodeled_residues if i in helical_residues]
    
    design_id = os.environ['SGE_TASK_ID']
    
    data = {}
    
    for root, dirs, files in os.walk(os.path.join(args.alphafold_directory, design_id)):
        for fname in files:
            for i in range(1,6):
                if fname.endswith('pdb') and 'rank_%d'%i in fname:
                    
                    # Make AF2 pdb pose 
                    alphafold_pdb_path = os.path.join(root, fname)
                    af_pose = pose_from_pdb(alphafold_pdb_path)
                    
                    # Align based on helical residues and get af_pose CA coords
                    align(af_pose, state_2_pose, helical_residues)
                    af_coords = get_atom_xyz(af_pose)
                    af_coords_overall = [af_coords[i-1] for i in idx_range]
                    af_coords_bb_remodeled = [af_coords[i-1] for i in bb_remodeled_residues if i in helical_residues]
                    
                    # Calculate RMSDs
                    af_rmsd_overall = RMSD(af_coords_overall, state_2_coords_overall)
                    af_rmsd_bb_remodeled = RMSD(af_coords_bb_remodeled, state_2_coords_bb_remodeled)
                    data['state_2_rmsd_overall_%d'%i] = af_rmsd_overall
                    data['state_2_rmsd_bb_remodeled_%d'%i] = af_rmsd_bb_remodeled

                    align(af_pose, state_1_pose, helical_residues)
                    af_coords = get_atom_xyz(af_pose)
                    af_coords_overall = [af_coords[i-1] for i in idx_range]
                    af_coords_bb_remodeled = [af_coords[i-1] for i in bb_remodeled_residues if i in helical_residues]
                    
                    af_rmsd_overall = RMSD(af_coords_overall, state_1_coords_overall)
                    af_rmsd_bb_remodeled = RMSD(af_coords_bb_remodeled, state_1_coords_bb_remodeled)
                    data['state_1_rmsd_overall_%d'%i] = af_rmsd_overall
                    data['state_1_rmsd_bb_remodeled_%d'%i] = af_rmsd_bb_remodeled
                    
                    # Get pLDDT data
                    with open(alphafold_pdb_path, 'r') as rf:
                        lines = rf.readlines()
                    pLDDT = []
                    for line in lines:
                        if 'CA' in line: 
                            pLDDT.append(float(line.split()[-2]))
                    data['pLDDT_overall_%d'%i] = np.average([pLDDT[i] for i in idx_range])
                    data['pLDDT_bb_remodeled_%d'%i] = np.average([pLDDT[i] for i in bb_remodeled_residues if i in helical_residues])

    with open(os.path.join(args.alphafold_directory, design_id, 'structure_pred_eval.json'), 'w') as wf:
        json.dump(data, wf)
    
