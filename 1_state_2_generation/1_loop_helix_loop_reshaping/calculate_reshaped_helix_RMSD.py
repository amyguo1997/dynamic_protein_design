# -*- coding: utf-8 -*-
'''
Sort by length and geometric variance
'''
import os 
import json
import numpy as np
import sys
import time

from pyrosetta import *
from pyrosetta.teaching import *

def get_helical_residues(pose, reshaped_residues):
    
    '''
    LHL_residues should be for one LHL subunit 
    '''
    
    # Non-helix residues to exclude (both loops, or one loop/one sheet):
    helical_residues = rosetta.utility.vector1_unsigned_long()
    
    for residue in reshaped_residues:
        if (pose.secstruct(residue) == 'H'):
            helical_residues.append(residue)
                  
    return helical_residues

def RMSD(pose1, pose2, helical_residues):
    '''Calcualte RMSD between two lists of numpy points.'''
    
    points1 = []
    points2 = []
    
    for residue in helical_residues:
        points1.append(pose1.residue(residue).xyz("CA"))
        points2.append(pose2.residue(residue).xyz("CA"))
        
    diff = []
    for i in range(0, len(points1)):
        diff.append(points1[i] - points2[i])
        
    rmsd = np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))  
    
    return(rmsd)

def align(pose1, pose2, bb_remodeled_residues):
    
    '''
    Align constant region atoms
    '''
    
    all_residues = range(1,pose1.total_residue()+1)
    
    # Get constant residues 
    constant_residues = rosetta.utility.vector1_unsigned_long()
    for residue in all_residues:
        if residue not in bb_remodeled_residues:
            constant_residues.append(residue)
            
    # Align CA atoms of constant region
    pyrosetta.rosetta.protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(pose1, pose2, constant_residues)
       
    
if __name__ == '__main__':
    
    pyrosetta.init("-mute all")
    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
    
    output_file = os.path.expanduser('~/data/rmsd.json') # Desired output file name and location         
    
    reference_pdb_path = os.path.expanduser('~/src/loop-helix-loop-reshaping/data/1smg_cleaned.pdb') # Reference pose
    reference_pose = pose_from_pdb(reference_pdb_path)
    DSSP.apply(reference_pose)

    # Get reshaped helix residues
    insertion_points_file = os.path.expanduser('~/src/loop-helix-loop-reshaping/data/insertion_points.json') # Insertion point file for LUCS
    with open(insertion_points_file, 'r') as rf:
        LHL_subunits = json.load(rf)

    start = LHL_subunits[0]['start']
    stop = LHL_subunits[0]['stop']
    reshaped_residues = list(range(start+1, stop))
    
    helical_residues = get_helical_residues(reference_pose, reshaped_residues)
        
    model_data = []
    data_path = os.path.expanduser('~/lucs_models') # Directory with loop-helix-loop reshaping output models model_<model_id>.pdb.gz
    os.chdir(data_path)
    
    for root, dirs, files in os.walk(data_path):
        for fname in files:
            if fname.endswith('pdb.gz'):
                model_id = fname.split('.')[0].split('_')[-1]
                reshaped_pose = pose_from_pdb(fname)
                DSSP.apply(reshaped_pose)
                
                # Superimpose on non-reshaped backbone
                align(reference_pose, reshaped_pose, reshaped_residues)
                
                # Calculate RMSD for all reshaped LHL subunits
                rmsd = RMSD(reference_pose, reshaped_pose, helical_residues)
                model_data[model_id] = rmsd

            
    with open(output_file,'w') as wf:
        json.dump(data, wf)