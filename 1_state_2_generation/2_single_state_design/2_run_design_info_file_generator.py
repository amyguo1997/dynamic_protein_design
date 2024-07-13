# -*- coding: utf-8 -*-
"""
implement sequence_design
"""
  
#!/usr/bin/env python3

import os
import time
import json

import pyrosetta
from pyrosetta import rosetta

import design_info_file_generator as LPSD

def get_bb_remodeled_residues_for_LHL_designs(file_for_insertion_points):
    '''Get the backbone remodeled residues for models made by
    LHL reshaping. The two achoring residues on the flanking
    secondary structures are included.
    '''
    with open(file_for_insertion_points, 'r') as f:
        insertion_points = json.load(f)

    bb_remodeled_residues = []

    for ip in insertion_points:
        start = ip['start'] + 1
        stop = ip['stop'] - 1 

        for i in range(start, stop + 1):
            bb_remodeled_residues.append(i)

    return bb_remodeled_residues

def design(input_dir, model_id, model_path, pre_moved_bb_pdb, file_for_pre_moved_bb_insertion_points, residues_to_fix):
         
    # Get the pose for the pre-moved structure for the remodeled residues. 

    pre_moved_pose_whole = rosetta.core.import_pose.pose_from_file(pre_moved_bb_pdb)
    pre_moved_bb_insertion_points = get_bb_remodeled_residues_for_LHL_designs(file_for_pre_moved_bb_insertion_points)
                
    v_pre_moved_bb_insertion_points = rosetta.utility.vector1_unsigned_long()
    for p in pre_moved_bb_insertion_points:
        v_pre_moved_bb_insertion_points.append(p)
               
    pre_moved_bb_pose = rosetta.core.pose.Pose()
    rosetta.core.pose.pdbslice(pre_moved_bb_pose, pre_moved_pose_whole, v_pre_moved_bb_insertion_points)

    # Get the bb remodeled residues

    bb_remodeled_residues = get_bb_remodeled_residues_for_LHL_designs(file_for_pre_moved_bb_insertion_points)
    LPSD.make_one_design(model_path, os.path.join(input_dir,'model_%s.pdb.gz'%model_id), bb_remodeled_residues, pre_moved_bb_pose=pre_moved_bb_pose, residues_to_fix=residues_to_fix)

if __name__ == '__main__':
    
    pyrosetta.init()

    residues_to_fix = [66,68,69,70,71,72,77]
    pre_moved_bb_pdb = os.path.expanduser('~/src/loop-helix-loop-reshaping/data/1smg_cleaned.pdb')
    file_for_pre_moved_bb_insertion_points = os.path.expanduser('~/src/loop-helix-loop-reshaping/data/insertion_points.json')
    
    # Paths
    model_rmsd_path = os.path.expanduser('~/data/rmsd.json') # File with RMSD values
    output_path = os.path.expanduser('~/data/outputs')
    lucs_models_path = os.path.join('~/lucs_models')
    
    ## Get the model number for parallelization on CPUs
    model_ids = []
    with open(model_rmsd_path,'r') as rf:
        model_rmsd_data = json.load(rf)
        
    for model_id, rmsd in model_rmsd_data.items():
        if rmsd < 10 and rmsd > 3:
            model_ids.append(model_id)
    
    task_id = int(os.environ['SGE_TASK_ID'])
    model_id = model_ids[task_id - 1]
    model_path = os.path.join(output_path, model_id)
    
    design(lucs_models_path, model_id, model_path, pre_moved_bb_pdb, file_for_pre_moved_bb_insertion_points, residues_to_fix)
