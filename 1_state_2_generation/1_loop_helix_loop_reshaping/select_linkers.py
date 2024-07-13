#!/usr/bin/env python3

import os
import sys
import json
import copy

import pyrosetta
from pyrosetta import rosetta

import loop_helix_loop_reshaping as LHLR


def select_and_dump_linkers(input_pdb, input_database, output_database, linker_length, insertion_points, insertion_id, front_linker):
    # Prepare the pose
    
    pose = rosetta.core.import_pose.pose_from_file(input_pdb)
    LHLR.simple_pose_moves.remove_insertion_residues(pose, insertion_points)
    LHLR.select_linkers.prepare_linker_selection(pose, linker_length, 
            insertion_points[insertion_id]['start'], insertion_points[insertion_id]['stop'], front_linker)

    print('select linkers from', input_database)

    with open(input_database, 'r') as f:
        candidate_linkers = json.load(f)

    # Select linkers

    selected_linkers = LHLR.select_linkers.select_non_clashing_linkers(pose, candidate_linkers, insertion_points[insertion_id]['start'])
    
    with open(output_database, 'w') as f:
        json.dump(selected_linkers, f)
    
    print('dump selected linkers to', output_database)

def select_linkers(data_path, input_pdb, input_insertion_points_file, insertion_ids_for_selection=None):
    pyrosetta.init()
    
    # Load insertion points
    
    with open(input_insertion_points_file, 'r') as f:
        insertion_points = json.load(f)

    if insertion_ids_for_selection is None:
        insertion_ids_for_selection = range(len(insertion_points))

    for insertion_id in insertion_ids_for_selection:
        start_ss = insertion_points[insertion_id]['start_ss']
        stop_ss = insertion_points[insertion_id]['stop_ss']

        for linker_length in [2, 3, 4, 5]:

            # Select front linkers

            input_database = 'database/linker_database/linker_{0}_helix_{1}_non_redundant.json'.format(start_ss, linker_length)
            output_database = os.path.join(data_path, 'selected_linkers_{0}_{1}_front.json'.format(insertion_id, linker_length))
            select_and_dump_linkers(input_pdb, input_database, output_database, linker_length, copy.deepcopy(insertion_points), insertion_id, True)
    
            # Select back linkers  

            input_database = 'database/linker_database/linker_helix_{0}_{1}_non_redundant.json'.format(stop_ss, linker_length)
            output_database = os.path.join(data_path, 'selected_linkers_{0}_{1}_back.json'.format(insertion_id, linker_length))
            select_and_dump_linkers(input_pdb, input_database, output_database, linker_length, copy.deepcopy(insertion_points), insertion_id, False)
            

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1
  
    input_pdb = 'data/1smg_cleaned.pdb'
    input_insertion_points_file = 'data/insertion_points.json'
    select_linkers(data_path, input_pdb, input_insertion_points_file, insertion_ids_for_selection=None)
