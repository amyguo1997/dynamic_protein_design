#!/usr/bin/env python3

import os
import sys
import json

import pyrosetta
from pyrosetta import rosetta

import loop_helix_loop_reshaping as LHLR

def screen(data_path, linker_database_path, input_pdb, input_insertion_points_file, num_jobs, job_id, max_num_success_each_db_pair=None, insertion_ids_to_screen=None,
        num_res_clashes_tolerance=0):
  
    # Load insertion points
    
    with open(input_insertion_points_file, 'r') as f:
        insertion_points = json.load(f)
    
    # Load and clean the pose
    
    pose = rosetta.core.import_pose.pose_from_file(input_pdb)
    LHLR.simple_pose_moves.remove_insertion_residues(pose, insertion_points)

    if insertion_ids_to_screen is None:
        insertion_ids_to_screen = range(len(insertion_points))
        
    for insertion_id in insertion_ids_to_screen:
        
        # Load linker databases
        
        front_linker_dbs = []
        back_linker_dbs = []

        for f in os.listdir(linker_database_path):
            if f.startswith('selected_linkers_{0}'.format(insertion_id)) and f.endswith('_front.json'):
                with open(os.path.join(linker_database_path, f), 'r') as f:
                    front_linker_dbs.append(json.load(f))

            elif f.startswith('selected_linkers_{0}'.format(insertion_id)) and f.endswith('_back.json'):
                with open(os.path.join(linker_database_path, f), 'r') as f:
                    back_linker_dbs.append(json.load(f))

        selected_lhl_units = LHLR.build_loop_helix_loop_unit.screen_all_loop_helix_loop_units(data_path, pose, insertion_points[insertion_id]['start'], 
                insertion_points[insertion_id]['stop'], front_linker_dbs, back_linker_dbs, num_jobs, job_id,
                max_num_success_each_db_pair=max_num_success_each_db_pair, num_res_clashes_tolerance=num_res_clashes_tolerance)

        # Dump the selected LHL units 

        with open(os.path.join(data_path, 'selected_lhl_units_{0}_{1}.json'.format(insertion_id, job_id)), 'w') as f:
            json.dump(selected_lhl_units, f)

if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1
   
    pyrosetta.init(options='-mute all')
    input_pdb = 'data/1smg_cleaned.pdb'
    linker_database_path = 'data/select_linkers'
    input_insertion_points_file = 'data/insertion_points.json'

    screen(data_path, linker_database_path, input_pdb, input_insertion_points_file, 
            num_jobs, job_id, max_num_success_each_db_pair=None, insertion_ids_to_screen=None, num_res_clashes_tolerance=0)
