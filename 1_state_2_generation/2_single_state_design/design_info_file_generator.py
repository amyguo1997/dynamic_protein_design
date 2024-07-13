# -*- coding: utf-8 -*-
"""
local_protein_sequence_design
"""

import os
import json

import pyrosetta
from pyrosetta import rosetta

def find_surrounding_seqposes_noGP(pose, central_residue_ids, cutoff_distance=10, pre_moved_bb_pose=None):
    '''Return the residue ids that surround a given list of
    central residues and the pose of unmoved backbone.
    The selected residues are within a cutoff CA distance and
    the CA-CB vectors are pointing to the central residues.
    GLYs and PROs are ignored
    '''
    rest_of_residues = [i for i in range(1, pose.size() + 1) if not i in central_residue_ids]
    surrounding_residues = set()
    
    central_residues = [pose.residue(i) for i in central_residue_ids]
    if not (pre_moved_bb_pose is None):
        for i in range(1, pre_moved_bb_pose.size() + 1):
            central_residues.append(pre_moved_bb_pose.residue(i))

    for res1_id in rest_of_residues:
        if pose.residue(res1_id).name3() in ['GLY', 'PRO']: 
            continue 
        ca1 = pose.residue(res1_id).xyz('CA')
        cb1 = pose.residue(res1_id).xyz('CB')

        for res2 in central_residues:
            nbra2 = res2.nbr_atom_xyz()

            if ca1.distance(nbra2) > cutoff_distance:
                continue

            if (cb1 - ca1).normalize().dot((nbra2 - ca1).normalize()) > 0.5:
                surrounding_residues.add(res1_id)

    return list(surrounding_residues)

def select_designable_residues(pose, bb_remodeled_residues, ignore_GP=True, pre_moved_bb_pose=None):
    '''Select residues that should be set to designable.
    Return a list of residue ids. 
    '''
    raw_designable_residues = bb_remodeled_residues + find_surrounding_seqposes_noGP(pose, bb_remodeled_residues, cutoff_distance=10, pre_moved_bb_pose=pre_moved_bb_pose)
    
    designable_residues = []
    
    for r in raw_designable_residues:
        if ignore_GP and pose.residue(r).name3() in ['GLY', 'PRO']:
            continue

        designable_residues.append(r)

    return designable_residues

def find_burial_state(pose, layer=None):

    if layer == 'core':
        xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        '''
        <RESIDUE_SELECTORS>
            <Layer name="core" select_core="true" ball_radius="2.0" />  
        </RESIDUE_SELECTORS>
        ''')
        core = xmlobj.get_residue_selector('core')
        core_residues = core.apply(pose)
        print(core_residues)
        core_list = []
        for idx, bool_ in enumerate(core_residues):
            if bool_:
                core_list.append(idx+1)
        return core_list
    
    elif layer == 'surface':
        xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        '''
        <RESIDUE_SELECTORS>
            <Layer name="surface" select_surface="true" ball_radius="2.0" />  
        </RESIDUE_SELECTORS>
        ''')
        surface = xmlobj.get_residue_selector('surface')
        surface_residues = surface.apply(pose)
        surface_list = []
        for idx, bool_ in enumerate(surface_residues):
            if bool_:
                surface_list.append(idx+1)
        return surface_list
    
def make_one_design(output_path, input_pdb, bb_remodeled_residues, designable_residues=None, repackable_residues=None, core_residues=None, surface_residues=None,
        pre_moved_bb_pose=None, residues_to_fix=[]):
    '''Make one design and dump the relative information.
    Args:
        output_path: path for the outputs
        input_pdb: path to the input pdb file
        bb_remodeled_residues: a list of sequence positions for backbone
            remodeled residues
    '''
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)
    sequence = pose.sequence()
    
    length = int(pose.total_residue())
    
    # Find designable and repackable residues

    if (designable_residues is None):
        designable_residues_tmp = select_designable_residues(pose, bb_remodeled_residues, pre_moved_bb_pose=pre_moved_bb_pose)
        designable_residues = [r for r in designable_residues_tmp if not (r in residues_to_fix)]
        
    if repackable_residues is None:    
        repackable_residues = find_surrounding_seqposes_noGP(pose, designable_residues, cutoff_distance=8)
    
    # Find burial state of residues
    
    if core_residues is None:
        core_residues = find_burial_state(pose, layer='core')
    
    if surface_residues is None:
        surface_residues = find_burial_state(pose, layer='surface')
        
    # Dump information

    info_dict = {
            'bb_remodeled_residues' : bb_remodeled_residues,
            'designable_residues' : designable_residues,
            'repackable_residues' : repackable_residues,
            'length' : length,
            'sequence' : sequence,
            'core_residues' : core_residues,
            'surface_residues' : surface_residues
            }

    with open(os.path.join(output_path, 'design_info.json'), 'w') as f:
        json.dump(info_dict, f)
