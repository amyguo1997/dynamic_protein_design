# -*- coding: utf-8 -*-
"""
Single-state state 2 design with Rosetta
"""

import os
import json

import pyrosetta
from pyrosetta import rosetta


def get_task_factory(pose, designable_residues, repackable_residues, extra_rotamers=True, limit_aro_chi2=True, layered_design=True,
        designable_aa_types=None):
    '''Get a task factory given the designable and repackable residues.'''
    def list_to_str(l):
        return ','.join(list(str(i) for i in l))

    task_factory = rosetta.core.pack.task.TaskFactory()

    if len(designable_residues) > 0:
        for i in range(len(designable_residues)):
            designable_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str([designable_residues[i]])) 
            racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()

            if designable_aa_types is None or len(designable_residues) != len(designable_aa_types):
                racaa.aas_to_keep('GAPVILMFYWSTKRDENQ') # No CYS or HIS
            else:
                racaa.aas_to_keep(designable_aa_types[i])

            designable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                    racaa, designable_selector)
            task_factory.push_back(designable_operation)

    if len(repackable_residues) > 0:
        repackable_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(repackable_residues)) 
        repackable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.RestrictToRepackingRLT(), repackable_selector)
        task_factory.push_back(repackable_operation)

    natro_residues = [i for i in range(1, pose.size() + 1) if not i in designable_residues + repackable_residues]
    if len(natro_residues) > 0:
        natro_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(natro_residues)) 
        natro_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.PreventRepackingRLT(), natro_selector)
        task_factory.push_back(natro_operation)

    if extra_rotamers:
        ers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        ers.ex1(True)
        ers.ex2(True)
        ers.extrachi_cutoff(18)

        task_factory.push_back(ers)

    if limit_aro_chi2:
        lac = rosetta.protocols.task_operations.LimitAromaChi2Operation()

        task_factory.push_back(lac)

    if layered_design:
        ld = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_task_operation(
            '''<LayerDesign name="layer_all" layer="core_boundary_surface_Nterm_Cterm" use_sidechain_neighbors="True">
    		<Nterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Nterm>
    		<Cterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Cterm>
        </LayerDesign>''')
        task_factory.push_back(ld)

    return task_factory

def get_move_map(bb_movable_residues, sc_movable_residues):
    '''Get a move map given the backbone movable residues, side chain
    movable residues and movable jumps. 
    Note that the side chains of bb_movable_residues are also movable.
    '''
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(False)

    for i in bb_movable_residues:
        mm.set_bb(i, True)
        mm.set_chi(i, True)

    for i in sc_movable_residues:
        mm.set_chi(i, True)
    
    return mm

def fast_design(pose, designable_residues, repackable_residues, fast_design_repeats=1, move_map=None, extra_rotamers=False):
    '''Do fast design.'''
    xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <MOVERS>
        <FastDesign name="fastdes" clear_designable_residues="0" repeats="{0}" ramp_down_constraints="1"/>
    </MOVERS>
    '''.format(fast_design_repeats))
    fast_design = xmlobj.get_mover('fastdes')

    # Use the rosettcon2018 relax script for better design

    script = rosetta.std.vector_std_string()
    script.append('repeat {0}'.format(fast_design_repeats))
    script.append('ramp_repack_min 0.079 0.01     1.0')
    script.append('ramp_repack_min 0.295 0.01     0.5')
    script.append('ramp_repack_min 0.577 0.01     0.0')
    script.append('ramp_repack_min 1     0.00001  0.0')
    script.append('accept_to_best')
    script.append('endrepeat')

    fast_design.set_script_from_lines(script)

    # Set score function

    sfxn = rosetta.core.scoring.get_score_function() 
    fast_design.set_scorefxn(sfxn)

    # Set task factory

    task_factory = get_task_factory(pose, designable_residues, repackable_residues, extra_rotamers=extra_rotamers)
    fast_design.set_task_factory(task_factory)

    # Set move_map
   
    if not move_map is None:
        fast_design.set_movemap(move_map)

    fast_design.apply(pose)

if __name__ == '__main__':
    
    pyrosetta.init()
    sfxn = rosetta.core.scoring.get_score_function()
    fastrelax = pyrosetta.rosetta.protocols.relax.FastRelax(sfxn, 1) 
    
    # Paths
    model_rmsd_path = os.path.expanduser('~/data/rmsd.json') # File with RMSD values
    output_path = os.path.join('~/data/outputs')

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
    design_info_path = os.path.join(model_path, 'design_info.json')
    pdb_path = os.path.join(model_path,'model_mut_%s_0001.pdb'%model_id)
   
    # Define the positions that are designable on the pose
    
    with open(design_info_path, 'r') as rf:
        data = json.load(rf)
    designable_residues = data['designable_residues']
    repackable_residues = data['repackable_residues']

    # Make the pose
    
    pose = pyrosetta.pose_from_file(pdb_path)

    # Define the move map
    
    mm = get_move_map(designable_residues + repackable_residues, [])
    os.makedirs(os.path.join(model_path, 'designs'), exist_ok=True)
      
    for i in range(20):         
        if os.path.exists(os.path.join(model_path,'designs','design_%s_%s.pdb.gz'%(model_id, str(i)))): continue   
        else:
            
            fast_design(pose, designable_residues, repackable_residues, move_map=mm)
            fastrelax.apply(pose)               
                             
            with open(os.path.join(model_path, 'designs', 'design_%s_info.json'%str(i)), 'w') as f:
                json.dump({'score':pose.energies().total_energy(), 'sequence':pose.sequence()}, f)
                        
            pose.dump_pdb(os.path.join(model_path, 'designs', 'design_%s_%s.pdb.gz'%(model_id,str(i))))
       
