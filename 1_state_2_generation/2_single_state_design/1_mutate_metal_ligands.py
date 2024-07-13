# -*- coding: utf-8 -*-
"""
Mutate residues in reshaped protein to metal binding ligands
"""

import os
import json
import subprocess

from pyrosetta import *
from pyrosetta import rosetta

init() 

def mutate(pose):
    '''Mutate residues directly.'''
    xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <MOVERS>
        <MutateResidue name="mutate_1" target="66" new_res="ASP" />
        <MutateResidue name="mutate_2" target="68" new_res="ASP" />
        <MutateResidue name="mutate_3" target="69" new_res="GLY" />
        <MutateResidue name="mutate_4" target="70" new_res="SER" />
        <MutateResidue name="mutate_5" target="71" new_res="GLY" />
        <MutateResidue name="mutate_6" target="72" new_res="THR" />
        <MutateResidue name="mutate_7" target="77" new_res="GLU" />
    </MOVERS>
    ''')
    mutate_1 = xmlobj.get_mover('mutate_1')
    mutate_2 = xmlobj.get_mover('mutate_2')
    mutate_3 = xmlobj.get_mover('mutate_3')
    mutate_4 = xmlobj.get_mover('mutate_4')
    mutate_5 = xmlobj.get_mover('mutate_5')
    mutate_6 = xmlobj.get_mover('mutate_6')
    mutate_7 = xmlobj.get_mover('mutate_7')

    mutate_1.apply(pose)
    mutate_2.apply(pose) 
    mutate_3.apply(pose)
    mutate_4.apply(pose)
    mutate_5.apply(pose) 
    mutate_6.apply(pose)
    mutate_7.apply(pose)
    

if __name__ == '__main__':
    
    relax_path = os.path.expanduser('~/rosetta_src_2019.35.60890_bundle/main/source/bin/relax.default.linuxgccrelease')
    db_path = os.path.expanduser('~/rosetta_src_2019.35.60890_bundle/main/database')
    lucs_models_path = os.path.expanduser('~/lucs_models')
    
    # Get the model number for parallelization on CPUs
    model_rmsd_path = os.path.expanduser('~/data/rmsd.json') # File with RMSD values
    output_path = os.path.expanduser('~/data/outputs/')
    os.makedirs(output_path, exist_ok = True)
    
    model_ids = []
    with open(model_rmsd_path,'r') as rf:
        model_rmsd_data = json.load(rf)
        
    for model_id, rmsd in model_rmsd_data.items():
        if rmsd < 10 and rmsd > 3:
            model_ids.append(model_id)
    
    task_id = int(os.environ['SGE_TASK_ID'])
    model_id = model_ids[task_id - 1]
    
    # Make and navigate to design folder
    model_path = os.path.join(output_path, model_id)  
    os.mkdir(model_path, exist_ok = True)
    os.chdir(model_path)
    
    cwd = os.getcwd()
    
    # Mutate metal ligands
    pdb_path = os.path.join(lucs_models_path, 'model_%s.pdb.gz'%model_id)
    pose = pose_from_pdb(pdb_path)
    mutate(pose)
    pose.dump_pdb('model_mut_%s.pdb'%model_id)

    # Relax
    relax_call_cmds = [relax_path, '-database', db_path, 
                       '-relax:constrain_relax_to_start_coords',
                       '-ex1','-ex2','-use_input_sc','-flip_HNQ',
                       '-no_optH false','-relax:ramp_constraints',
                       'false', '-out:path:pdb','%s','-s', 'model_mut_%s.pdb']
    
    relax_call = ' '.join(relax_call_cmds)%(model_path,model_id)
    subprocess.call(relax_call, shell=True)
    os.remove('./model_mut_%s.pdb'%model_id)
