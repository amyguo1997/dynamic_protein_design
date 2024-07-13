''' 
Iterates through designs and implements filter script
'''
import sys
sys.path.append('..')

import os
import json

from pyrosetta import *
from pyrosetta import rosetta

init() 

from rosetta_abinitio import ab_initio_frag_dependencies

sfxn = get_fa_scorefxn()

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
os.chdir(model_path)

# Find best score 
scores = []
for root, dirs, files in os.walk('.'):
    for fname in files:
        if fname.startswith('design') and fname.endswith('.pdb.gz'):
            pose = pose_from_pdb(fname)
            score = sfxn(pose)
            scores.append([fname.split('.')[0].split('_')[-1],score])

best_design_num = sorted(scores, key = lambda x: x[1])[0][0]
with open('design_info.json','r') as rf:
    data = json.load(rf)
data['best_design_num'] = best_design_num

with open('design_info.json', 'w') as wf:
    json.dump(data,wf)

pose = pose_from_pdb('design_{0}_{1}.pdb.gz'.format(model_id,best_design_num))
ab_initio_frag_dependencies.get_fragment_quality_scores(pose, model_path, best_design_num, model_id)