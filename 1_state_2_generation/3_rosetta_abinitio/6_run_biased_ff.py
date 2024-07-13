import subprocess
import os
import json

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

with open('design_info.json', 'r') as rf:
    data = json.load(rf)
best_design_num = data['best_design_num']

call_cmds = ["$ROSETTA/main/source/bin/minirosetta.static.linuxgccrelease",
             "-ex1 1 -ex2 1", "-abinitio::fastrelax 1", "-relax::default_repeats 15",
             "-abinitio::use_filters false", "-abinitio::increase_cycles 10", 
             "-frag3 frags.3t200.3mers", "-frag9 frags.3t200.9mers", 
             "-abinitio::number_3mer_frags 1", "-abinitio::number_9mer_frags 1", 
             "-abinitio::rsd_wt_loop 0.5", "-abinitio::rsd_wt_helix 0.5", 
             "-abinitio::rg_reweight 0.5", "-in:file:native design_%s_%s.pdb.gz", 
             "-silent_gz 1", "-out:file:silent", "fold_design.out", "-out:file:scorefile fold_design.sc",
             "-out:path %s -nstruct 30"]

call = ' '.join(call_cmds)                                                                                                                                                                                                                
subprocess.call(call%(model_id, best_design_num, model_path), shell = True)