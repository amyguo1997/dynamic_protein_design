import subprocess
import os
import json

if __name__ == '__main__':
    
    # Get the model number for parallelization on CPUs
    script_path = 'biased_frag_generator.py'
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
    
    # Navigate to the design folder

    model_path = os.path.join(output_path, model_id)  
    os.chdir(model_path)
    
    call_3mers = 'python %s -frag_qual frags.fsc.200.3mers -ntop 3 -fullmer frags.200.3mers -out frags.3t200.3mers'
    call_9mers = 'python %s -frag_qual frags.fsc.200.9mers -ntop 3 -fullmer frags.200.9mers -out frags.3t200.9mers'
    
    subprocess.call(call_3mers%script_path, shell=True)
    subprocess.call(call_9mers%script_path, shell=True)