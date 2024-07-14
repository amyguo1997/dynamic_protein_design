import subprocess
import os
import sys

if __name__ == '__main__':

    pdb_id = sys.argv[1]

    designs = [f for f in os.listdir(os.path.join('./inputs', pdb_id)) if not f.startswith('.')]
    design_id = designs[int(os.environ['SGE_TASK_ID']) - 1].split('.')[0]

    fasta_file = os.path.join('./inputs', pdb_id, '%s.fasta'%design_id)
    output_dir = os.path.join('./outputs', pdb_id, design_id)
    print('Running colabfold for %s'%design_id)
    os.makedirs(output_dir, exist_ok=True)

    # Check if results already exist
    if not os.path.exists(os.path.join(output_dir, '_%s_plddt.png'%design_id)):
        subprocess.check_call('colabfold_batch --num-recycle 5 --num-models 5 --msa-mode single_sequence --rank plddt %s %s'%(fasta_file, output_dir), shell=True)

