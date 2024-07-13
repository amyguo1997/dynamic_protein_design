import os
import json
import numpy as np

from pyrosetta import *
from pyrosetta import rosetta

init()

from biased_ff import IO 
from biased_ff import ab_initio_frag_func_lib 

# Fill in paths!
ss = {
      "runpsipred_single" : "$PSIPRED/runpsipred_single",
      "fragment_picker" : "$ROSETTA/main/source/bin/fragment_picker.static.linuxgccrelease",
      "vall" : "$ROSETTA/tools/fragment_tools/vall.jul19.2011.gz",
      "csblast" : "$CSBLAST",
      "blastpgp" : "$BLAST/bin/blastpgp",
      "placeholder_seqs" : "~/frag_qual_analysis_files/placeholder/placeholder_seqs",
      "sparksx_path" : "$SPARKSX/sparks-x",
      "rosetta_database_fragment_picking" : "$ROSETTA/main/database",
      "dalphaball" : "$ROSETTA/main/source/external/DAlpahBall/DAlphaBall.gcc",
      "fragment_quality_analysis_weights" : "~/frag_qual_analysis_files/standard.wghts"
      }

def get_fragment_quality_scores(pose, design_path, design_num, model_id):
    ''' Get the frag quality score in the remodeled region.
    Returns the fragment with the highest rmsd and the avg rmsd.'''
    
    # Create FASTA file
    sequence = pose.sequence()
    IO.sequence_to_fasta_file('design_{0}.fasta'.format(design_num),'design',sequence)
    
    # Fragment picking    
    fqa = ab_initio_frag_func_lib.FragmentQualityAnalyzer(
           ss['runpsipred_single'], ss['csblast'], ss['blastpgp'], ss['placeholder_seqs'], ss['sparksx_path'],
           ss['fragment_picker'], ss['vall'], ss['fragment_quality_analysis_weights'],
           rosetta_database=ss['rosetta_database_fragment_picking'])

    fqa.pick_fragments('design_{0}_{1}.pdb.gz'.format(model_id,design_num), 'design_{0}.fasta'.format(design_num), design_path)    