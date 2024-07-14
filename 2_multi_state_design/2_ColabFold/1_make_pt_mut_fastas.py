#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make fasta files with point mutations
"""
import os

def make_reversion_fastas(seq_state_1, seq_state_2, output_dir):
    
    '''
    Outputs fasta files containing all point mutants from state_2 aa to state_1 aa
    '''
    os.makedirs(output_dir, exist_ok=True)
    for idx, aa_state_2 in enumerate(list(seq_state_2)):
        
        original_state_2 = list(seq_state_2).copy()
        aa_state_1 = list(seq_state_1)[idx]
        
        if aa_state_2 != aa_state_1:
            original_state_2[idx] = aa_state_1
            fname = '%s%d%s'%(aa_state_2, (idx+1), aa_state_1)
            f = open(os.path.join(output_dir, '%s.fasta'%fname),'w')
            f.write('> %s\n'%fname)
            f.write(''.join(original_state_2))
            f.close()

    
def make_DMS_fastas(seq, idx_list, output_dir):
    '''
    Make fastas for all amino acids at a list of positions given by idx_list (0-indexed)
    '''
    os.makedirs(output_dir, exist_ok=True)
    aas = ['A','D','E','F','I','K','L','M','N','Q','R','S','T','V','W','Y']
    seq = list(seq)
    for idx in idx_list:
        tmp = seq.copy()
        for aa in aas:
            tmp[idx] = aa
        
            f = open(os.path.join(output_dir, '%s%d%s.fasta'%(seq[idx],(idx+1),aa)),'w')
            f.write('> %s%d%s\n'%(seq[idx],(idx+1),aa))
            f.write(''.join(tmp))
            f.close()
        
if __name__ == '__main__':
    
    seq_state_1 = 'ASMTDQQAEARAFLSEEMIAEFKAAFDMFDADGGGDISTKAFGTVMRMLGQNPDKEEQDAIIEEVDEDGSGTIDFEEFLVMMVRQMKEDA'
    seq_state_2 = 'ASMEDLQAEARAFLSEEMIAEFKAAFDMFDADGGGDISYKAVGTVFRMLGINPSKEVLDYLKEKIDVDGSGTIDFEEFLVLMVYIMKQDA'
    
    fasta_dir = '$output_path'
   
