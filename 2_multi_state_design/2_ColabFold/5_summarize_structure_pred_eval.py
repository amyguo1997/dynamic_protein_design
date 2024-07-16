#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarirze af structure data output from structure_pred_eval.py and visualize
"""
import os
import json

def summarize(af_dir):
    
    summarized_data = {}
    
    for root, dirs, files in os.walk(af_dir):
        for design_id in dirs:
            with open(os.path.join(af_dir, design_id, 'structure_pred_eval.json'), 'r') as rf:
                data = json.load(rf)
            summarized_data[design_id] = data
        
    
    with open(os.path.join(af_dir, 'summarized_structure_pred_eval.json'), 'w') as wf:
        json.dump(summarized_data, wf)
        
    
if __name__ == '__main__':
    
    af_dir = '$AF2_prediction_path' # contains subdirectories containing AF2 *.pdb predictions for each design
    summarized_data_path = os.path.join(af_dir,'summarized_structure_pred_eval.json')
