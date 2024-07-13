#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate change in local environment
"""

import numpy as np

def get_coordinates(PDB):

    with open(PDB, 'r') as rf:
        lines = rf.readlines()
    
    coordinates = {}
    for line in lines:
        
        if ' CB' in line:
            
            residue = int(line.split()[5])
            x = float(line.split()[6])
            y = float(line.split()[7])
            z = float(line.split()[8])
            coordinates[residue] = np.array((x,y,z))
        if ' CA' in line and 'GLY' in line:
            
            residue = int(line.split()[5])
            x = float(line.split()[6])
            y = float(line.split()[7])
            z = float(line.split()[8])
            coordinates[residue] = np.array((x,y,z))
    
    return coordinates

def get_median_RMSD(state_1_coords, state_2_coords, dist_cutoff=10.0, long_range_cutoff = 4):
    
    residues = state_1_coords.keys()
    
    median_rmsd = {}
    
    for i in residues:
        rmsds = []
        
        for j in residues: 

            i_xyz_1 = state_1_coords[i]
            j_xyz_1 = state_1_coords[j]
            dist_1 = np.linalg.norm(i_xyz_1-j_xyz_1)
            
            i_xyz_2 = state_2_coords[i]
            j_xyz_2 = state_2_coords[j]
            dist_2 = np.linalg.norm(i_xyz_2-j_xyz_2)
            
            if (dist_1 < dist_cutoff or dist_2 < dist_cutoff) and abs(i - j) > long_range_cutoff:
                rmsds.append(abs(dist_1-dist_2))
        
        if len(rmsds) == 0:
            median_rmsd[i] = 0
        else:
            median_rmsd[i] = np.median(rmsds)

    return median_rmsd
    
if __name__ == "__main__": 
    
    state_1_PDB = './I85_coord.pdb'
    state_2_PDB = './S85_coord.pdb'
    
    state_1_coords = get_coordinates(state_1_PDB)
    state_2_coords = get_coordinates(state_2_PDB)

    
    average_rmsd = get_median_RMSD(state_1_coords, state_2_coords)
    