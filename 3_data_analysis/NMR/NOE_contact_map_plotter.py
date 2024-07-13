#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read and plot long range NOE assignments
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def get_structure_contacts(PDB):

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
    
    i_ = []
    j_ = []
    distances = []

    for i,i_xyz in coordinates.items():
        for j,j_xyz in coordinates.items():
            if i > j: 
                dist = np.linalg.norm(i_xyz-j_xyz)
                if dist < 10.0 and abs(i - j) > 4:
                    i_.append(i)
                    j_.append(j)
                    distances.append(dist)
    
    return i_, j_, distances

def get_structure_differences(PDB_1, PDB_2,offset=0):

    with open(PDB_1, 'r') as rf:
        lines_1 = rf.readlines()
    
    coordinates_1 = {}
    for line in lines_1:
        if ' CB' in line:
            residue = int(line.split()[5])
            x = float(line.split()[6])
            y = float(line.split()[7])
            z = float(line.split()[8])
            coordinates_1[residue] = np.array((x,y,z))
            
    with open(PDB_2, 'r') as rf:
        lines_2 = rf.readlines()
    
    coordinates_2 = {}
    for line in lines_2:
        if ' CB' in line:
            residue = int(line.split()[5])
            x = float(line.split()[6])
            y = float(line.split()[7])
            z = float(line.split()[8])
            coordinates_2[residue] = np.array((x,y,z))
    
    i_ = []
    j_ = []
    distances = []

    for i,i_xyz in coordinates_1.items():
        for j,j_xyz in coordinates_1.items():
            if i>j: 
                try:
                    dist_1 = np.linalg.norm(i_xyz-j_xyz)
                    dist_2 = np.linalg.norm(coordinates_2[i]-coordinates_2[j])
                    if (dist_1 < 10.0 or dist_2 < 10.0) and abs(i - j) > 4:
                        
                        difference = (dist_1-dist_2)
    
                        i_.append(i+offset)
                        j_.append(j+offset)
                        
                        distances.append(difference)
                except:continue
             
    return i_, j_, distances

def get_upl(upl_file):
    
    i = []
    j = []
    upl = []
    
    with open(upl_file, 'r') as rf:
        lines = rf.readlines()
    for line in lines:
        line = line.split()
        i_ = int(line[0])
        j_ = int(line[3])
        upl_ = float(line[6])
        

        if i_ < j_: 
            i.append(i_)
            j.append(j_)
            upl.append(upl_)
         
        elif i_ > j_: 
            i.append(j_)
            j.append(i_)
            upl.append(upl_)
                  
    return i,j,upl
        
if __name__ == '__main__':
    
    # Get upper limit distance restraints and plot
    project = '$CYANA_OUTPUT_PATH'
    upl_file = os.path.join(project,'structure.upl')
    i,j,upl,noe_array,noe_i,noe_j,num = get_upl(upl_file,length=94)
    plt.scatter(i,j,color='lightgrey', s=60,edgecolors='black')
    
    # Get contact map for each state and plot difference
    state_1_PDB = './I85_coord.pdb'
    state_2_PDB = './S85_coord.pdb'
    
    i_1, j_1, dist_1 = get_structure_contacts(state_1_PDB)
    i_2, j_2, dist_2 = get_structure_contacts(state_2_PDB)

    i, j, diff = get_structure_differences(state_1_PDB,state_2_PDB)

    plt.scatter(i, j, c=diff, s=60, cmap='coolwarm_r',vmin=-10,vmax=10, edgecolor='grey')
    
   