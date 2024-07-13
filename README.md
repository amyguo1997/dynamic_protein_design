# Deep learning guided design of dynamic proteins
Code accompanying "Deep learning guided design of dynamic proteins" by Amy B. Guo, Deniz Akpinaroglu, Mark J.S. Kelly, and Tanja Kortemme.

Data, scripts, and molecular dynamics trajectories can be downloaded here: https://doi.org/10.5061/dryad.m37pvmdbm.

![Design approach](https://github.com/user-attachments/assets/3925e8fd-165a-4fa8-8238-b0753d32e522)

## Generating non-native alternative states (single-state design)
Given a starting functional backbone (state 1), we first generate alternative conformations (state 2) using loop-helix-loop unit combinatorial sampling. Detailed documentation on this Python-based method can be found at: https://github.com/Kortemme-Lab/loop_helix_loop_reshaping. Job scripts, analysis scripts, and inputs specific for this study can be found in state_2_generation/loop_helix_loop_reshaping.

To evaluate the designability of each backbone, we then perform single-state design on each generated backbone. Scripts for (1) mutating residue positions in the alternative state essential for forming the functional motif (Ca<sup>2+</sup> binding site), creating a design information file specifying designable/repackable residues, and running Rosetta LayerDesign can be found in state_2_generation/single_state_design. The lowest energy design is then evaluated using Rosetta _ab initio_ structure prediction. Scripts for fragment generation and biased forward folding can be found in state_2_generation/rosetta_abinitio. 

## Deep-learning guided multi-state design
A high-sequence identity design predicted to fold into the state 2 backbone was identified by evaluating _in silico_ mutations increasing sequence identity to state 1 with [ColabFold](https://github.com/sokrypton/ColabFold). Designable residues sampled during [multi-state design](https://github.com/dauparas/ProteinMPNN) were restricted to only positions with differing amino acids between states and their neighbors. The job submission script used in this study is included in multi_state_design/ProteinMPNN_scripts. The resulting designs were evaluated by ColabFold.

## Analysis of dynamic designs
Dynamic designs were characterized using nuclear magnetic resonance (NMR) and molecular dynamics simulations. Scripts for analyzing NOESY-derived distance restraints and approximating the change in local chemical environment at each residue position can be found in data_analysis/NMR. Detailed documentation on mutual information analysis can be found at https://github.com/stefdoerr/mutinf and https://simtk.org/projects/mutinf. 
