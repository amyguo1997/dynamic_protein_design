#!/bin/bash
#$ -S /bin/bash
#$ -o job_outputs
#$ -e job_outputs
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=32G
#$ -l scratch=1G
#$ -l h_rt=24:00:00

source ~/mpnn/bin/activate

folder_with_pdbs="../PDB_homooligomers/inputs/af2_6306_min_apo_1smg"

output_dir="../PDB_homooligomers/outputs/1smg_WT_6306/af2_6306_min_apo_1smg_1.0_1.0_noCHX_bias_62_key_only"

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_positions.jsonl"
path_for_bias=$output_dir"/bias_by_res.jsonl"

chains_to_design="A B"
pos_neg_chain_list="A,B"
chain_betas="1.0,1.0"

fixed_positions="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45 46 47 48 49 50 52 53 55 56 59 61 63 65 66 68 69 70 71 72 73 74 75 76 77 78 79 80 82 83 84 86 87 88 89 90, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45 46 47 48 49 50 52 53 55 56 59 61 63 65 66 68 69 70 71 72 73 74 75 76 77 78 79 80 82 83 84 86 87 88 89 90"

python ../helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python ../helper_scripts/make_pos_neg_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --homooligomer 1 --pos_neg_chain_list $pos_neg_chain_list --pos_neg_chain_betas $chain_betas

python ../helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python ../helper_scripts/make_bias_per_res_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_bias

python ../protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --tied_positions_jsonl $path_for_tied_positions \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --bias_by_res_jsonl $path_for_bias \
        --num_seq_per_target 10000 \
        --sampling_temp "0.3" \
        --batch_size 1 \
        --omit_AAs 'CHX' \
        --save_probs 1