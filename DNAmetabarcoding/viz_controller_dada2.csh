#!/bin/tcsh

# activate the conda environment
conda activate /usr/local/usrapps/trnL_blast/jrabasc/meta_and_qiime

# set the input list folder
set file_dir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S
set run_name = 16S_viz
set outdir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2

set fasta_loc = ${file_dir}/16S_rep_seqs.fna
set output_loc = ${file_dir}/qiime_outputs/16S_rep_seqs.qza

./assets/qiime_vizualizer_pt1.py -i ${file_dir}
qiime tools import --input-path ${fasta_loc} --output-path ${output_loc} --type 'FeatureData[Sequence]'

# deactivate the conda environment
#conda deactivate


