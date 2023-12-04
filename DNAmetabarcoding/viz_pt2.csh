#!/bin/tcsh

conda activate /usr/local/usrapps/trnL_blast/jrabasc/meta_and_qiime

# set the input list folder
set file_dir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S
set run_name = 16S_viz
set outdir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/16S_rep_seqs.qza \
  --o-alignment /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/aligned-rep-seqs.qza \
  --o-masked-alignment /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/masked-aligned-rep-seqs.qza \
  --o-tree /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/unrooted-tree.qza \
  --o-rooted-tree /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/rooted-tree.qza \
  --i-table /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S/qiime_outputs/core-metrics-results

# deactivate the conda environment
#conda deactivate


