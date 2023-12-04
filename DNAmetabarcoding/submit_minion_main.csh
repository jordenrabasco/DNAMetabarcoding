#!/bin/tcsh

# activate the conda environment
conda activate /usr/local/usrapps/trnL_blast/jrabasc/meta_and_qiime


# set the input list folder, the primer fasta file, and the output directory
set file_dir = /gpfs_common/share02/test2/jrabasc/Workflow/Samples/minion_data
set outdir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S
set TMP = /gpfs_common/share02/test2/jrabasc/tmp/dnametabarcoding/seq_2/16S

#set runtime options for script
set cutoff = 200
set threads = 8

#set taxanomiic assingnment parameters
set taxmethod = BLAST
set blastdatabase = /usr/local/usrapps/trnL_blast/test_db/nt
set dada2database = UNITE_eukaryote

#set unique run id 
set unique_run_id = test

# create the output directory
mkdir -p $outdir

# submit a job to the cluster for each input file
# the job name and outputs are named based on the filename (with .fastq.gz removed)
# the app.py command can be modified to run either BLAST or DADA2 taxonomy classification
foreach file ("$file_dir"/*)
  set filename = `basename $file .fastq.gz`
  bsub \
      -J $filename \
      -W 500 \
      -n $threads \
      -x \
      -R "span[hosts=1]" \
      -o $outdir/$filename.out.%J \
      -e $outdir/$filename.err.%J \
      "./minion_main.py -i ${file} -o ${outdir}/${filename}.csv --temp_folder ${TMP} --taxmethod ${taxmethod} --blastdatabase ${blastdatabase} --taxreference ${dada2database} --threads ${threads} --cutoff ${cutoff}"
end

# deactivate the conda environment
conda deactivate


