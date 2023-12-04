#!/bin/tcsh

# activate the conda environment
conda activate /usr/local/usrapps/trnL_blast/jrabasc/meta_and_qiime


# set the input list folder, the primer fasta file, and the output directory
set file_dir = /gpfs_common/share02/test2/jrabasc/Workflow/Samples/expanded_barcodes/seqrun2_test/16S
set forward_primers = /gpfs_common/share02/test2/jrabasc/Workflow/Primers/fwd_all.fasta
set reverse_primers = /gpfs_common/share02/test2/jrabasc/Workflow/Primers/rev_all.fasta
set outdir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S
set TMP = /gpfs_common/share02/test2/jrabasc/tmp/dnametabarcoding/seq_2/16S

#set cutadapt cutoff
set cutoff = 50

#set runtime options for script
set threads = 8

#set taxanomiic assingnment parameters
set taxmethod = BLAST
set blastdatabase = /usr/local/usrapps/trnL_blast/test_db/nt
set dada2database = UNITE_eukaryote

#set unique run id 
set unique_run_id = test

#set whether to trim fwd/rev/both_rev/both_fwd primers, and which ones to primers to filter by  
set trim_both = both_fwd
set prim_three_or_five_primer = three
set reverse_comp = True


# create the output directory
mkdir -p $outdir

# submit a job to the cluster for each input file
# the job name and outputs are named based on the filename (with .fastq.gz removed)
# the app.py command can be modified to run either BLAST or DADA2 taxonomy classification
cd assets

foreach file ("$file_dir"/*)
  set filename = `basename $file .fastq.gz`
  bsub \
      -J $filename \
      -W 500 \
      -n $threads \
      -R "rusage[mem=200GB]" \
      -q shared_memory \
      -o $outdir/$filename.out.%J \
      -e $outdir/$filename.err.%J \
      "./main_cutadapt.py -i ${file} -o ${outdir}/${filename}.csv --temp_folder ${TMP} --forward_primers ${forward_primers} --reverse_primers ${reverse_primers} --taxmethod ${taxmethod} --blastdatabase ${blastdatabase} --taxreference ${dada2database} --threads ${threads} --cutoff ${cutoff} --unique_run_id ${unique_run_id} --trim_both ${trim_both} --prim_three_or_five_primer ${prim_three_or_five_primer} --reverse_comp ${reverse_comp}"
end

cd ..

# deactivate the conda environment
conda deactivate


