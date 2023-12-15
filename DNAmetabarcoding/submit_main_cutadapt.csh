#!/bin/tcsh

# activate the conda environment
conda activate /usr/local/usrapps/trnL_blast/jrabasc/meta_and_qiime

#set the input folder that has the fastq files that will be run through the pipeline
set file_dir = /gpfs_common/share02/test2/jrabasc/Workflow/Samples/expanded_barcodes/seqrun2_test/16S

#set location of primer fasta files
set forward_primers = /gpfs_common/share02/test2/jrabasc/Workflow/Primers/fwd_all.fasta
set reverse_primers = /gpfs_common/share02/test2/jrabasc/Workflow/Primers/rev_all.fasta

#set output (final output) and temporary (intermediate files) directories for all files
set outdir = /gpfs_common/share02/test2/jrabasc/experiment_output/seq_2/16S
set TMP = /gpfs_common/share02/test2/jrabasc/tmp/dnametabarcoding/seq_2/16S

#set cutadapt -m option (removes reads smaller than the cutoff variable)
set cutoff = 50

#set runtime options for script
set threads = 8

#set taxanomiic assingnment parameters
set taxmethod = BLAST
set blastdatabase = /usr/local/usrapps/trnL_blast/test_db/nt
set dada2database = UNITE_eukaryote

#set unique run id 
set unique_run_id = test

#set whether to trim both fwd and rev primers or just one of them. If a "both_" option is specified the algorithm will filter by the primer specified
set trim_both = both_fwd

#set which end of the read the specified primer in the trim_both option is located at
set prim_three_or_five_primer = three

#Select what orientation the primers will be in. False=fwd primer in 3`-5` and rev in 5`-3`, True=fwd primer in 5`-3` and rev in 3`-5`
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


