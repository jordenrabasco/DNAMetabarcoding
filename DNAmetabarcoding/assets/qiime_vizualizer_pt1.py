#!/usr/bin/env python

import click
import glob
import os
import pandas as pd
import subprocess
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import gzip
import random


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True,
              help='fastq file path for sequence reads')


def main(input):
    #generates dictionary for iding expected taxa and removing all others
    exp_taxa_dict = {'16S': ['superkingdom','Bacteria'], 'ITS1': ['kingdom', 'Fungi'], 'ITS2F': ['kingdom', 'Viridiplantae'], 'ITSp3': ['kingdom', 'Viridiplantae'], 'mLepCOI': ['kingdom', 'Metazoa'], 'ZBJCOI': ['kingdom', 'Metazoa'], 'trnL': ['kingdom', 'Viridiplantae']}
    exp_names_arr = []
    
    #sets up pathing and extracts which files to merge
    pipe_output_dir = input
    
    #matches expected_taxa_dict key with primer indicator in pathing
    looking_for = ''
    for key in exp_taxa_dict.keys():
        for f in os.listdir(pipe_output_dir):
            if key in str(os.path.join(pipe_output_dir, f)):
                looking_for = key
            break
        if looking_for != '':
            break
    
    #extracts expected criteria column and what the expected criteria is
    exp_taxa_id_col = exp_taxa_dict[looking_for][0]
    exp_taxa_ident = exp_taxa_dict[looking_for][1]
    
    #sets up df and output file pathing
    df = pd.DataFrame()
    writer_full_address = pipe_output_dir + '/' + str(looking_for) + '_rep_seqs.fna'
    output_dest = pipe_output_dir + '/' + str(looking_for) + '_exps.tsv'
    
    #iterates through all .csv files in a folder and merges them
    path_index = 0
    for f in os.listdir(pipe_output_dir):
        if(('.csv' in str(f)) and ('untrimmed' not in str(f))):
            exp_name = str(f).split('.csv')[0]
            full_path = os.path.join(pipe_output_dir, f)
            temp_df = pd.read_csv(full_path, sep=',')
            exp_names_arr.append(exp_name)
            del temp_df["stitle"]
            bad_indexes = temp_df.index[temp_df[exp_taxa_id_col] != exp_taxa_ident].tolist()
            temp_df.drop(bad_indexes, inplace=True)
            temp_df = temp_df[['sequence','abundance']]
            temp_df = temp_df.rename(columns={'abundance': exp_name})
            if path_index == 0:
                df = temp_df
            else:
                df = pd.merge(df, temp_df, how='outer', on=['sequence', 'sequence'])
            path_index = path_index + 1
    df = df.fillna(0)
    df = df.set_index('sequence')
    df.index.name = None
    
    #list exps that didn't make it
    bad_exp_names = []
    for file in os.listdir(pipe_output_dir):
        temp_name = str(file).split('_')[0]
        if '.csv' not in temp_name:
            if((temp_name not in exp_names_arr) and (temp_name not in bad_exp_names)):
                bad_exp_names.append(temp_name)
    
    if len(bad_exp_names) == 0:
        print("Every Experiment Made it")
    else:
        for bad_boy in bad_exp_names:
            print("Exp " + str(bad_boy) + ' analysis failed. Please see corresponding error file for more information')
    
    #uses df to generate rep_Sequences to be used in qiime2 analysis
    rep_sequences = df.index
    handle_out = open(writer_full_address, "wt")
    for seq in rep_sequences:
        temp_rec = SeqRecord(id=str(seq), seq=Seq(seq), name='', description='')
        SeqIO.write(temp_rec, handle_out, "fasta")
    
    df.to_csv(output_dest, sep='\t')
    
    qiime_output_dest = pipe_output_dir + '/qiime_outputs'
    os.mkdir(qiime_output_dest)
    
    #file wraggling for qiime input
    #biom_table = f'{qiime_output_dest}' + '/' + str(looking_for) +'_biom-table.txt'
    #cmd_biom = 'echo -n "#OTU Table" | cat - ' + str(output_dest) + ' > ' + str(biom_table)
    #subprocess.check_call(cmd_biom, shell=True)
 
    #os.chdir(qiime_output_dest)

    table_biom_format = f'{qiime_output_dest}' + '/' + str(looking_for) +'_table.biom'    
    cmd_biom_convert = 'biom convert -i ' + str(output_dest) + ' -o ' + str(table_biom_format) + ' --table-type="OTU table" --to-hdf5 '
    subprocess.check_call(cmd_biom_convert, shell=True)
    
    table_qza_address = f'{qiime_output_dest}' + '/' + str(looking_for) + '_feature_table.qza'
    cmd_table_qiime_import = 'qiime tools import --input-path ' + str(table_biom_format) + " --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path " + str(table_qza_address)
    subprocess.check_call(cmd_table_qiime_import, shell=True)
        

if __name__ == '__main__':
    main()
