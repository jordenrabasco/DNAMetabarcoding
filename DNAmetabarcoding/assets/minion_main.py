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
import numpy as np
from Bio import pairwise2
import json


def main(match_point, mismatch_point, start_gap_point, ext_gap_point, rand_multi, iteration):
    print("Processing Parameter Iteration: " + str(iteration))


    primer_dir = '/gpfs_common/share02/test2/jrabasc/Workflow/Primers/'

    primer_data_dict = {}
    for filename in os.listdir(primer_dir):
        if (".fasta" in filename):
            f = os.path.join(primer_dir, filename)
            # print name of file being worked on
            for (record_one) in zip(SeqIO.parse(f, "fasta")):
                primer_id = record_one[0].id
                temp_arr = str(primer_id).split('_')
                primer_key = temp_arr[0]
                primer_seq = record_one[0].seq
                if primer_key not in primer_data_dict.keys():
                    primer_data_dict[primer_key] = [primer_seq]
                else:
                    primer_data_dict[primer_key].append(primer_seq)



    univ_adapt_illumina ='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    directory = '/gpfs_common/share02/test2/jrabasc/Workflow/Samples/minion_data'
    #expected_dict = json.load(open("C://Users//jorde//OneDrive//Desktop//test_file_json.txt"))
    assessment_dic = {}

    test_index = 0
    rand_base = ['A', 'T', 'G', 'C']
    adapter_ided = True

    for filename in os.listdir(directory):
        if("test_file.fastq.gz" in filename):
            f = os.path.join(directory, filename)
            with gzip.open(f, "rt") as handle_r_one:
                for (record_one) in zip(SeqIO.parse(handle_r_one, "fastq")):
                    test_index = test_index + 1
                    fwd_seq = record_one[0].seq
                    key = record_one[0].id
                    total_score = 0
                    read_desig = ''
                    new_key = ''
                    rand_seq = ''
                    for keyer in primer_data_dict.keys():
                        primer_seqs = primer_data_dict[keyer]
                        fwd_primer = primer_seqs[0]
                        rev_primer = primer_seqs[1]
                        fwd_align_score = pairwise2.align.globalms(fwd_seq, fwd_primer, match_point, mismatch_point, start_gap_point, ext_gap_point,
                                                               penalize_end_gaps=False, score_only=True)
                        rev_align_score = pairwise2.align.globalms(fwd_seq, rev_primer, match_point, mismatch_point, start_gap_point, ext_gap_point,
                                                               penalize_end_gaps=False, score_only=True)
                        fwd_comp_fwd_align_score = pairwise2.align.globalms(fwd_seq, fwd_primer.reverse_complement(), match_point, mismatch_point, start_gap_point, ext_gap_point,
                                                               penalize_end_gaps=False, score_only=True)
                        rev_comp_rev_align_score = pairwise2.align.globalms(fwd_seq, rev_primer.reverse_complement(), match_point, mismatch_point, start_gap_point, ext_gap_point,
                                                               penalize_end_gaps=False, score_only=True)
                        read_total_score = fwd_align_score + rev_align_score + fwd_comp_fwd_align_score + rev_comp_rev_align_score
                        if read_total_score > total_score:
                            total_score = read_total_score
                            read_desig = keyer
                    if len(rand_seq) < len(fwd_primer):
                        new_base = rand_base[random.randint(0, 3)]
                        rand_seq = rand_seq + new_base
                    rand_align_score = pairwise2.align.globalms(fwd_seq, rand_seq, 3, -1.5, -2, -1.9,
                                                            penalize_end_gaps=False, score_only=True)*0.5
                    if((rand_align_score*rand_multi) > read_total_score):
                        read_desig = 'no_primer'
                        if(('full' not in key) and ('Partial' not in key)):
                            new_key = read_desig + '_correct_id'
                        else:
                            new_key = read_desig + '_' + key + '_mismatch'
                    else:
                        if read_desig in key:
                            new_key = read_desig + '_correct_id'
                        else:
                            new_key = read_desig + '_' + key + '_mismatch'

                    if new_key in assessment_dic.keys():
                        assessment_dic[new_key] = assessment_dic[new_key] + 1
                    else:
                        assessment_dic[new_key] = 1

                #print(assessment_dic.keys())
                #print(assessment_dic)
                num_correct = 0
                num_incorrect = 0
                for key in assessment_dic.keys():
                    if 'correct' in key:
                        num_correct = num_correct + assessment_dic[key]
                    else:
                        num_incorrect = num_incorrect + assessment_dic[key]
                #print(assessment_dic['no_primer_correct_id'])
                #print("Number Correctly ided: " + str(num_correct))
                #print("Number Incorrectly ided: " + str(num_incorrect))
                return num_correct
test_match_point = list(np.linspace(0.1, 5, 9))
test_mismatch_point = list(np.linspace(-3.0, 1.0, 9))
test_start_gap_point = list(np.linspace(-5, -0.5, 9))
test_ext_gap_point = list(np.linspace(-5, -0.5, 9))
test_rand_multiplyer = list(np.linspace(1.0, 6.0, 9))
best_params_indic = 0
best_params = []
iteration = 0
for rand_multi in test_rand_multiplyer:
    for match_point in test_match_point:
        for mismatch_point in test_mismatch_point:
            for start_gap_point in test_start_gap_point:
                for ext_gap_point in test_ext_gap_point:
                    iteration = iteration + 1
                    print('Processing Parameter Iteration: ' + str(iteration) + ' of' + str((len(test_match_point)*len(test_mismatch_point)*len(test_start_gap_point)*len(test_ext_gap_point)*len(test_rand_multiplyer))))
                    if ext_gap_point > start_gap_point:
                        temp_best_params_indic = main(match_point, mismatch_point, start_gap_point, ext_gap_point, rand_multi, iteration)
                        if temp_best_params_indic > best_params_indic:
                            print('Best Params Updated ' + str(best_params))
                            print('% ided' + str(((best_params_indic/4000)*100)))
                            best_params_indic = temp_best_params_indic
                            best_params = [match_point, mismatch_point, start_gap_point, ext_gap_point, rand_multi]
print(best_params_indic)
print(best_params)



if __name__ == '__main__':
    main()
