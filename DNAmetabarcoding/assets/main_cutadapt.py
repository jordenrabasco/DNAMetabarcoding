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

# set blast filtering thresholds
EVALUE = 0.001
IDENTITY = 95
COVERAGE = 90

# set taxonomy ranks to include (in order of increasing specificity)
ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# get username for determining file paths
user = os.environ['USER']


# set dada2 taxonomy database paths to include all fasta files
# in the dada2_taxonomy directory, keyed by the filename
base_dada2_tax = f'/usr/local/usrapps/trnL_blast/dada2_taxonomy'
dada2_tax_dbs = {
    os.path.basename(file).replace('.fasta', ''): file
    for file in glob.glob(f'{base_dada2_tax}/*.fasta')
}

def run_dada2(input, output, feat_table, rep_seqs, TMP):
    """Run DADA2 on the input sequence reads to identify ASVs.

    Arguments:
    input -- the fastq file of sequencing reads
    output -- the file path for the DADA2 results
    feat_table -- the file path for the feature table to be turned into a Qiime
    rep_seqs -- the file oath for the representative sequences from DADA2 output to be turned into a Qiime object
    TMP -- directory for temporary, intermediate files

    Returns:
    a dataframe of ASVs with the columns: index, sequence, and abundance
    """
    subprocess.check_call(['Rscript', 'dada2.R', input, output, feat_table, rep_seqs, TMP])
    if(os.path.isfile(output) == True):
        return pd.read_csv(output, index_col=0).reset_index()
    else:
        return pd.DataFrame()


def dada2_taxonomy(input, output, reference):
    """
    Run DADA2's assignTaxonomy method for taxonomic classification.

    Arguments:
    input -- the file of ASVs output from DADA2
    output -- the file path for the DADA2 taxonomy results
    reference -- the file path for the reference database to use for
      taxonomic assignment

    Returns:
    a dataframe containing the sequence and the taxonomy assignment
    """
    subprocess.check_call([
        'Rscript', 'dada2_taxonomy.R', input, output, reference
    ])
    # read taxonomy results
    taxa = (
        pd.read_csv(output, index_col=0)
        .reset_index()
        .rename(columns={'index': 'sequence'})
        .insert(loc=1, column='superkingdom', value='')
    )
    # change column names to lowercase
    taxa.rename(columns={c: c.lower() for c in taxa.columns}, inplace=True)
    # remove the prefix from the taxon names (ie. 'g__' for genus names)
    for c in ranks:
        taxa[c] = taxa[c].str.split('__').str[1]
    return taxa


def write_blast_fasta(data, output):
    """Create a fasta file from a dataframe for use in BLAST.

    Arguments:
    data -- a dataframe containing sequence and index (identifier) columns
      for the ASVs
    output -- the file path for the fasta file
    """
    records = [
        SeqRecord(Seq(x.sequence), id=str(x.index), description='')
        for x in data.itertuples()
    ]
    SeqIO.write(records, output, 'fasta')


def run_blast(input, output, database, threads):
    """Run BLAST and filter the results.
    The BLAST results are filtered by evalue, percent identity, and coverage.
    The result(s) with maximum percent identity are returned for each input
    sequence.

    Arguments:
    input -- the fasta file to be queried with BLAST
    output -- the file path for the BLAST results
    database -- the BLAST database to be searched
    threads -- the threads (or cpus) to use for BLAST

    Returns:
    a dataframe containing the filtered BLAST results
    """
    subprocess.check_call([
        'blastn',
        '-db', database,
        '-query', input,
        '-max_target_seqs', '10',
        '-evalue', str(EVALUE),
        '-perc_identity', str(IDENTITY),
        '-num_threads', str(threads),
        '-outfmt', '6 qacc sacc qlen slen pident length qcovs staxid stitle sseqid',
        '-out', output
    ])
    # filter blast results
    blast_header = ['qacc', 'sacc', 'qlen', 'slen', 'pident', 'length', 'qcovs', 'staxid', 'stitle', 'sseqid']
    blast_data = pd.read_csv(output, header=None, names=blast_header, quoting=csv.QUOTE_NONE, sep='\t')
    blast_data = blast_data[blast_data['qcovs']>=COVERAGE].reset_index(drop=True)
    blast_data = pd.merge(
    blast_data.groupby('qacc')[['pident']].max().reset_index(),
    blast_data,
    on=['qacc', 'pident'],
    how='inner'
    )
    return blast_data


def run_taxize(input, output):
    """Run taxizedb to look up the taxa for the input taxids.
    taxizedb uses a local database to look up the NCBI taxonomy.

    Arguments:
    input -- the file containing the list of the taxids
    output -- the file path for the taxizedb results

    Returns:
    a dataframe containing the taxids and corresponding taxa at all ranks
    """
    subprocess.check_call(['Rscript', 'taxizedb.R', input, output])
    if(os.path.isfile(output) == True):
        taxa = pd.read_csv(output)
        taxa = (
            taxa[taxa['rank'].isin(ranks)]
            .pivot(index='taxid', columns='rank', values='name')
            .reset_index()
        )
        # remove the genus portion that is included in the species name
        col_check = list(taxa.columns)
        if((col_check.count('species') > 0) and (col_check.count('genus') > 0)):
            taxa['species'] = taxa.apply(
                lambda x:
                str(x.species).replace(str(x.genus), '').strip(),
                axis=1
            )
        #fill null kingdom values from superkingdom
        taxa.fillna('', inplace=True)
        return taxa
    else:
        return pd.DataFrame()

    


def consistent_taxa(x):
    """Determine a consistent taxon at each rank.
    If all results are the same at a given rank, that name is returned for the
    taxon; otherwise, a null value is returned.

    Arguments:
    x -- a dataframe containing the taxa identified for one query sequence

    Returns:
    a dictionary of the consistent taxa names
    """
    new_taxa = {'qacc': x['qacc'].unique().tolist()[0]}
    stop = 0
    for r in ranks:
        if r in x.columns:
            r_taxa = x[r].dropna().unique().tolist()
            if len(r_taxa) == 1:
                new_taxa[r] = r_taxa[0]
            elif len(r_taxa) == 0 and r in ['superkingdom', 'kingdom']:
                new_taxa[r] = None
            else:
                stop += 1
            if stop > 0:
                new_taxa[r] = None
        else:
            new_taxa[r] = None
    return new_taxa


def combine_blast_taxize(blast_data, taxa):
    """Combine the BLAST and taxizedb data.
    For each query sequence, a single consistent set of taxon names are
    returned, or a null value for ranks where the taxon names are not
    consistent.

    Arguments:
    blast_data -- the dataframe containing the BLAST results
    taxa -- the dataframe containing the taxa assignments from taxizedb

    Returns:
    a dataframe containing the resulting taxa assignments
    """
    # merge blast results with taxize information
    merged = (
        pd.merge(
            blast_data[['qacc', 'staxid','stitle']],
            taxa,
            left_on='staxid',
            right_on='taxid',
            how='left'
        )
        .drop(['staxid', 'taxid'], axis=1)
        .drop_duplicates()
        .reset_index(drop=True)
    )
    # keep taxa from BLAST only down to the last consistent rank
    final_taxa = pd.DataFrame.from_dict(
        merged.groupby('qacc')
        .apply(lambda x: consistent_taxa(x))
        .to_list()
    )
    return final_taxa

def parse_primers(input_fw, input_rev, fwd, rev):
    """Prepare the primer sequence for use in trimming.
    A ^ symbol is added to the start of the sequence for use as a forward
    primer, requiring that it is at the start of the read. An X is added
    to the end of the reverse complement sequence for use as a reverse primer,
    requiring that it is at the end of the read and that it can be a partial
    primer sequence.

    Arguments:
    input_fw -- a fasta file containing forward primer sequences
    input_rev -- a fasta file containing reverse primer sequences
    fwd -- output file path for forward primer sequences
    rev -- output file path for reverse complement primer sequences
    """
    fwd_records = []
    rev_records = []
    for record in SeqIO.parse(input_fw, 'fasta'):
        fwd_records.append(
            SeqRecord(record.seq, id=record.id, description='')
        )
    for record in SeqIO.parse(input_rev, 'fasta'):
        rev_records.append(
            SeqRecord(record.seq.reverse_complement(), id=record.id, description='')
        )
    SeqIO.write(fwd_records, fwd, 'fasta')
    SeqIO.write(rev_records, rev, 'fasta')


def total_reads(fastq_dest):
    """Counts the total reads in a fastq file

    Arguments:
    fastq_dest-- a fasta file containing reads
    """
    fastq = fastq_dest
    with open(fastq, 'r') as f:
        reads = []
        for line in f:
            reads.append(line)
        return (len(reads)/4)

def trim_primers(input, trimmed, primerless, tooshort, prim_primers,  sec_primer, output_namer, cutoff, prim_three_or_five_primer):
    """Trim primer sequences from the sequencing reads.
    The reads without a primer are excluded and output to a separate file. The reads that are too short (less than 50
    bases) after trimming are also excluded and output to a separate file.

    Arguments:
    input -- the fastq file containing the sequencing reads
    trimmed -- the file path for the trimmed reads
    primerless -- the file path for the reads that do not have a primer
    tooshort -- the file path for the reads that are too short after trimming
    prim_primers -- 
    sec_primer -- 
    output_namer -- 
    cutoff -- 
    prim_three_or_five_primer -- 
    """
    prim_primer_call = str(prim_primers)
    sec_primer_call = str(sec_primer)


    if prim_three_or_five_primer == 'three':
       prim_orientation = '-a'
       sec_orientation = '-g'
    elif prim_three_or_five_primer == 'five':
       prim_orientation = '-g'
       sec_orientation = '-a'

    else:
       prim_orientation = '-b'
       sec_orientation = '-b'
  

    subprocess.check_call([
        'cutadapt',
        prim_orientation, prim_primer_call,
        '-e', '0.25', 
        '--untrimmed-output', primerless,
        '-o', output_namer,
        input
    ])

    subprocess.check_call([
        'cutadapt',
        sec_orientation, sec_primer_call,
        '-e', '0.25',
        '-m', str(cutoff),
        '--too-short-output', tooshort,
        '-o', trimmed,
        output_namer
    ])


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True,
              help='fastq file path for sequence reads')
@click.option('-o', '--output', type=click.Path(exists=False), required=True,
              help='output CSV file path')
@click.option('--temp_folder', type=click.Path(exists=False), required=True,
              help='temp folder for intermediate files')
@click.option('--forward_primers', type=click.Path(exists=True), required=True,
              help='forward primers fasta file path')
@click.option('--reverse_primers', type=click.Path(exists=True), required=True,
              help='reverse primers fasta file path')
@click.option('--taxmethod',
              type=click.Choice(['BLAST', 'DADA2'], case_sensitive=False),
              required=True, help='the method for taxonomy assignment')
@click.option('--taxreference', show_default=True, default='UNITE_eukaryote',
              type=click.Choice(dada2_tax_dbs.keys(), case_sensitive=False),
              help='the taxonomy reference database to be used if assigning '
              'taxonomy using DADA2')
@click.option('--blastdatabase', type=click.Path(), show_default=True,
              default='/gpfs_partners/databases/ncbi/blast/nt/nt',
              help='the path for the BLAST database to be used if assigning '
              'taxonomy using BLAST')
@click.option('--threads', type=int, default=4, show_default=True,
              help='the threads/cpus to be used for running multithreaded tasks')
@click.option('--cutoff', type=int, default=50, show_default=True,
              help='cutoff length for sequences after primer trimming')
@click.option('--unique_run_id', type=str, default='', show_default=True,
              help='id to differentiate files in different runs')
@click.option('--trim_both', type=str, default='both_rev', show_default=True,
              help='id to differentiate files in different runs')
@click.option('--prim_three_or_five_primer', type=str, default='three', show_default=True,
              help='id to differentiate files in different runs')
@click.option('--reverse_comp', type=str, default='True', show_default='True',
              help='id to differentiate files in different runs')



def main(input, output, temp_folder, forward_primers, reverse_primers, taxmethod, taxreference, blastdatabase, threads, cutoff, unique_run_id, trim_both, prim_three_or_five_primer, reverse_comp):
    """Identify ASVs and assign taxonomy.
    This is the main function for the app and it's arguments are specified
    from the command line. The output is a CSV file containing the ASVs,
    their abundance, and their taxonomy assignment.\f

    Arguments:
    input -- fastq file path for sequence reads
    output -- output CSV file path
    primers -- primers fasta file path
    taxmethod -- the method for taxonomy assignment
    taxreference -- the taxonomy reference database to be used if assigning
      taxonomy using DADA2
    blastdatabase -- the path for the BLAST database to be used if assigning
      taxonomy using BLAST
    threads -- the threads/cpus to be used for running multithreaded tasks
    cutoff -- cutoff length for sequences after primer trimming
    """
    TMP = temp_folder
    os.makedirs(TMP, exist_ok=True)


    # set input file base name and get output directory
    base = os.path.basename(input).replace('.gz', '').replace('.fastq', '')
    trunc_base_arr = base.split('_')
    trun_base = trunc_base_arr[0]
    outdir = os.path.dirname(output)

    # set file names for primers and primer trimming
    fwdprimers = f'{TMP}/{base}_forwardprimerstemp.fasta'
    revprimers = f'{TMP}/{base}_reverseprimerstemp.fasta'
    trimmed = f'{TMP}/{base}_trimmed.fastq'
    primerless = f'{outdir}/{base}_untrimmed.fastq'
    tooshort = f'{outdir}/{base}_tooshort.fastq'
    
    parse_primers(forward_primers, reverse_primers, fwdprimers, revprimers)

    temp_input=input
    temp_output = output
    primer_inputer = open(f'{revprimers}')
    primer_fwd_aller = open(f'{fwdprimers}')
    fw_primer_list = list(SeqIO.parse(primer_fwd_aller,'fasta'))
    rev_primer_list = list(SeqIO.parse(primer_inputer,'fasta'))
    trim_land=[]
    too_short_land = []
    primerless_land = []

    primary_primer_arr=[]
    secondary_primer_arr =[]

    for x in range(0,len(rev_primer_list)):
        #temp variables
        sing_primer = rev_primer_list[x]
        sec_primer = rev_primer_list[x]

        rev_sing_primer = rev_primer_list[x]
        fwd_sing_primer = fw_primer_list[x]
       
        if trim_both == "both_rev":
            sing_primer = rev_sing_primer
            sec_primer = fwd_sing_primer
        elif trim_both == "both_fwd":
            sing_primer = fwd_sing_primer
            sec_primer = rev_sing_primer
        elif trim_both == "fwd":
            sing_primer = fwd_sing_primer
            sec_primer = fwd_sing_primer
        else:
            sing_primer = rev_sing_primer
            sec_primer = rev_sing_primer

        primer_indic_arr = sing_primer.id.split("_")
        primer_indic = primer_indic_arr[0]
        if primer_indic in base:
                
            trimmed_temp = f'{TMP}/{base}_trimmed_'   + sing_primer.id +'_' + unique_run_id + '.fastq'
            trimmed = f'{trimmed_temp}'
            trim_land.append(trimmed)

            primerless_temp = f'{outdir}/{base}_untrimmed_after_' + sing_primer.id + '_removal' + '_' + unique_run_id+ '.fastq'
            primerless = f'{primerless_temp}'
            primerless_land.append(primerless)

            tooshort_temp = f'{outdir}/{base}_tooshort_' + sing_primer.id +'_' + unique_run_id + '.fastq'
            tooshort = f'{tooshort_temp}'
            too_short_land.append(tooshort)

            output_namer = f'{TMP}/{base}_temp_primary_primer_trimmed_' +  sing_primer.id +'_' + unique_run_id + '.fastq'

            curr_prim_primer = sing_primer.seq
            curr_sec_primer = sec_primer.seq

            if reverse_comp == "True":
                
                if trim_both == "both_rev":
                    curr_prim_primer = sing_primer.seq.reverse_complement()
                    curr_sec_primer = sec_primer.seq.reverse_complement()
                elif trim_both == "both_fwd":
                    curr_prim_primer = sing_primer.seq.reverse_complement()
                    curr_sec_primer = sec_primer.seq.reverse_complement()
                elif trim_both == "fwd":
                    curr_prim_primer = sing_primer.seq.reverse_complement()
                    curr_sec_primer = "GGGGGGGGGGGGGGGGGGG"
                else:
                    curr_prim_primer = sing_primer.seq.reverse_complement()
                    curr_sec_primer = "GGGGGGGGGGGGGGGGGGG"
            else:
                if trim_both == "both_rev":
                    curr_prim_primer = sing_primer.seq
                    curr_sec_primer = sec_primer.seq
                elif trim_both == "both_fwd":
                    curr_prim_primer = sing_primer.seq
                    curr_sec_primer = sec_primer.seq
                elif trim_both == "fwd":
                    curr_prim_primer = sing_primer.seq
                    curr_sec_primer = "GGGGGGGGGGGGGGGGGGG"
                else:
                    curr_prim_primer = sing_primer.seq
                    curr_sec_primer = "GGGGGGGGGGGGGGGGGGG"

            primary_primer_arr.append(sing_primer)

            secondary_primer_arr.append(sec_primer)

            trim_primers(temp_input, trimmed, primerless, tooshort, curr_prim_primer, curr_sec_primer, output_namer, cutoff, prim_three_or_five_primer)
            temp_input = input

    for y in range(0,len(trim_land)):
        sing_primer = primary_primer_arr[y]
        trimmed = trim_land[y]
        primerless = primerless_land[y]

        if((os.path.isfile(trimmed) == True) and (int(total_reads(trimmed))>100)):

            # dada2 asv identification
            print('Identifying ASVs with DADA2 ...')
            asv = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_asv.csv'
            feat_table = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_seqtab-nochim.txt'
            rep_seqs = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_rep-seqs.fna'
            asv_data = run_dada2(trimmed, asv, feat_table, rep_seqs, TMP)

            if((os.path.isfile(asv) == True) and (len(asv_data) > 0)):

                #file wraggling for qiime input
                biom_table = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_biom-table.txt'
                table_biom_format = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_table.biom'
                cmd_biom = 'echo -n "#OTU Table" | cat - ' + str(feat_table) + ' > ' + str(biom_table) 
                subprocess.check_call(cmd_biom, shell=True)

                cmd_biom_convert = 'biom convert -i ' + str(biom_table) + ' -o ' + str(table_biom_format) + ' --table-type="OTU table" --to-hdf5 '
                #biom convert -i dada2-analysis/biom-table.txt -o dada2-analysis/table.biom --table-type="OTU table" --to-hdf5
                subprocess.check_call(cmd_biom_convert, shell=True)

                table_qza_address = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_table.qza'
                cmd_table_qiime_import = 'qiime tools import --input-path ' + str(table_biom_format) + " --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path " + str(table_qza_address)
                subprocess.check_call(cmd_table_qiime_import, shell=True)
                
                # taxonomic classification
                # via blast
                if taxmethod == 'BLAST':

                    # prep input fasta file for blast
                    asv_fasta = f'{TMP}/{base}_'+  sing_primer.id +'_' + unique_run_id + '_asv.fasta'
                    write_blast_fasta(asv_data, asv_fasta)
         
                    # run blast
                    print('Running BLAST ...')
                    blast_results = f'{TMP}/{base}_'+  sing_primer.id + '_' + unique_run_id + '_blast.tsv'
                    blast_data = run_blast(asv_fasta, blast_results, blastdatabase, threads)

                    if(os.path.isfile(blast_results) == True and (len(blast_data) > 0)):
                        # look up taxa
                        print('Running Taxize ...')
                        taxid_list = f'{TMP}/{base}_' +  sing_primer.id + '_' + unique_run_id + '_taxids.csv'
                        blast_data[['staxid']].drop_duplicates().to_csv(taxid_list, header=False, index=False)
                        taxize_results = f'{TMP}/{base}_'+  sing_primer.id + '_' + unique_run_id + '_taxa.csv'

                        if(os.path.isfile(taxid_list) == True):
                            taxa = run_taxize(taxid_list, taxize_results)
                            if((os.path.isfile(taxize_results) == True) and (len(taxa) > 0)):
                                # join blast and taxonomy information
                                final_taxa = combine_blast_taxize(blast_data, taxa)

                                # join taxize taxonomy results with ASV data
                                temp_output = pd.merge(
                                    asv_data,
                                    blast_data[['qacc','stitle']].rename(columns={'qacc': 'index'}),
                                    on='index',
                                    how='left'
                                 )

                                output_data = pd.merge(
                                    temp_output,
                                    final_taxa.rename(columns={'qacc': 'index'}),
                                    on='index',
                                    how='left'
                                 ).drop('index', axis=1)
                                temp_output = f'{outdir}/{trun_base}' + '.csv'

                                temp_seq = ""
                                row_ind_to_drop = []
                                for index, row in output_data.iterrows():
                                    if temp_seq != row['sequence']:
                                        temp_seq = row['sequence']
                                    else:
                                        row_ind_to_drop.append(index)
                                output_data = output_data.drop(row_ind_to_drop)
         
                                # output table with ASVs, abundance, and taxonomy
                                output_data.to_csv(temp_output, index=False)
                                print(f'Analysis complete. Find results at {temp_output}')
                            else:
                                err_output = 'ERROR: No taxize output, analysis failure at taxize step for primer set ' +  sing_primer.id
                                print(err_output)
                        else:
                            err_output = 'ERROR: No taxize list output, analysis failure at taxize list generation step for primer set ' +  sing_primer.id
                            print(err_output)
                    else:
                        err_output = 'ERROR: No blast output, analysis failure at ASV blast step for primer set ' +  sing_primer.id
                        print(err_output)
                else:
                    print("ERROR: That is not a valid Taxmethod option")
            else:
                err_output = 'ERROR: No dada2 output, analysis failure at ASV generation step for primer set ' +  sing_primer.id
                print(err_output)
        else:
            print(('ERROR: No reads passed the filter for primer ' + sing_primer.id))


        if((os.path.isfile(primerless) == True) and (int(total_reads(primerless))>100)):

            # dada2 asv identification
            print('Identifying ASVs with DADA2 ...')
            asv = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_asv.csv'
            feat_table = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_seqtab-nochim.txt'
            rep_seqs = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_rep-seqs.fna'
            asv_data = run_dada2(primerless, asv, feat_table, rep_seqs, TMP)


            if((os.path.isfile(asv) == True) and (len(asv_data) > 0)):

                #file wraggling for qiime input
                biom_table = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_biom-table.txt'
                table_biom_format = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_table.biom'
                cmd_biom = 'echo -n "#OTU Table" | cat - ' + str(feat_table) + ' > ' + str(biom_table) 
                subprocess.check_call(cmd_biom, shell=True)

                cmd_biom_convert = 'biom convert -i ' + str(biom_table) + ' -o ' + str(table_biom_format) + ' --table-type="OTU table" --to-hdf5 '
                #biom convert -i dada2-analysis/biom-table.txt -o dada2-analysis/table.biom --table-type="OTU table" --to-hdf5
                subprocess.check_call(cmd_biom_convert, shell=True)

                table_qza_address = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_table.qza'
                cmd_table_qiime_import = 'qiime tools import --input-path ' + str(table_biom_format) + " --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path " + str(table_qza_address)
                subprocess.check_call(cmd_table_qiime_import, shell=True)
                
               # taxonomic classification
                # via blast
                if taxmethod == 'BLAST':

                    # prep input fasta file for blast
                    asv_fasta = f'{TMP}/{base}_untrimmed_'+  sing_primer.id +'_' + unique_run_id + '_asv.fasta'
                    write_blast_fasta(asv_data, asv_fasta)
         
                    # run blast
                    print('Running BLAST ...')
                    blast_results = f'{TMP}/{base}_untrimmed_'+  sing_primer.id + '_' + unique_run_id + '_blast.tsv'
                    blast_data = run_blast(asv_fasta, blast_results, blastdatabase, threads)

                    if((os.path.isfile(blast_results) == True) and (len(blast_data) > 0)):
                        # look up taxa
                        print('Running Taxize ...')
                        taxid_list = f'{TMP}/{base}_untrimmed_' +  sing_primer.id + '_' + unique_run_id + '_taxids.csv'                        
                        blast_data[['staxid']].drop_duplicates().to_csv(taxid_list, header=False, index=False)
                        taxize_results = f'{TMP}/{base}_untrimmed_'+  sing_primer.id + '_' + unique_run_id + '_taxa.csv'
                        if(os.path.isfile(taxid_list) == True):
                            taxa = run_taxize(taxid_list, taxize_results)
                            if((os.path.isfile(taxize_results) == True) and (len(taxa) > 0)):
                                # join blast and taxonomy information
                                final_taxa = combine_blast_taxize(blast_data, taxa)

                                # join blast taxonomy results with ASV data
                                output_data = pd.merge(
                                    asv_data,
                                    final_taxa.rename(columns={'qacc': 'index'}),
                                    on='index',
                                    how='left'
                                 ).drop('index', axis=1)
                           
                                temp_output = f'{outdir}/{trun_base}_untrimmed' + '.csv'
                                # output table with ASVs, abundance, and taxonomy

                                temp_seq = ""
                                row_ind_to_drop = []
                                for index, row in output_data.iterrows():
                                    if temp_seq != row['sequence']:
                                        temp_seq = row['sequence']
                                    else:
                                        row_ind_to_drop.append(index)
                                output_data = output_data.drop(row_ind_to_drop)


                                output_data.to_csv(temp_output, index=False)
                                print(f'Analysis complete. Find results at {temp_output}')
                            else:
                                err_output = 'ERROR: No taxize output, analysis failure at taxize step for primer set ' +  sing_primer.id
                                print(err_output)
                        else:
                            err_output = 'ERROR: No taxize list output, analysis failure at taxize list generation step for primer set ' +  sing_primer.id
                            print(err_output)
                    else:
                        err_output = 'ERROR: No blast output, analysis failure at ASV blast step for primer set ' +  sing_primer.id
                        print(err_output)
                else:
                    print("ERROR: That is not a valid Taxmethod option")
            else:
                err_output = 'ERROR: No dada2 output, analysis failure at ASV generation step for primer set ' +  sing_primer.id
                print(err_output)
        else:
            print(('ERROR: No reads passed the filter for primer ' + sing_primer.id))



if __name__ == '__main__':
    main()
