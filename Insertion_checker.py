import argparse
from collections import defaultdict
import gzip
import logging
import os
import pandas as pd
import re
import subprocess as sb

def main():
    parser = argparse.ArgumentParser(prog='Insertion Checker')
    parser.add_argument('-s', '--spacer_seq', type=str, help='Spacer sequence', required=True)
    parser.add_argument('-c', '--cut_offset', type=int, help='Cut offset in spacer where insertions are to be counted', default=-3)
    parser.add_argument('-b','--barcode_file', type=str, help='Barcode sequence file', required=True)
    parser.add_argument('-o', '--output', type=int, help='Output file', default=None)
    parser.add_argument('-r1', '--fastq_r1', type=str, help='Input fastq r1 file', required=True)
    parser.add_argument('-r2', '--fastq_r2', type=str, help='Input fastq r2 file', default=None)
    parser.add_argument('--num_bp_pre', type=int, help='Number of bases from the spacer pre-cut to match', default=6)
    parser.add_argument('--num_bp_post', type=int, help='Number of bases from the spacer post-cut to match', default=6)
    parser.add_argument('--min_paired_end_reads_overlap', type=int, help='Minimum bp overlap between paired end reads for read merging', default=20)
    parser.add_argument('--debug', help='Whether to output debug information', action='store_true')
    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


    # create console handler and set level to debug
    ch = logging.StreamHandler()

    # create formatter
    formatter = logging.Formatter('%(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    # check if input files exist
    if not os.path.exists(args.fastq_r1):
        logger.critical('Input fastq r1 file does not exist')
        return

    if not os.path.exists(args.barcode_file):
        logger.critical('Barcode file does not exist')
        return

    # set output file
    output_file = args.output
    if args.output is None:
        output_file = args.fastq_r1 + '.barcode_counts.txt'

    # merge paired end reads
    input_fastq = args.fastq_r1
    if args.fastq_r2 is not None:
        logger.debug('Merging paired end reads...')
        input_fastq = output_file+".merged.fq"
        unmerged_r1 = output_file+".merged_r1.fq"
        input_fastq = output_file+".merged_r1.fq"
        merge_cmd = '{command} -i {r1} -I {r2} --merge --merged_out {out_merged} --unpaired1 {unpaired1} --unpaired2 {unpaired2} --overlap_len_require {min_overlap} --html {html_report}'.format(
                command='fastp',
                r1=args.fastq_r1,
                r2=args.fastq_r2,
                out_merged=input_fastq,
                unpaired1=output_file+".unmerged_r1.fq",
                unpaired2=output_file+".unmerged_r2.fq",
                min_overlap=args.min_paired_end_reads_overlap,
                html_report=output_file + '.fastp_report.html',
            )

        fastp_status = sb.call(merge_cmd, shell=True)
        if fastp_status:
            logger.critical('Fastp failed to run, please check the log file.')
            return

        if not os.path.isfile(input_fastq):
            logger.critical('Fastp merging failed, please check the log file.')

    #read barcodes from barcode file
    if args.barcode_file.endswith('.xlsx') or args.barcode_file.endswith('.xls'):
        barcode_df = pd.read_excel(args.barcode_file)
    elif args.barcode_file.endswith('.csv'):
        barcode_df = pd.read_csv(args.barcode_file)
    elif args.barcode_file.endswith('.txt'):
        barcode_df = pd.read_csv(args.barcode_file, sep='\t')

    if 'Barcode' not in barcode_df.columns:
        logger.critical('Barcode column not found in barcode file. Make sure your barcode file has a column named "Barcode".')
        return

    barcode_lens = barcode_df['Barcode'].str.len()
    if min(barcode_lens) != max(barcode_lens):
        logger.critical(f'Barcodes are not of the same length ({min(barcode_lens)} vs {max(barcode_lens)})')
        return

    valid_barcodes = dict.fromkeys(barcode_df['Barcode'].tolist(),1)

    left_bases_to_check = args.spacer_seq[args.cut_offset-args.num_bp_pre:args.cut_offset]
    right_bases_to_check = args.spacer_seq[args.cut_offset:]
    right_bases_to_check = right_bases_to_check[:args.num_bp_post]

    logger.debug('spacer seq: ' + args.spacer_seq)
    logger.debug('left bases: ' + str(left_bases_to_check) + ' right bases: ' + str(right_bases_to_check))
    logger.debug('barcode regex: ' + f'{left_bases_to_check}([ATGC]{{{max(barcode_lens)}}}){right_bases_to_check}')

    barcode_regex = re.compile(f'{left_bases_to_check}([ATGC]{{{max(barcode_lens)}}}){right_bases_to_check}')

    #read fastq file
    no_barcode_count = 0
    invalid_barcode_count = 0 #barcode is correct length but not in barcode list
    valid_barcode_count = 0

    read_seq_count = 0
    barcode_counts = defaultdict(int)
    invalid_barcode_counts = defaultdict(int)
    if input_fastq.endswith('.gz'):
        fin = gzip.open(input_fastq, 'rt')
    else:
        fin = open(input_fastq, 'r')

    while(True):
        info_line = fin.readline()
        seq_line = fin.readline()
        plus_line = fin.readline()
        qual_line = fin.readline()
        if not info_line:
            break
        read_seq_count += 1

        #check to see if sequence contains barcode_regex
        match = barcode_regex.search(seq_line)
        if match:
            barcode = match.group(1)
            barcode_counts[barcode] += 1
            if barcode in valid_barcodes:
                valid_barcode_count += 1
            else:
                invalid_barcode_count += 1
                invalid_barcode_counts[barcode] += 1

        else:
            no_barcode_count += 1
    fin.close()

    barcode_df['Count'] = [barcode_counts[x] for x in barcode_df['Barcode']]

    barcode_df.to_csv(output_file, sep='\t', index=False)
    invalid_barcode_file = output_file + '.invalid_barcodes.txt'
    with open(invalid_barcode_file, 'w') as fout:
        for barcode in invalid_barcode_counts:
            fout.write(f'{barcode}\t{invalid_barcode_counts[barcode]}\n')


    invalid_counts = len(invalid_barcode_counts.keys())
    valid_counts = len(barcode_counts.keys())
    logger.info(f'Finished processing {read_seq_count} reads (valid barcodes: {valid_barcode_count}, invalid barcodes: {invalid_barcode_count}, unidentified reads: {no_barcode_count})')
    logger.info(f'Printed {valid_counts} valid barcode counts to {output_file}')
    logger.info(f'Printed {invalid_counts} invalid barcodes to {invalid_barcode_file}')


if __name__ == '__main__':
    main()
