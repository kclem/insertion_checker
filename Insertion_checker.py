import argparse
import os
import pandas as pd

def main():
    parser = argparse.ArgumentParser(prog='Insertion Checker'))
    parser.add_argument('-s', '--spacer_seq', type=str, help='Spacer sequence', required=True)
    parser.add_argument('-c', '--cut_offset', type=int, help='Cut offset in spacer where insertions are to be counted', default=-3)
    parser.add_argument('-b','--barcode_file', type=str, help='Barcode sequence file', required=True)
    parser.add_argument('-o', '--output', type=int, help='Output file', default=None)
    parser.add_argument('--num_bp_pre', type=int, help='Number of bases from the spacer pre-cut to match', default=4)
    parser.add_argument('--num_bp_post', type=int, help='Number of bases from the spacer post-cut to match', default=4)
    args = parser.parse_args()

    if not os.path.exists(args.barcode_file):
        print('Barcode file does not exist')
        return

    #read barcodes from barcode file which is an excel file
    barcode_df = pd.read_excel(args.barcode_file, header=None)
    barcode_df.columns = ['barcode']
    


if __name__ == '__main__':
    main()
