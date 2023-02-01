# Insertion Checker
This tool identifies and counts barcodes that have been inserted into genomic positions.


## Usage:
```
usage: Insertion Checker [-h] -s SPACER_SEQ [-c CUT_OFFSET] -b BARCODE_FILE [-o OUTPUT] -r1 FASTQ_R1 [-r2 FASTQ_R2]
                         [--num_bp_pre NUM_BP_PRE] [--num_bp_post NUM_BP_POST]
                         [--min_paired_end_reads_overlap MIN_PAIRED_END_READS_OVERLAP] [--debug]

optional arguments:
  -h, --help            show this help message and exit
  -s SPACER_SEQ, --spacer_seq SPACER_SEQ
                        Spacer sequence
  -c CUT_OFFSET, --cut_offset CUT_OFFSET
                        Cut offset in spacer where insertions are to be counted
  -b BARCODE_FILE, --barcode_file BARCODE_FILE
                        Barcode sequence file with the column 'Barcode' containing the barcodes. This can be a .txt tab-separated, .csv comma-separated, or xlsx file (the data must be in the first tab)
  -o OUTPUT, --output OUTPUT
                        Output file
  -r1 FASTQ_R1, --fastq_r1 FASTQ_R1
                        Input fastq r1 file
  -r2 FASTQ_R2, --fastq_r2 FASTQ_R2
                        Input fastq r2 file
  --num_bp_pre NUM_BP_PRE
                        Number of bases from the spacer pre-cut to match
  --num_bp_post NUM_BP_POST
                        Number of bases from the spacer post-cut to match
  --min_paired_end_reads_overlap MIN_PAIRED_END_READS_OVERLAP
                        Minimum bp overlap between paired end reads for read merging
  --debug               Whether to output debug information
  ```
