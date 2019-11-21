#!/usr/bin/env python
'''
DESCRIPTION

    Bam input, add "chr" to chromosome names

FOR HELP

    python reformat_bam_header.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2019-09-08
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam

def reformat_header(header, add_prefix="chr"):
    new_header = header.to_dict()
    contigs = new_header['SQ']
    new_contigs = []
    for contig in contigs:
        contig['SN'] = ''.join([add_prefix, contig['SN']])
        new_contigs.append(contig)
    new_header['SQ'] = new_contigs
    return pysam.AlignmentHeader.from_dict(new_header)


def main():
    parser = argparse.ArgumentParser(description='Bam input, add "chr" to chromosome names')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam file')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output bam file')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    with pysam.AlignmentFile(args.infile, "rb") as inbam:
        new_header = reformat_header(inbam.header, add_prefix = "chr")
        with pysam.AlignmentFile(args.outfile, "wb", header = new_header) as outbam:
            for read in inbam:
                outbam.write(read)
        pysam.index(args.outfile)


if __name__ == '__main__':
    main()
