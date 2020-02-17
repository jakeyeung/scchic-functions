#!/usr/bin/env python
'''
DESCRIPTION

    Merge huge motevo outputs (across chunks and motifs)
    For Mara_DHS we dont do across chromosomes because they are stil too large

FOR HELP

    python merge_motevo_output.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2015-06-03
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse
import os
import csv

def get_list_from_file(fname, col_i = 0):
    '''
    Read first column in fname and return everything as a list
    '''
    lst = []
    with open(fname, 'rb') as readfile:
        reader = csv.reader(readfile, delimiter = '\t')
        for row in reader:
            lst.append(row[col_i])
    return(lst)

def get_site_files(indir, wmnames, dirs):
    site_files = {}
    for wm in wmnames:
        # open output file in outdir
        outname = '.'.join([wm, "combined.sites"])
        if wm not in site_files:
            site_files[wm] = []
        else:
            print 'Warning: %s already has %s' % (site_files, wm)
        for jdir in dirs:
            fname = '.'.join([wm, "sites"])
            fpath = os.path.join(indir, jdir, fname)
            if os.path.isfile(fpath):
                site_files[wm].append(fpath)
            else:
                print 'Warning: could not find %s' % fpath
    return site_files

class MotRow:
    '''
    Class for handling motevo output row
    '''
    def __init__(self, row):
        self.region = row[0]
        self.strand = row[1]
        self.prob = float(row[2])
        self.motif = row[3]
        self.motevo_id = row[4]

    def merge_row(self, Row):
        '''
        Update parameters by merging with another Row object
        '''
        self.region = ','.join([self.region, Row.region])
        self.strand = ','.join([self.strand, Row.strand])
        self.prob += Row.prob

    def write_to_file(self, writer):
        writer.writerow([self.region, self.strand, self.prob, self.motif, self.motevo_id])

def merge_path(reader, writer):
    '''
    If row id matches next row id, then merge it together by adding up the rows then
    write it to outfile
    '''
    for rowcount, row in enumerate(reader):
        if rowcount == 0:
            PrevRow = MotRow(row)
            try:
                row = reader.next()
            except StopIteration:
                PrevRow.write_to_file(writer)

        Row = MotRow(row)
        if PrevRow.motevo_id == Row.motevo_id:
            # merge rows because they are same id
            PrevRow.merge_row(Row)
        else:
            # write row because they are different
            PrevRow.write_to_file(writer)
            # reinitiate PrevRow with current row
            PrevRow = Row

def write_rows(reader, writer):
    '''
    If row id matches next row id, then merge it together by adding up the rows then
    write it to outfile
    '''
    for rowcount, row in enumerate(reader):
        Row = MotRow(row)
        Row.write_to_file(writer)

def main():
    parser = argparse.ArgumentParser(description='Merge huge motevo outputs (across chromosomes and across motifs)')
    parser.add_argument('indir', metavar='INDIR',
                        help='Input directory containing chunk dirs and motevo output')
    parser.add_argument('outdir', metavar='OUTDIR',
                        help='Output directory containing merged motevo output genome wide')
    parser.add_argument('--wmnamesfile', metavar='FILE', required = True,
                        help='File containing names of weight matrices.'\
                            'Sanity check that all the motifs are present')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--merge', '-m', action='store_true',
                       help='Merge rows with common IDs and add up their probabilities')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    dirs = os.listdir(args.indir)
    wmnames = get_list_from_file(args.wmnamesfile)

    if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

    site_files = get_site_files(args.indir, wmnames, dirs)

    for wm, paths in site_files.iteritems():
        outname = '.'.join([wm, 'combined.sites'])
        outpath = os.path.join(args.outdir, outname)
        with open(outpath, 'wb') as outfile:
            writer = csv.writer(outfile, delimiter = ' ')
            for p in paths:
                with open(p, 'rb') as readfile:
                   reader = csv.reader(readfile, delimiter = ' ')
                   if args.merge:
                       merge_path(reader, writer)
                   else:
                       write_rows(reader, writer)


if __name__ == '__main__':
    main()
