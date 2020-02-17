#!/usr/bin/env python
'''
DESCRIPTION

    Assign closest bed rather than overwriting. Appends to end of columns.

FOR HELP

    python overwrite_nearest_gene_bed.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2016-04-02
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, csv
import cPickle as pickle


class BedRow():
    def __init__(self, row, indx=False, motevo_coord=False):
        '''
        is_closest contains additional columns.
        Columns 6 to 11 contain closestbed info
        closestchri: column index containing chromo of closest bed
        '''
        self.chromo = row[0]
        self.start = row[1]
        self.end = row[2]
        self.startend = '-'.join([self.start, self.end])
        if not motevo_coord:
            self.coord = ':'.join([self.chromo, self.startend])
        else:
            # coord is found in 4th column must parse
            # REST.p3;mm10_chr1:68258053-68258553 -> chr1:68258053-68258553
            self.motevo_id = row[3]
            self.motevo_id = self.motevo_id.split(";")[1]
            self.coord = self.motevo_id.split("_")[1]
            # check that coord is: M098_0.6;mm10_chr1:36320203-36324029(+)
            # or M098_0.6;mm10_chr1:36320203-36324029(-)
            # if so, remove the (+) or (-)
            if self.coord.endswith("(+)") or self.coord.endswith("(-)"):
                # remove last 3 characters, which corresponds to (+) or (-)
                self.coord = self.coord[:-3]

        if indx:
            # if bed is an indx bed
            self.name = row[3]
            self.value = row[4]

        if not self.coord.startswith("chr"):
            print("Exiting, coords not correct!")
            print(coord)
            sys.exit()


def index_bed(inbed):
    '''
    Index bed {chr1:start:end: genes, dists}
    '''
    annots = {}
    with open(inbed, 'rb') as readfile:
        bedreader = csv.reader(readfile, delimiter = '\t')
        for row in bedreader:
            Row = BedRow(row, indx=True)
            annots[Row.coord] = [Row.name, Row.value]
    return(annots)


def assign_to_bed(bedin, new_genes, bedout):
    '''
    Overwrite annotations of closestBed to new_genes
    '''
    with open(bedin, 'rb') as readfile, open(bedout, 'wb') as writefile:
        bedreader = csv.reader(readfile, delimiter = '\t')
        bedwriter = csv.writer(writefile, delimiter = '\t')
        for row in bedreader:
            Row = BedRow(row, indx=False, motevo_coord=True)
            try:
                new_annot = new_genes[Row.coord]
                # print(Row.coord)
                # raw_input()
            except KeyError:
                continue
            # overwrite
            for i in [0, 1]:  # gene, dist
                row.append(new_annot[i])
            bedwriter.writerow(row)


def main():
    parser = argparse.ArgumentParser(description='Given closestBed output for closest gene, modify to closest genes within interval')
    parser.add_argument('bedin', metavar='BEDFILE IN',
                        help='Bed input. Default first 3 columns are chromo, start, end')
    parser.add_argument('new_annots_bed', metavar='BEDFILE IN',
                        help='Bed file containing new annotations')
    parser.add_argument('bedout', metavar='BEDFILE OUT',
                        help='Bed output, gene and distance assignments append at end')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--has_motevo_id', '-m', action='store_true',
                        help='Flag if coord to match is in 4th column as STRING;mm10_coord')
    parser.add_argument('--save_pickle', '-s', action='store_true',
                        help='Save pickle to nohup tmp')
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

    print("Indexing bed")
    new_genes = index_bed(args.new_annots_bed)
    if args.save_pickle:
        pklfile = ''.join([args.bedout, ".pkl"])
        pickle.dump(new_genes, open(pklfile, "wb"))
    print("Assigning index to new bed")
    assign_to_bed(args.bedin, new_genes, args.bedout)

if __name__ == '__main__':
    main()
