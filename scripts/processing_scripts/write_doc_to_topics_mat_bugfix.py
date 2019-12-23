#!/usr/bin/env python
'''
DESCRIPTION

    

FOR HELP

    python write_doc_to_topics_mat_bugfix.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2019-07-20
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import gensim
import scipy.io
import numpy as np
import csv

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infilemm', metavar='INFILE',
                        help='mm file')
    parser.add_argument('infilelda', metavar='INFILE',
                        help='LDA file')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='doc to topics output')
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

    scipy_sparse_matrix=scipy.io.mmread(args.infilemm)
    corpus = gensim.matutils.Sparse2Corpus(scipy_sparse_matrix)
    lda = gensim.models.LdaModel.load(args.infilelda)

    all_topics = lda.get_document_topics(corpus, per_word_topics=False, minimum_probability=0, minimum_phi_value=0)  # zeros still omitted 

    with open(args.outfile, mode="w") as wf:
        writer = csv.writer(wf, delimiter = "\t")
        for d in all_topics:
            # document list length is not necessarily num_topics
            topic_i = 0  # track topics
            outrow = []
            for tup in d:
                topic, weight = tup[0], tup[1]
                assert topic_i <= topic  # otherwise something is wrong!
                while topic > topic_i:  # means we skipped some zero-weight topics
                    outrow.append(0)
                    topic_i += 1
                if topic == topic_i:
                    outrow.append(weight)
                    topic_i += 1
            while topic_i < lda.num_topics:
                outrow.append(0)
                topic_i += 1
            # topic_i should always finish at lda.num_topics - 1
            assert topic_i == lda.num_topics
            assert len(outrow) == lda.num_topics
            assert(abs(sum(outrow) - 1) < 0.00001)  # check things sum up to ~1
            writer.writerow(outrow)


if __name__ == '__main__':
    main()
