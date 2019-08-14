#!/usr/bin/env python
'''
DESCRIPTION

    Quality control plots after running LDA on gensim

FOR HELP

    python plot_likelihood_across_iterations_gensim.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2019-08-06
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

# https://askubuntu.com/questions/1045720/what-is-a-good-default-backend-for-matplotlib
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
plt.ioff() #http://matplotlib.org/faq/usage_faq.html (interactive mode)
plt.show._needmain=False
import sys, argparse, datetime
import gensim
import re
import scipy.io
from gensim import corpora
from gensim.test.utils import datapath



def main():
    parser = argparse.ArgumentParser(description='Quality control plots after running LDA using gensim VI')
    parser.add_argument('infile', metavar='INFILE',
                        help='Path to error file containing likelihood and perplexity statements every couple of lines')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output main. One is .png output and one .txt output.')
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

    plotfname = ''.join([args.outfile, ".png"])
    txtfname = ''.join([args.outfile, ".txt"])

    p = re.compile("(-*\d+\.\d+) per-word .* (\d+\.\d+) perplexity")
    matches = [p.findall(l) for l in open(args.infile)]
    matches = [m for m in matches if len(m) > 0]
    tuples = [t[0] for t in matches]
    perplexity = [float(t[1]) for t in tuples]
    liklihood = [float(t[0]) for t in tuples]
    iter = list(range(0,len(tuples)*10,10))
    if len(iter) == 0:
        print("No text matches, exiting")
        sys.exit()

    plt.subplot(1, 2, 1)

    plt.plot(iter,liklihood,c="black")
    plt.ylabel("log liklihood")
    plt.xlabel("iteration")
    plt.title("Topic Model Convergence")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(iter,perplexity,c="black")
    plt.ylabel("perplexity")
    plt.xlabel("iteration")
    plt.title("Topic Model Convergence")
    plt.grid()
    plt.savefig(plotfname)

    # write to output
    with open(txtfname, "w") as outf:
        for i in range(len(liklihood)):
            rowint = [iter[i], liklihood[i], perplexity[i]]
            rowstr = [str(x) for x in rowint]
            rowout = '\t'.join(rowstr)
            print(rowout)
            # outf.write('\t'.join([iter[i], liklihood[i], perplexity[i]]))
            outf.write(rowout)

if __name__ == '__main__':
    main()
