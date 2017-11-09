#! /usr/bin/env python

headerhelp = \
'''
View and analyze NMR spectra.

Details on actions:

viewhsqc    Show a single HSQC spectrum
dualhsqc    Show overlay of two HSQC spectra
peakhsqc    Find and edit peaks in HSQC spectrum
titrhsqc    Manually process set of NMR titrations
autocros    Auto-process a set of spectra looking for shifts
peakcros    Match a set of peaks looking for shifts

'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter, REMAINDER
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('-a', '--action', action='append',
                    choices =   [   
                                    'viewhsqc',
                                    'dualhsqc',
                                    'peakhsqc',
                                    'titrhsqc',
                                    'autocros',
                                    'peakcros',
                                    'peakmatch',
                                    'peakmxy',
                                    'dualpdfs',
                                ],
                    default = [],
                    metavar = '', help='Action to perform')
parser.add_argument('-i', '--input_file', help='Input data file (HSQC in .nv format)')
parser.add_argument('--peakfile', default="", help='Peak data file')
parser.add_argument('--auxpeakfile', help='Auxillary peak data file')
parser.add_argument('-n', '--num-peaks', default=50, type=int, help='Number of peaks to detect')
parser.add_argument('--bright',  type=int, help='Right edge of the box for peak search')
parser.add_argument('--bleft',  type=int, help='Left edge of the box for peak search')
parser.add_argument('--btop',  type=int, help='Top edge of the box for peak search')
parser.add_argument('--bbottom',  type=int, help='Bottom of the box for peak search')
parser.add_argument('--aponum',  type=int, help='apo-protein spectrum number')
parser.add_argument('--holonums', help='holo-protein spectrum number range, as python expression')
parser.add_argument('--controls', help='Control spectra number range, as python expression')
parser.add_argument('--folder', default='', help='FOlder with nv format HSQC numbered spectra files')
parser.add_argument('--vcutoff', default=0.1, type=float, help="Relative amplitude cutoff in oeak validation over multiple controls")
parser.add_argument('--dcutoff', default=0.1, type=float, help="Peak matching cutoff, ppm")
parser.add_argument('--scutoff', default=2, type=float, help="Next peak separation cutoff, relative units")
parser.add_argument('--zcutoff', default=2.0, type=float, help="Peak matching cutoff, Z-score")
parser.add_argument('--print-peakmatch', action='store_true')
parser.add_argument('--print-full', action='store_true')
parser.add_argument('--pdfile', default="epynmr.pdf", help='Output PDF file name')

parser.add_argument('xargs', nargs=REMAINDER)

args = parser.parse_args()

import epynmractions

for action in args.action:
    epynmractions.__getattribute__(action)(args)

