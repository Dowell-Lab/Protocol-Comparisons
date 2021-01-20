import argparse
import main

parser = argparse.ArgumentParser(description = 'Running DWT on nascent sequencing data.')
requiredNamed = parser.add_argument_group('required named arguments')

parser.add_argument('-s', '--samples', dest="samp", help = 'file with paths to bedgraphs', metavar="FILE")
parser.add_argument('-o', '--outdirectory', dest="outdir", help = 'directory for output', metavar="FILE")
parser.add_argument('-a', '--annotations', dest="anno", help = 'annotation file (.bed, gff, gff3, gtf)', metavar="FILE")
parser.add_argument('-w', '--wavelet', dest="wave", help = 'wavelet function to scan', metavar="STR")
parser.add_argument('-p', '--plot', dest="plt", help = "plot PCA plots for all genes", action='store_true')

args = parser.parse_args()


main.run(args.samp, args.outdir, args.anno, args.wave, args.plt)

