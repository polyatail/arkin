#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/26/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean, log2
from common import load_barcodes, load_plate, map_bc_to_sample, multiplate_sorter
import subprocess
import sys
import os

def writetable(table, out_fname):
  taxa = list(set(sum([table[x].keys() for x in table], [])))

  with open(os.path.join(options.output_dir, out_fname), "w") as fp:
    fp.write("\t".join(["sample"] + taxa) + "\n")

    for sample in sorted(table.keys(), cmp=multiplate_sorter):
      row = [sample]

      for taxon in taxa:
        try:
          row.append(table[sample][taxon])
        except KeyError:
          row.append(0)

      fp.write("\t".join(map(str, row)) + "\n")

def normalize(table, std_org, log_ratio=False):
  new_table = {}

  for sample in table:
    new_table[sample] = {}

    for taxon in table[sample]:
      if std_org not in table[sample]:
        continue

      if log_ratio:
        new_table[sample][taxon] = log2(table[sample][taxon] / float(table[sample][std_org]))
      else:
        new_table[sample][taxon] = table[sample][taxon] / float(table[sample][std_org])

  return new_table

def percentages(table):
  new_table = {}

  for sample in table:
    new_table[sample] = {}

    for taxon in table[sample]:
      new_table[sample][taxon] = table[sample][taxon] * 100 / float(sum(table[sample].values()))

  return new_table

def mktable(iterable, barcode_key):
  table = dict([(x, {}) for x in barcode_key.values()])

  for l in iterable:
    l = l.strip().split("\t")

    barcode = "_".join(l[0].split("_")[:-1])

    try:
      sample = barcode_key[barcode]
    except KeyError:
      continue

    try:
      table[sample][l[1]] += 1
    except KeyError:
      table[sample][l[1]] = 1

  return table

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <assignments.txt> <barcodes.txt> <plate_layout.txt>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./tables_figs]",
                    default="./tables_figs",
                    help="write output files to this directory")

  parser.add_option("--std-org",
                    dest="std_org",
                    metavar="[10F2]",
                    default=False,
                    help="divide samples by counts for this organism")

  options, args = parser.parse_args(arguments)

  if len(args) <> 3:
    print "Error: Incorrect number of arguments"
    parser.print_help()
    sys.exit(1)

  if not os.path.isfile(args[0]):
    print "Error: Specified assignments file does not exist"
    parser.print_help()
    sys.exit(1)

  if not os.path.isfile(args[1]):
    print "Error: Specified barcode file does not exist"
    parser.print_help()
    sys.exit(1)

  if not os.path.isfile(args[2]):
    print "Error: Specified plate layout file does not exist"
    parser.print_help()
    sys.exit(1)

  if not options.std_org:
    print "Warning: Will not normalize without --std-org"

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specified path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

def main():
  parse_options(sys.argv[1:])

  # load barcodes
  fwd_bcs, rev_bcs = load_barcodes(args[1])

  # load plate layout
  plate = load_plate(args[2])

  # map barcodes to samples
  barcode_to_sample = map_bc_to_sample(plate, fwd_bcs, rev_bcs)

  # generate table
  sys.stderr.write("Tallying assigned reads...")
  table = mktable(open(args[0], "r"), barcode_to_sample)
  writetable(table, "tally.txt")

  sys.stderr.write("\n\nWriting table, percentages of total reads...")
  percentage_table = percentages(table)
  writetable(percentage_table, "tally.percentages.txt")

  if options.std_org:
    sys.stderr.write("\n\nWriting table, counts normalized to %s..." % options.std_org)
    normalized_table = normalize(table, options.std_org)
    writetable(normalized_table, "tally.normalized.txt")

    sys.stderr.write("\nWriting table, log2(counts normalized to %s)..." % options.std_org)
    log2_norm_table = normalize(table, options.std_org, log_ratio=True)
    writetable(log2_norm_table, "tally.normalized.log2.txt")

  sys.stderr.write("\n")

if __name__ == "__main__":
  main()
