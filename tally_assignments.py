# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/26/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean
import sys
import os

def load_barcodes(bc_file):
  fwd_bcs = {}
  rev_bcs = {}

  header = False

  for l in open(bc_file, "r"):
    if l.startswith("#"):
      if header:
        continue
      elif l.startswith("#name"):
        header = l[1:].strip().split("\t")
    else:
      l = dict(zip(header, l.strip().split("\t")))

      if l["orientation"] == "forward":
        fwd_bcs[l["name"]] = l["barcode"]

      if l["orientation"] == "reverse":
        rev_bcs[l["name"]] = l["barcode"]

  return (fwd_bcs, rev_bcs)

def load_plate(plate_file):
  sample_data = {}
  header = False

  for l in open(plate_file, "r"):
    if l.startswith("#"):
      if header:
        continue
      elif l.startswith("#sample"):
        header = l[1:].strip().split("\t")
    else:
      l = dict(zip(header, l.strip().split("\t")))

      sample_data[l["sample"]] = l

  return sample_data

def crunch(plate, fwd_bcs, rev_bcs):
  barcode_to_sample = {}

  for sample in plate:
    barcode = "%s_%s" % (fwd_bcs[plate[sample]["fwd_barcode"]],
                         rev_bcs[plate[sample]["rev_barcode"]])

    if barcode in barcode_to_sample:
      raise ValueError("Sample %s shares barcode (%s) with %s" % (sample, barcode, barcode_to_sample[barcode]))

    barcode_to_sample[barcode] = sample

  return barcode_to_sample

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <assignments.txt> <barcodes.txt> <plate_layout.txt>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./miseq_out]",
                    default="./miseq_out",
                    help="write output files to this directory")

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

def main():
  parse_options(sys.argv[1:])

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specified path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

  # load barcodes
  fwd_bcs, rev_bcs = load_barcodes(args[1])

  # load plate layout
  plate = load_plate(args[2])

  # map barcodes to samples
  barcode_to_sample = crunch(plate, fwd_bcs, rev_bcs)

  table = dict([(x, {}) for x in barcode_to_sample.values()])
  taxa = set()
 
  for l in open(args[0], "r"):
    l = l.strip().split("\t")

    barcode = "_".join(l[0].split("_")[:-1])

    try:
      sample = barcode_to_sample[barcode]
    except KeyError:
      continue

    taxa.add(l[1])

    try:
      table[sample][l[1]] += 1
    except KeyError:
      table[sample][l[1]] = 1

  taxa_order = list(taxa)

  print "\t".join(["sample"] + taxa_order)

  for sample in sorted(table.keys()):
    row = [sample]

    for taxon in taxa_order:
      try:
        row.append(table[sample][taxon])
      except KeyError:
        row.append(0)

    print "\t".join(map(str, row))

if __name__ == "__main__":
  main()
