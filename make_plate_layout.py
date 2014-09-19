#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "3/5/2014"
__version__ = 1.0

from common import multiplate_sorter
from optparse import OptionParser, OptionGroup
import sys
import os

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] -b arrayed_barcodes.txt -m plate_num [-m ...]",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./data]",
                    default="./data",
                    help="write output files to this directory")

  parser.add_option("-b",
                    dest="barcodes",
                    metavar="[arrayed_barcodes.txt]",
                    default=False,
                    help="multi-plate mode, use this barcode layout")

  parser.add_option("-m",
                    dest="multi_plates",
                    action="append",
                    metavar="[plate_num]",
                    default=[],
                    help="multi-plate mode, make this plate number")

  options, args = parser.parse_args(arguments)

  if len(args) <> 0:
    print "Error: Incorrect number of arguments"
    parser.print_help()
    sys.exit(1)

  if not (options.barcodes and options.multi_plates):
    print "Error: both -b and -m must be specified"
    parser.print_help()
    sys.exit(1)

  options.multi_plates = map(int, options.multi_plates)

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specified path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

  if os.path.isfile(os.path.join(options.output_dir, "plate_layout.txt")):
    print "Error: File plate_layout.txt already exists in output directory"
    parser.print_help()
    sys.exit(1)

def load_multi_barcode(iterable):
  rev_barcodes = {}
  fwd_barcodes = {}
  plates = {}

  action = False
  counter = False

  for l in iterable:
    l = l.strip()

    if l.startswith("#reverse"):
      action = "reverse"
      counter = 1
      continue
    elif l.startswith("#forward"):
      action = "forward"
      counter = 1
      continue
    elif l.startswith("#plate"):
      action = "plate"
      continue

    if action == "reverse":
      rev_barcodes[counter] = ["R%02d" % int(x) for x in l.split()]
      counter += 1
    elif action == "forward":
      fwd_barcodes[counter] = ["F%02d" % int(x) for x in l.split()]
      counter += 1
    elif action == "plate":
      l = map(int, l.split())

      plates[l[0]] = {}

      for row_num, row_letter in enumerate("ABCDEFGH"):
        for col_num, col_letter in enumerate(range(1, 13)):
          well = "%s_%s%s" % (l[0], row_letter, col_letter)

          plates[l[0]][well] = {"forward": fwd_barcodes[l[1]][row_num],
                                "reverse": rev_barcodes[l[2]][col_num]}

  return plates

def writeplate(plate, fname, append = False):
  writemode = "a" if append else "w"

  with open(os.path.join(options.output_dir, fname), writemode) as fp:
    if not append:
      fp.write("\t".join(["#sample", "fwd_barcode", "rev_barcode"]) + "\n")

    for well in sorted(plate.keys(), cmp=multiplate_sorter):
      row = [well, plate[well]["forward"], plate[well]["reverse"]]
      fp.write("\t".join(row) + "\n")

def main():
  parse_options(sys.argv[1:])

  barcodes = load_multi_barcode(open(options.barcodes, "r"))

  append = False

  for plate_num in options.multi_plates:
    writeplate(barcodes[plate_num], "plate_layout.txt", append)
    append = True

if __name__ == "__main__":
  main()
