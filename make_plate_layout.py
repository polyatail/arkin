#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "3/5/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
import sys
import os

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <samplevols.txt>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./data]",
                    default="./data",
                    help="write output files to this directory")

  options, args = parser.parse_args(arguments)

  if len(args) <> 1:
    print "Error: Incorrect number of arguments"
    parser.print_help()
    sys.exit(1)

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

def mkplate(iterable):
  plate = {}
  unknowns = {}
  barcodes = {}

  new_plate = False
  first_cols = False
  first_rows = False

  cols = []
  rows = []
 
  for l_num, l in enumerate(iterable):
    if l.startswith("#"):
      if l.strip().lower() == "#newplate":
        new_plate = True
      else:
        continue
    else:
      l_split = l.strip().split("\t")

      if new_plate:
        contains = l_split[0]

        if contains.lower().endswith("barcode"):
          contains_split = contains.lower().split("_")

          if contains_split[0] == "fwd":
            contains = "forward"
          elif contains_split[0] == "rev":
            contains = "reverse"
          else:
            raise ValueError("Line %s: Invalid barcode type specified: %s" % (l_num, contains))
        elif contains.lower() == "media":
          contains = "media"
        elif contains.lower() == "unknowns":
          contains = "unknowns"

        if first_cols:
          if len(first_cols) != len(l_split[1:]):
            raise ValueError("Line %s: Plate header doesn't have %s columns" % (l_num, len(first_cols)))
          elif set(l_split[1:]).symmetric_difference(first_cols):
            raise ValueError("Line %s: Plate has different column names from last plate" % l_num)
        else:
          first_cols = l_split[1:]

        if rows and not first_rows:
          first_rows = rows

        cols = l_split[1:]
        new_plate = False
      else:
        if len(l_split[1:]) != len(cols):
          raise ValueError("Line %s: Row doesn't have %s columns" % (l_num, len(cols)))

        if first_rows:
          if l_split[0] not in first_rows:
            raise ValueError("Line %s: Row %s not present in first plate" % (l_num, l_split[0]))
        else:
          rows.append(l_split[0])

        well_vals = [("%s%s" % (l_split[0], x), y) for x, y in zip(cols, l_split[1:])]

        if contains in ("forward", "reverse"):
          for well, val in well_vals:
            try:
              barcodes[well][contains] = val
            except KeyError:
              barcodes[well] = {contains: val}
        elif contains == "unknowns":
          for well, val in well_vals:
            try:
              unknowns[well][contains] = bool(val)
            except KeyError:
              unknowns[well] = {contains: bool(val)}
        else:
          for well, val in well_vals:
            try:
              plate[well][contains] = float(val)
            except KeyError:
              plate[well] = {contains: float(val)}

  return plate, barcodes, unknowns

def writeplate(plate, barcodes, unknowns, fname):
  with open(os.path.join(options.output_dir, fname), "w") as fp:
    fp.write("\t".join(["#sample", "fwd_barcode", "rev_barcode", "type", "expected", "description"]) + "\n")

    for well in sorted(plate.keys()):
      row = [well, barcodes[well]["forward"], barcodes[well]["reverse"]]

      # what type of row is this?
      if well in unknowns[well]:
        if sum(plate[well].values()) > 0:
          raise ValueError("Well %s: Listed as unknown but has non-zero known volumes" % well)

        row.append("Unknown_Std1")
      elif sum(plate[well].values()) == 0:
        row.append("Empty")
      else:
        row.append("Standard_Std1")

        # figure out what we expect to see
        denom_no_media = float(sum([plate[well][x] for x in plate[well] if x != "media"]))
        row.append(",".join(["%s:%s" % (x, y * 100 / denom_no_media) for x, y in plate[well].items() if x != "media" and y > 0]))

      fp.write("\t".join(row) + "\n")

def main():
  parse_options(sys.argv[1:])

  plate, barcodes, unknowns = mkplate(open(args[0], "r"))
  writeplate(plate, barcodes, unknowns, "plate_layout.txt")

if __name__ == "__main__":
  main()
