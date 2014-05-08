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

  parser = OptionParser(usage="%prog [options] <-v samplevols.txt|-b barcodes.txt -m plate_num [-m ...]>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./data]",
                    default="./data",
                    help="write output files to this directory")

  parser.add_option("-v",
                    dest="samplevols",
                    metavar="[samplevols.txt]",
                    default=False,
                    help="expected volume mode, use these sample volumes")

  parser.add_option("-b",
                    dest="barcodes",
                    metavar="[barcodes.txt]",
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

  if options.barcodes and not options.multi_plates:
    print "Error: -m must be specified with -b"
    parser.print_help()
    sys.exit(1)

  if options.samplevols and (options.barcodes or options.multi_plates):
    print "Error: -v cannot be used with -m or -b"
    parser.print_help()
    sys.exit(1)

  if not (options.samplevols or options.barcodes):
    print "Error: Must specify one of -b or -v"
    parser.print_help()
    sys.exit(1)

  if options.multi_plates:
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
  all_plates = {}

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

      plate = {}
      barcodes = {}
      unknowns = {}

      for row_num, row_letter in enumerate("ABCDEFGH"):
        for col_num, col_letter in enumerate(range(1, 13)):
          well = "%s_%s%s" % (l[0], row_letter, col_letter)

          unknowns[well] = {"unknowns": True}
          barcodes[well] = {"forward": fwd_barcodes[l[1]][row_num],
                            "reverse": rev_barcodes[l[2]][col_num]}
          plate[well] = {}

      all_plates[l[0]] = (plate, barcodes, unknowns)

  return all_plates

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

def writeplate(plate, barcodes, unknowns, fname, append = False):
  writemode = "a" if append else "w"

  with open(os.path.join(options.output_dir, fname), writemode) as fp:
    if not append:
      fp.write("\t".join(["#sample", "fwd_barcode", "rev_barcode", "type", "expected", "description"]) + "\n")

    for well in sorted(plate.keys()):
      row = [well, barcodes[well]["forward"], barcodes[well]["reverse"]]

      # what type of row is this?
      if well in unknowns:
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

  if options.samplevols:
    plate, barcodes, unknowns = mkplate(open(options.samplevols, "r"))
    writeplate(plate, barcodes, unknowns, "plate_layout.txt")
  elif options.barcodes:
    barcodes = load_multi_barcode(open(options.barcodes, "r"))

    append = False

    for plate_num in options.multi_plates:
      writeplate(barcodes[plate_num][0], barcodes[plate_num][1],
                 barcodes[plate_num][2], "plate_layout.txt", append)
      append = True

if __name__ == "__main__":
  main()
