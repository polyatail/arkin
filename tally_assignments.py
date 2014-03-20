#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/26/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean, log2
from common import load_barcodes, load_plate, map_bc_to_sample
import subprocess
import sys
import os

PATH_TO_SHUF = "shuf"

def writetable(table, out_fname):
  taxa = list(set(sum([table[x].keys() for x in table], [])))

  with open(os.path.join(options.output_dir, out_fname), "w") as fp:
    fp.write("\t".join(["sample"] + taxa) + "\n")

    for sample in sorted(table.keys()):
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

def crunch(table, plate, log_ratio=False):
  # contains actual vs expected values
  data = {"standards": {},
          "unknowns": {}}

  for sample in plate:
    if plate[sample]["type"].lower() == "empty":
      #TODO: add some checks that the well actually is empty
      continue
    elif plate[sample]["type"].lower().startswith("unknown"):
      std_name = plate[sample]["type"].split("_")[1]
      #TODO: keep track of raw values
      continue
    elif plate[sample]["type"].lower().startswith("standard"):
      std_name = plate[sample]["type"].split("_")[1]
      expected = dict([(x.split(":")[0], float(x.split(":")[1])) for x in plate[sample]["expected"].split(",")])

      try:
        data["standards"][std_name]
      except KeyError:
        data["standards"][std_name] = {}

      for orgname, expval in expected.items():
        if log_ratio:
          try:
            actual = log2(table[sample][orgname] / float(table[sample][options.std_org]))
          except KeyError:
            sys.stderr.write("Warning: %s: Expected %s = %s%% of reads, but found 0\n" % (sample, orgname, expval))
            continue

          # if an expected value for the standard is there, use it
          if options.std_org in expected:
            denom = expected[options.std_org]
          # otherwise assume it occupies the rest of the sample
          else:
            denom = (100 - sum([float(x[1]) for x in expected.items()]))

          expval = log2(expval / denom)
        else:
          try:
            actual = table[sample][orgname]
          except KeyError:
            sys.stderr.write("Warning: %s: Expected %s = %s%% of reads, but found 0\n" % (sample, orgname, expval))
            continue

        try:
          data["standards"][std_name][orgname]
        except KeyError:
          data["standards"][std_name][orgname] = {}

        try:
          data["standards"][std_name][orgname][expval]
        except KeyError:
          data["standards"][std_name][orgname][expval] = []

        data["standards"][std_name][orgname][expval].append(actual)

  #TODO: do regression (linear??) and quantitate unknowns

  # write data for plotting
  suffix = ".log2.txt" if log_ratio else ".percentages.txt"

  with open(os.path.join(options.output_dir, "plot_titles" + suffix), "w") as titles_fp, \
       open(os.path.join(options.output_dir, "plot_data" + suffix), "w") as data_fp:
    for std_name in data["standards"]:
      for org_name in data["standards"][std_name]:
        if org_name == options.std_org:
          continue

        header, row = [], []

        for expval, actval in data["standards"][std_name][org_name].items():
          header.extend([expval] * len(actval))
          row.extend(actval)

        titles_fp.write("Standard '%s': %s\n" % (std_name, org_name))
        data_fp.write("%s\n%s\n" % ("\t".join(map(str, header)), "\t".join(map(str, row))))

  return data

def rms_error(crunched):
  errors = []

  for std_name in crunched["standards"].keys():
    for org_name in crunched["standards"][std_name]:
      for expval in crunched["standards"][std_name][org_name]:
        errors.extend([expval - x for x in crunched["standards"][std_name][org_name][expval]])

  return rms(errors)

def rms(vec):
  return (sum([x ** 2 for x in vec]) / float(len(vec))) ** 0.5

def subsample(fname, size):
  shuf = subprocess.Popen([PATH_TO_SHUF],
                          stdin=open(fname, "r"),
                          stdout=subprocess.PIPE)

  sample = []

  for l in shuf.stdout:
    sample.append(l)

    if len(sample) == size:
      break

  shuf.kill()

  return sample

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

  parser.add_option("--rarefaction",
                    dest="rarefaction",
                    action="store_true",
                    default=False,
                    help="skip tally and do rarefaction analysis")

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

  if options.rarefaction:
    sys.stderr.write("Performing rarefaction analysis...\n")
    def _countlines(fname):
      for i, _ in enumerate(open(fname)):
        pass

      return i

    maxlines = _countlines(args[0])

    sys.stderr.write("  Total assigned reads: %s\n" % maxlines)

    # 100 to 100M or maxlines reads
    with open(os.path.join(options.output_dir, "rarefaction.txt"), "w") as fp:
      fp.write("reads_sampled\trms_error\n")

      for exp in range(2, 9):
        i = 10 ** exp

        if i > maxlines:
          i = maxlines
          trials = 1
        else:
          trials = 10

        for j in range(trials):
          sys.stderr.write("\r  Sampling %s reads trial %s/%s" % (i, j+1, trials))
          ss_reads = subsample(args[0], i)
          ss_table = mktable(ss_reads, barcode_to_sample)
          ss_perc_table = percentages(ss_table)
          ss_log_crunch = crunch(ss_perc_table, plate, log_ratio=True)

          fp.write("%s\t%s\n" % (i, rms_error(ss_log_crunch)))

        sys.stderr.write("\n")

        if i == maxlines:
          break
  else:
    # generate table
    sys.stderr.write("Tallying assigned reads...")
    table = mktable(open(args[0], "r"), barcode_to_sample)
    writetable(table, "tally.txt")

    sys.stderr.write("\n\nWriting table, percentages of total reads...")
    percentage_table = percentages(table)
    writetable(percentage_table, "tally.percentages.txt")

    sys.stderr.write("\nCrunching percentage data...")
    perc_crunch = crunch(percentage_table, plate)

    if options.std_org:
      sys.stderr.write("\n\nWriting table, counts normalized to %s..." % options.std_org)
      normalized_table = normalize(table, options.std_org)
      writetable(normalized_table, "tally.normalized.txt")

      sys.stderr.write("\nWriting table, log2(counts normalized to %s)..." % options.std_org)
      log2_norm_table = normalize(table, options.std_org, log_ratio=True)
      writetable(log2_norm_table, "tally.normalized.log2.txt")

      sys.stderr.write("\nCrunching log2 ratio data...")
      log_crunch = crunch(percentage_table, plate, log_ratio=True)

  sys.stderr.write("\n")

if __name__ == "__main__":
  main()
