#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/26/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean, log2
import subprocess
import sys
import os

PATH_TO_SHUF = "shuf"

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

def map_bc_to_sample(plate, fwd_bcs, rev_bcs):
  barcode_to_sample = {}

  for sample in plate:
    barcode = "%s_%s" % (fwd_bcs[plate[sample]["fwd_barcode"]],
                         rev_bcs[plate[sample]["rev_barcode"]])

    if barcode in barcode_to_sample:
      raise ValueError("Sample %s shares barcode (%s) with %s" % (sample, barcode, barcode_to_sample[barcode]))

    barcode_to_sample[barcode] = sample

  return barcode_to_sample

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
      expected = [x.split(":") for x in plate[sample]["expected"].split(",")]

      try:
        data["standards"][std_name]
      except KeyError:
        data["standards"][std_name] = {}

      for orgname, expval in expected:
        expval = float(expval)

        if log_ratio:
          try:
            actual = log2(table[sample][orgname] / float(table[sample][options.std_org]))
          except KeyError:
            #sys.stderr.write("Warning: %s: Expected %s = %s%% of reads, but found 0\n" % (sample, orgname, expval))
            continue

          expval = log2(expval / (100 - sum([float(x[1]) for x in expected])))
        else:
          try:
            actual = table[sample][orgname]
          except KeyError:
            #sys.stderr.write("Warning: %s: Expected %s = %s%% of reads, but found 0\n" % (sample, orgname, expval))
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
                    metavar="[./]",
                    default="./",
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
    sys.stderr.write("Tallying assigned reads...\n")
    table = mktable(open(args[0], "r"), barcode_to_sample)
    writetable(table, "tally.txt")

    sys.stderr.write("\nWriting table, percentages of total reads...\n")
    percentage_table = percentages(table)
    writetable(percentage_table, "tally.percentages.txt")

    sys.stderr.write("Crunching percentage data...\n")
    perc_crunch = crunch(percentage_table, plate)

    if options.std_org:
      sys.stderr.write("\nWriting table, counts normalized to %s...\n" % options.std_org)
      normalized_table = normalize(table, options.std_org)
      writetable(normalized_table, "tally.normalized.txt")

      sys.stderr.write("Writing table, log2(counts normalized to %s)...\n" % options.std_org)
      log2_norm_table = normalize(table, options.std_org, log_ratio=True)
      writetable(log2_norm_table, "tally.normalized.log2.txt")

      sys.stderr.write("Crunching log2 ratio data...\n")
      log_crunch = crunch(percentage_table, plate, log_ratio=True)

if __name__ == "__main__":
  main()
