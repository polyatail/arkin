#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/24/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean
from itertools import izip
from fastq import fast_fastq
from string import maketrans, translate
from common import load_barcodes, load_plate, map_bc_to_sample
import sys
import time
import os
import subprocess
import tempfile
import gzip
import zipfile

PATH_TO_PEAR = "pear"

hamming_dist = lambda a, b: sum([1 for x, y in zip(a, b) if x != y])
revcomp_trans = maketrans("ATGCN", "TACGN")
revcomp = lambda x: translate(x[::-1], revcomp_trans)

def fetch_unzip(zipfile_obj, filename):
  # accepts zipfile object and filename contained within, returns
  # temporary file containing raw reads

  gzipped = tempfile.NamedTemporaryFile(delete=False)
  gzipped.write(zipfile_obj.open(filename, "r").read())
  gzipped.flush()

  ungzipped = tempfile.NamedTemporaryFile(delete=False)
  ungzipped.write(gzip.GzipFile(gzipped.name).read())
  ungzipped.flush()

  os.unlink(gzipped.name)

  return ungzipped

def pear(fwd_fname, rev_fname, mem_size, num_threads):
  # accepts fwd and rev paired end reads, returns temporary file containing
  # presumably overlapping reads merged into single reads and percentage of
  # reads unable to be merged (high % indicates a problem)

  outfile = tempfile.NamedTemporaryFile(delete=False)

  pear = subprocess.Popen([PATH_TO_PEAR,
                           "-f", fwd_fname, "-r", rev_fname,
                           "-o", outfile.name,
                           "-y", mem_size, "-j", str(num_threads),
                           "-p", "0.01"],
                          stdout=subprocess.PIPE)

  while pear.poll() == None:
    time.sleep(0.1)

  if pear.returncode <> 0:
    raise ValueError("Pear exited with an error (%s)" % pear.returncode)

  raw_pear_log = []
  stats = {}

  for line in pear.stdout:
    raw_pear_log.append(line)
    line = line.strip()

    if line.endswith("%)"):
      perc = float(line.split(" ")[-1][1:-2])

      if line.startswith("Assembled reads"):
        stats["assembled_reads"] = perc
      elif line.startswith("Discarded reads"):
        stats["discarded_reads"] = perc
      elif line.startswith("Not assembled reads"):
        stats["unassembled_reads"] = perc

  return outfile, stats, "".join(raw_pear_log)

def find_pairs(zipfile_obj):
  # accepts zipfile object and returns tuples of paired-end reads

  filenames = zipfile_obj.namelist()
  used_filenames = []
  pairs = []

  for fname in filenames:
    mate = fname.split("_")

    if mate[-2] == "R1":
      mate[-2] = "R2"
    elif mate[-2] == "R2":
      mate[-2] = "R1"
    else:
      raise ValueError("Invalid filename format (%s)" % fname)

    mate = "_".join(mate)

    if mate in filenames:
      used_filenames.extend((fname, mate))
      pairs.append(tuple(sorted((fname, mate))))

  unused_filenames = set(filenames).difference(used_filenames)

  if unused_filenames:
    raise ValueError("Found %s samples without pairs" % len(unused_filenames))

  return list(set(pairs))

def match_bc(read, bcs, max_mismatch):
  # greedy barcode matching algorithm returns first barcode with a match
  # hamming distance <= max mismatches away

  for barcode in bcs:
    if hamming_dist(read[:len(barcode)], barcode) <= max_mismatch:
      # barcode match
      return barcode
  else:
    return False

def demultiplex(fastq_fwd, fwd_bcs, fastq_rev, rev_bcs, max_mismatch):
  # takes input reads, tries to assign them bins based on barcodes, then
  # returns trimmed reads and bin name

  fwd_match = match_bc(fastq_fwd.sequence, fwd_bcs, max_mismatch)

  if fastq_rev == None:
    # merged read
    rev_match = match_bc(revcomp(fastq_fwd.sequence), rev_bcs, max_mismatch)
  else:
    rev_match = match_bc(fastq_rev.sequence, rev_bcs, max_mismatch)

  if fwd_match == False or rev_match == False:
    # couldn't match to both fwd and rev barcodes
    return False
  else:
    # trim reads and quality scores
    if fastq_rev == None:
      # merged read
      fastq_fwd = fastq_fwd[len(fwd_match):-len(rev_match)]
    else:
      fastq_fwd = fastq_fwd[len(fwd_match):]
      fastq_rev = fastq_rev[len(rev_match):]

    return (fastq_fwd, fastq_rev, fwd_match, rev_match)

def qual_filter(read, min_qual, phred_offset, max_errors):
  votes = []

  if min_qual:
    if mean([ord(x) - phred_offset for x in read.quals]) < min_qual:
      votes.append(False)
    else:
      votes.append(True)
  else:
    votes.append(True)

  if max_errors:
    if sum([10 ** -((ord(x) - phred_offset)/10.0) for x in read.quals]) > max_errors:
      votes.append(False)
    else:
      votes.append(True)
  else:
    votes.append(True)

  return reduce(lambda x, y: x & y, votes)

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <-f/-r/-m/-z/-s filename> <barcodes.txt> <plate_layout.txt>",
                        version="%prog " + str(__version__))

  group1 = OptionGroup(parser, "Input Files")
  group2 = OptionGroup(parser, "Demultiplexing")
  group3 = OptionGroup(parser, "Read Merging")
  group4 = OptionGroup(parser, "Quality Filtering")

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./reads]",
                    default="./reads",
                    help="write output files to this directory")

  parser.add_option("-p",
                    dest="num_threads",
                    type="int",
                    metavar="[1]",
                    default=1,
                    help="number of threads used during analysis")

  parser.add_option("--mem",
                    dest="mem_size",
                    metavar="[1G]",
                    default="1G",
                    help="amount of memory to use during read merging")

  group1.add_option("-f",
                    dest="fwd_fname",
                    type="str",
                    default=False,
                    help="path to FASTQ containing forward reads")

  group1.add_option("-r",
                    dest="rev_fname",
                    type="str",
                    default=False,
                    help="path to FASTQ containing reverse reads")

  group1.add_option("-m",
                    dest="merged_fname",
                    type="str",
                    default=False,
                    help="path to FASTQ containing merged reads")

  group1.add_option("-z",
                    dest="zip_fname",
                    type="str",
                    default=False,
                    help="path to zipfile from BaseSpace")

  group1.add_option("-s",
                    dest="sample_name",
                    type="str",
                    default=False,
                    help="name of sample within BaseSpace zipfile")

  group2.add_option("--max-mismatch",
                    dest="max_mismatch",
                    type="int",
                    metavar="[2]",
                    default=2,
                    help="maximum allowed mismatches in barcodes")

  group2.add_option("--no-plate-layout",
                    dest="use_plate",
                    action="store_false",
                    metavar="",
                    default=True,
                    help="name reads by barcodes, not by sample names")

  group3.add_option("--merge",
                    dest="merge",
                    action="store_true",
                    default=False,
                    help="merge forward and reverse reads with PEAR")

  group3.add_option("--min-merged-perc",
                    dest="min_merged_perc",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="warn if fewer than this % of reads can be merged")

  group4.add_option("--phred_offset",
                    dest="phred_offset",
                    type="int",
                    metavar="[33]",
                    default=33,
                    help="phred quality offset (default 33 for illumina 1.8+)")

  group4.add_option("--mean-qual",
                    dest="min_qual",
                    type="int",
                    metavar="[25]",
                    default=False,
                    help="discard reads with mean quality below this value")

  group4.add_option("--max-errors",
                    dest="max_errors",
                    type="float",
                    metavar="[2]",
                    default=False,
                    help="discard reads with more than this many expected errors")

  group4.add_option("--min-qual-perc",
                    dest="min_qual_perc",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="warn if fewer than this % of reads pass quality filter")

  parser.add_option_group(group1)
  parser.add_option_group(group2)
  parser.add_option_group(group3)
  parser.add_option_group(group4)

  options, args = parser.parse_args(arguments)

  if (options.use_plate == True and len(args) <> 2) or \
     (options.use_plate == False and len(args) <> 1):
    print "Error: Incorrect number of arguments"
    parser.print_help()
    sys.exit(1)

  if not (0 <= options.min_merged_perc <= 100):
    print "Error: --min-merged-perc must be in [0,100]"
    parser.print_help()
    sys.exit(1)

  if not (0 <= options.min_qual_perc <= 100):
    print "Error: --min-qual-perc must be in [0,100]"
    parser.print_help()
    sys.exit(1)

  options.min_merged_perc /= 100.0
  options.min_qual_perc /= 100.0

  if options.merged_fname and \
     (options.zip_fname or options.sample_name) and \
     (options.fwd_fname or options.rev_fname):
    print "Error: Options -m and -z/-s and -f/-r are mutually exclusive"
    parser.print_help()
    sys.exit(1)
  elif not (options.merged_fname or options.zip_fname or \
            options.sample_name or options.fwd_fname or \
            options.rev_fname):
    print "Error: Must specify either -m, -z/-s, or -f/-r"
    parser.print_help()
    sys.exit(1)
  elif bool(options.fwd_fname) ^ bool(options.rev_fname):
    print "Error: When using -f or -r, both must be specified"
    parser.print_help()
    sys.exit(1)
  elif bool(options.zip_fname) ^ bool(options.sample_name):
    print "Error: When using -z or -s, both must be specified"
    parser.print_help()
    sys.exit(1)

  if options.merged_fname and options.merge:
    print "Error: Option --merge doesn't make sense with -m"
    parser.print_help()
    sys.exit(1)

  for fname in (options.zip_fname, options.fwd_fname, options.rev_fname, options.merged_fname):
    if fname != False and \
       not os.path.isfile(fname):
      print "Error: Specified FASTQ or zip file (%s) does not exist" % fname
      parser.print_help()
      sys.exit(1)

  if not os.path.isfile(args[0]):
    print "Error: Specified barcode file does not exist"
    parser.print_help()
    sys.exit(1)

  if options.use_plate and not os.path.isfile(args[1]):
    print "Error: Specified plate layout file does not exist"
    parser.print_help()
    sys.exit(1)

  if not (options.min_qual or options.max_errors):
    print "Warning: Specify --min-qual or --max-errors to enable quality filtering"

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specified path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

def main():
  parse_options(sys.argv[1:])

  # load barcodes
  fwd_bcs, rev_bcs = load_barcodes(args[0])

  if options.use_plate:
    plate = load_plate(args[1])
    barcode_to_sample = map_bc_to_sample(plate, fwd_bcs, rev_bcs)

  if options.zip_fname:
    # open zipfile
    miseq_zip = zipfile.ZipFile(options.zip_fname)

    # find pairs
    pairs = dict([("_".join(x.split("/")[-1].split(".")[0].split("_")[:-2]), (x, y)) for x, y in find_pairs(miseq_zip)])

    if options.sample_name not in pairs:
      raise ValueError("Could not find %s in %s!" % (options.sample_name, options.zip_fname))

    # extract read files
    sys.stderr.write("Extracting reads from zipfile...\n")

    fwd = fetch_unzip(miseq_zip, pairs[options.sample_name][0])
    rev = fetch_unzip(miseq_zip, pairs[options.sample_name][1])
  elif options.fwd_fname and options.rev_fname:
    fwd = open(options.fwd_fname, "r")
    rev = open(options.rev_fname, "r")

  barcode_to_count = {}

  if options.merge or options.merged_fname:
    with open(os.path.join(options.output_dir, "merged_reads.assigned.fastq"), "w") as assigned, \
         open(os.path.join(options.output_dir, "merged_reads.unassigned.fastq"), "w") as unassigned:
      read_length_bins = {}

      total_reads = 0.0
      quality_reads = 0.0

      if options.merge:
        # run pear to merge reads
        sys.stderr.write("Merging reads...\n")
        merged, stats, raw_pear_log = pear(fwd.name, rev.name, options.mem_size, options.num_threads)

        open(os.path.join(options.output_dir, "pear.log"), "w").write(raw_pear_log)

        if stats["assembled_reads"] < options.min_merged_perc:
          sys.stderr.write("  Warning: only %.02f%% of reads assembled\n" % stats["assembled_reads"])

        merged_fastq = fast_fastq(open("%s.assembled.fastq" % merged.name, "r"))
      elif options.merged_fname:
        merged_fastq = fast_fastq(open(options.merged_fname, "r"))

      # read through fastq
      sys.stderr.write("Filtering and demultiplexing reads...\n")

      for merged_read in merged_fastq:
        total_reads += 1

        # check that read passes quality filter
        if not qual_filter(merged_read, options.min_qual, options.phred_offset, options.max_errors):
          continue

        # keep track of stats
        quality_reads += 1

        try:
          read_length_bins[len(merged_read.sequence)] += 1
        except KeyError:
          read_length_bins[len(merged_read.sequence)] = 1

        # demultiplex
        dm_out = demultiplex(merged_read, fwd_bcs.values(), None, rev_bcs.values(), options.max_mismatch)

        if dm_out == False:
          # strip pair info from read and write to unassigned file
          merged_read.id = merged_read.id.split(" ")[0]
          unassigned.write(merged_read.raw())
        else:
          # rename read to barcodes used and write to assigned file
          trimmed_read, _, f_bc, r_bc = dm_out

          bc_name = "%s_%s" % (f_bc, r_bc)

          try:
            barcode_to_count[bc_name] += 1
          except KeyError:
            barcode_to_count[bc_name] = 1

          # we want to rename to sample names, and a barcode was found to match
          # something in barcodes.txt, but that particular barcode is not in use
          # in our plate layout. in this case, ignore this read
          if options.use_plate and bc_name not in barcode_to_sample:
            # strip pair info from read and write to unassigned file
            merged_read.id = merged_read.id.split(" ")[0]
            unassigned.write(merged_read.raw())

            continue

          if options.use_plate:
            trimmed_read.id = "%s_%s" % (barcode_to_sample[bc_name], barcode_to_count[bc_name])
          else:
            trimmed_read.id = "%s_%s" % (bc_name, barcode_to_count[bc_name])

          assigned.write(trimmed_read.raw())

      if total_reads > 0:
        if quality_reads / total_reads < options.min_qual_perc:
          sys.stderr.write("  Warning: only %.02f%% of reads passed quality filter\n" % (quality_reads * 100 / total_reads))

      sys.stderr.write("\nSummary")
      sys.stderr.write("\n  Total reads:       %d" % total_reads)
      sys.stderr.write("\n  Quality reads:     %d" % quality_reads)
      sys.stderr.write("\n  Min read length:   %d" % min(read_length_bins.keys()))
      sys.stderr.write("\n  Mean read length:  %d" % mean(read_length_bins.keys()))
      sys.stderr.write("\n  Max read length:   %d\n" % max(read_length_bins.keys()))
      sys.stderr.write("\n  Assigned reads:    %d" % sum(barcode_to_count.values()))
      sys.stderr.write("\n  Unassigned reads:  %d" % (quality_reads - sum(barcode_to_count.values())))
      sys.stderr.write("\n  Avg reads/barcode: %d\n" % mean(barcode_to_count.values()))

      with open(os.path.join(options.output_dir, "read_lengths.log"), "w") as fp:
        for read_length, count in sorted(read_length_bins.items(), key=lambda x: x[0], reverse=True):
          fp.write("%s\t%s\n" % (read_length, count))

      if options.merge:
        # remove temporary files
        os.unlink(merged.name)
        os.unlink("%s.assembled.fastq" % merged.name)
        os.unlink("%s.discarded.fastq" % merged.name)
        os.unlink("%s.unassembled.forward.fastq" % merged.name)
        os.unlink("%s.unassembled.reverse.fastq" % merged.name)
  else:
    with open(os.path.join(options.output_dir, "fwd_reads.assigned.fastq"), "w") as f_assigned, \
         open(os.path.join(options.output_dir, "rev_reads.assigned.fastq"), "w") as r_assigned, \
         open(os.path.join(options.output_dir, "fwd_reads.unassigned.fastq"), "w") as f_unassigned, \
         open(os.path.join(options.output_dir, "rev_reads.unassigned.fastq"), "w") as r_unassigned:
      total_reads = 0.0
      quality_reads = 0.0

      # read through foward and reverse fastq simultaneously
      sys.stderr.write("Filtering and demultiplexing reads...\n")

      for f_read, r_read in izip(fast_fastq(open(fwd.name, "r")), fast_fastq(open(rev.name, "r"))):
        total_reads += 1

        # check that read passes quality filter
        if qual_filter(f_read, options.min_qual, options.phred_offset, options.max_errors) == False or \
           qual_filter(r_read, options.min_qual, options.phred_offset, options.max_errors) == False:
          continue

        quality_reads += 1

        # demultiplex
        dm_out = demultiplex(f_read, fwd_bcs.values(), r_read, rev_bcs.values(), options.max_mismatch)

        if dm_out == False:
          # strip pair info from read and write to unassigned file
          f_read.id = f_read.id.split(" ")[0]
          f_unassigned.write(f_read.raw())

          r_read.id = r_read.id.split(" ")[0]
          r_unassigned.write(r_read.raw())
        else:
          # rename reads to barcodes used and write to assigned file
          trimmed_f_read, trimmed_r_read, f_bc, r_bc = dm_out

          bc_name = "%s_%s" % (f_bc, r_bc)

          try:
            barcode_to_count[bc_name] += 1
          except KeyError:
            barcode_to_count[bc_name] = 1

          # we want to rename to sample names, and a barcode was found to match
          # something in barcodes.txt, but that particular barcode is not in use
          # in our plate layout. in this case, ignore this read
          if options.use_plate and bc_name not in barcode_to_sample:
            # strip pair info from read and write to unassigned file
            f_read.id = f_read.id.split(" ")[0]
            f_unassigned.write(f_read.raw())

            r_read.id = r_read.id.split(" ")[0]
            r_unassigned.write(r_read.raw())

            continue

          if options.use_plate:
            trimmed_f_read.id = "%s_%s" % (barcode_to_sample[bc_name], barcode_to_count[bc_name])
          else:
            trimmed_f_read.id = "%s_%s" % (bc_name, barcode_to_count[bc_name])

          f_assigned.write(trimmed_f_read.raw())

          if options.use_plate:
            trimmed_r_read.id = "%s_%s" % (barcode_to_sample[bc_name], barcode_to_count[bc_name])
          else:
            trimmed_r_read.id = "%s_%s" % (bc_name, barcode_to_count[bc_name])

          r_assigned.write(trimmed_r_read.raw())

      if total_reads > 0:
        if quality_reads / total_reads < options.min_qual_perc:
          sys.stderr.write("  Warning: only %.02f%% of reads passed quality filter\n" % (quality_reads * 100 / total_reads))

      sys.stderr.write("\nSummary")
      sys.stderr.write("\n  Total pairs:       %d" % total_reads)
      sys.stderr.write("\n  Quality pairs:     %d" % quality_reads)
      sys.stderr.write("\n  Assigned pairs:    %d" % sum(barcode_to_count.values()))
      sys.stderr.write("\n  Unassigned pairs:  %d" % (quality_reads - sum(barcode_to_count.values())))
      sys.stderr.write("\n  Avg pairs/barcode: %d\n" % int(mean(barcode_to_count.values())))

  sorted_barcode_to_count = sorted(barcode_to_count.items(), key=lambda x: x[1], reverse=True)

  with open(os.path.join(options.output_dir, "barcode_to_count.log"), "w") as fp:
    if options.use_plate:
      for barcode, count in sorted_barcode_to_count:
        try:
          sample_name = barcode_to_sample[barcode]
        except KeyError:
          sample_name = "BARCODE_NOT_IN_PLATE_LAYOUT"

        fp.write("%s\t%s\t%s\n" % (barcode, sample_name, count))
    else:
      for barcode, count in sorted_barcode_to_count:
        fp.write("%s\t%s\n" % (barcode, count))

  # remove temporary files
  if options.zip_fname:
    os.unlink(fwd.name)
    os.unlink(rev.name)

if __name__ == "__main__":
  main()
