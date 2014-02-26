# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/24/2014"
__version__ = 1.0

from optparse import OptionParser, OptionGroup
from numpy import mean
from itertools import izip
import sys
import time
import os
import subprocess
import tempfile
import gzip
import zipfile

PATH_TO_PEAR = "pear"

hamming_dist = lambda a, b: sum([1 for x, y in zip(a, b) if x != y])

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

def pear(fwd_fname, rev_fname):
  # accepts fwd and rev paired end reads, returns temporary file containing
  # presumably overlapping reads merged into single reads and percentage of
  # reads unable to be merged (high % indicates a problem)

  outfile = tempfile.NamedTemporaryFile(delete=False)

  pear = subprocess.Popen([PATH_TO_PEAR,
                           "-f", fwd_fname, "-r", rev_fname,
                           "-o", outfile.name,
                           "-y", options.mem_size, "-j", str(options.num_threads),
                           "-p", "0.01"],
                          stdout=subprocess.PIPE)

  while pear.poll() == None:
    time.sleep(0.1)

  if pear.returncode <> 0:
    raise ValueError("Pear exited with an error (%s)" % pear.returncode)

  stats = {}

  for line in pear.stdout:
    line = line.strip()

    if line.endswith("%)"):
      perc = float(line.split(" ")[-1][1:-2])

      if line.startswith("Assembled reads"):
        stats["assembled_reads"] = perc
      elif line.startswith("Discarded reads"):
        stats["discarded_reads"] = perc
      elif line.startswith("Not assembled reads"):
        stats["unassembled_reads"] = perc

  return outfile, stats

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

def match_bc(read, bcs, max_spacer = 7, max_hamming = 2):
  # greedy barcode matching algorithm returns first barcode with a match
  # of less than or equal to max_hamming distance away

  for i in range(0, max_spacer + 1):
    for barcode in bcs:
      if hamming_dist(read[i:i+len(barcode)], barcode) <= max_hamming:
        # barcode match
        return (i, barcode)
  else:
    return False

def demultiplex(fastq_fwd, fwd_bcs, rev_bcs, fastq_rev = None):
  # takes input reads, tries to assign them bins based on barcodes, then
  # returns trimmed reads and bin name

  fwd_match = match_bc(fastq_fwd.sequence, fwd_bcs)

  if fastq_rev == None:
    # merged read
    rev_match = match_bc(fastq_fwd.sequence[::-1], rev_bcs)
  else:
    rev_match = match_bc(fastq_rev.sequence, rev_bcs)

  if fwd_match == False or rev_match == False:
    # couldn't match to both fwd and rev barcodes
    return False
  else:
    if fastq_rev == None:
      # merged read
      fastq_fwd.sequence = fastq_fwd.sequence[fwd_match[0]+len(fwd_match[1]):-(rev_match[0]+len(rev_match[1]))]
    else:
      fastq_fwd.sequence = fastq_fwd.sequence[fwd_match[0]+len(fwd_match[1]):]
      fastq_rev.sequence = fastq_rev.sequence[rev_match[0]+len(rev_match[1]):]

    return (fastq_fwd, fastq_rev, fwd_match[1], rev_match[1])

def fast_fastq(fp_in):
    buf = []

    for line in fp_in:
        buf.append(line)

        if len(buf) == 4:
            yield fastq_record(buf, 0)
            buf = []

class fastq_record():
    def __init__(self, lines, offset):
        lines = [x.strip() for x in lines]

        self.id = lines[0][1:]
        self.sequence = lines[1]
        self.quals = lines[3]
        self.offset = offset

    def raw(self):
        return "\n".join(["@%s" % (self.id,), self.sequence, "+", self.quals, ""])

def qual_filter(read):
  if mean([ord(x) - options.phred_offset for x in read.quals]) < options.min_qual:
    return False
  else:
    return True

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
        fwd_bcs[l["name"]] = l["barcode"])

      if l["orientation"] == "reverse":
        rev_bcs[l["name"]] = l["barcode"])

  return (fwd_bcs, rev_bcs)

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <miseq_fastq.zip> <sample_name> <barcodes.txt>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./miseq_out]",
                    default="./miseq_out",
                    help="write output files to this directory")

  parser.add_option("-p",
                    dest="num_threads",
                    type="int",
                    metavar="[1]",
                    default=1,
                    help="number of threads used during analysis")

  parser.add_option("-m",
                    dest="mem_size",
                    metavar="[1G]",
                    default="1G",
                    help="amount of memory to use during read merging")

  parser.add_option("-q",
                    dest="min_qual",
                    type="int",
                    metavar="[25]",
                    default=25,
                    help="discard reads with mean quality below this value")

  parser.add_option("--merge",
                    dest="merge",
                    action="store_true",
                    default=False,
                    help="merge forward and reverse reads with PEAR")

  parser.add_option("--phred_offset",
                    dest="phred_offset",
                    type="int",
                    metavar="[33]",
                    default=33,
                    help="phred quality offset (default 33 for illumina 1.8+)")

  parser.add_option("--min-qual-perc",
                    dest="min_qual_perc",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="warn if fewer than this % of reads pass quality filter")

  parser.add_option("--min-merged-perc",
                    dest="min_merged_perc",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="warn if fewer than this % of reads can be merged")

  options, args = parser.parse_args(arguments)

  if len(args) <> 3:
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

  if not os.path.isfile(args[0]):
    print "Error: Specified sequencing zip file does not exist"
    parser.print_help()
    sys.exit(1)

  if not os.path.isfile(args[2]):
    print "Error: Specified barcode file does not exist"
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
  fwd_bcs, rev_bcs = load_barcodes(args[2])

  # open zipfile
  miseq_zip = zipfile.ZipFile(args[0])

  # find pairs
  pairs = dict([("_".join(x.split("/")[-1].split(".")[0].split("_")[:-2]), (x, y)) for x, y in find_pairs(miseq_zip)])

  if args[1] not in pairs:
    raise ValueError("Could not find %s in %s!" % (args[1], args[0]))

  barcode_to_count = {}

  if options.merge:
    sys.stderr.write("Merging reads...\n")

    with open(os.path.join(options.output_dir, "merged_reads.assigned.fastq"), "w") as assigned, \
         open(os.path.join(options.output_dir, "merged_reads.unassigned.fastq"), "w") as unassigned:
      min_quality_read_length = 1e10
      max_quality_read_length = 0
      total_quality_read_length = 0.0

      total_reads = 0.0
      quality_reads = 0.0

      # extract read files
      fwd = fetch_unzip(miseq_zip, pairs[args[1]][0])
      rev = fetch_unzip(miseq_zip, pairs[args[1]][1])
  
      # run pear to merge reads
      merged, stats = pear(fwd.name, rev.name)
  
      if stats["assembled_reads"] < options.min_merged_perc:
        sys.stderr.write("  Warning: only %.02f%% of reads assembled\n" % stats["assembled_reads"])
  
      # read through fastq
      for merged_read in fast_fastq(open("%s.assembled.fastq" % merged.name, "r")):
        total_reads += 1
  
        # check that read passes quality filter
        if not qual_filter(merged_read):
          continue

        # keep track of stats
        quality_reads += 1

        if len(merged_read.sequence) < min_read_length:
          min_quality_read_length = len(merged_read.sequence)

        if len(merged_read.sequence) > max_read_length:
          max_quality_read_length = len(merged_read.sequence)
  
        total_quality_read_length += len(merged_read.sequence)

        # demultiplex
        dm_out = demultiplex(merged_read, fwd_bcs.values(), rev_bcs.values(), None)
 
        if dm_out == False:
          # strip pair info from read and write to unassigned file
          merged_read.id = merged_read.id.split(" ")[0]
          unassigned.write(merged_read.raw())
        else:
          # rename read to barcodes used and write to assigned file
          f_read, r_read, f_bc, r_bc = dm_out

          bc_name = "%s_%s" % (f_bc, r_bc)

          try:
            barcode_to_count[bc_name] += 1
          except KeyError:
            barcode_to_count[bc_name] = 1

          f_read.id = "%s_%s" % (bc_name, barcode_to_count[bc_name])
          assigned.write(f_read.raw())
  
      # remove temporary files
      os.unlink(fwd.name)
      os.unlink(rev.name)
      os.unlink(merged.name)
      os.unlink("%s.assembled.fastq" % merged.name)
      os.unlink("%s.discarded.fastq" % merged.name)
      os.unlink("%s.unassembled.forward.fastq" % merged.name)
      os.unlink("%s.unassembled.reverse.fastq" % merged.name)
  
      if total_reads > 0:
        if quality_reads / total_reads < options.min_qual_perc:
          sys.stderr.write("  Warning: only %.02f%% of reads passed quality filter (mean(qv) > %s)\n" % (quality_reads * 100 / total_reads, options.min_qual))

      sys.stderr.write("\nSummary")
      sys.stderr.write("\n  Total reads:       %d" % total_reads)
      sys.stderr.write("\n  Quality reads:     %d" % quality_reads)
      sys.stderr.write("\n  Min read length:   %d" % min_quality_read_length)
      sys.stderr.write("\n  Mean read length:  %d" % (total_quality_read_length / quality_reads))
      sys.stderr.write("\n  Max read length:   %d\n" % max_quality_read_length)
      sys.stderr.write("\n  Assigned reads:    %d" % sum(barcode_to_count.values()))
      sys.stderr.write("\n  Unassigned reads:  %d" % quality_reads - sum(barcode_to_count.values()))
      sys.stderr.write("\n  Avg reads/barcode: %d\n" % mean(barcode_to_count.values()))
  else:
    sys.stderr.write("Filtering reads...\n")

    with open(os.path.join(options.output_dir, "fwd_reads.assigned.fastq"), "w") as f_assigned, \
         open(os.path.join(options.output_dir, "rev_reads.assigned.fastq"), "w") as r_assigned, \
         open(os.path.join(options.output_dir, "fwd_reads.unassigned.fastq"), "w") as f_unassigned, \
         open(os.path.join(options.output_dir, "rev_reads.unassigned.fastq"), "w") as r_unassigned:
      total_reads = 0.0
      quality_reads = 0.0

      # extract read files
      fwd = fetch_unzip(miseq_zip, pairs[args[1]][0])
      rev = fetch_unzip(miseq_zip, pairs[args[1]][1])
  
      # read through foward and reverse fastq simultaneously
      for f_read, r_read in izip(fast_fastq(open(fwd.name, "r")), fast_fastq(open(rev.name, "r"))):
        total_reads += 1
  
        # check that read passes quality filter
        if qual_filter(f_read) == False or qual_filter(r_read) == False:
          continue

        quality_reads += 1

        # demultiplex
        dm_out = demultiplex(f_read, fwd_bcs.values(), rev_bcs.values(), r_read)
 
        if dm_out == False:
          # strip pair info from read and write to unassigned file
          f_read.id = f_read.id.split(" ")[0]
          f_unassigned.write(f_read.raw())

          r_read.id = r_read.id.split(" ")[0]
          r_unassigned.write(r_read.raw())
        else:
          # rename reads to barcodes used and write to assigned file
          f_read, r_read, f_bc, r_bc = dm_out

          bc_name = "%s_%s" % (f_bc, r_bc)

          try:
            barcode_to_count[bc_name] += 1
          except KeyError:
            barcode_to_count[bc_name] = 1

          f_read.id = "%s_%s" % (bc_name, barcode_to_count[bc_name])
          f_assigned.write(f_read.raw())
  
          r_read.id = "%s_%s" % (bc_name, barcode_to_count[bc_name])
          r_assigned.write(r_read.raw())

      if total_reads > 0:
        if quality_reads / total_reads < options.min_qual_perc:
          sys.stderr.write("  Warning: only %.02f%% of reads passed quality filter (mean(qv) > %s)\n" % (quality_reads * 100 / total_reads, options.min_qual))

      sys.stderr.write("\nSummary")
      sys.stderr.write("\n  Total reads:       %d" % total_reads)
      sys.stderr.write("\n  Quality reads:     %d" % quality_reads)
      sys.stderr.write("\n  Assigned reads:    %d" % sum(barcode_to_count.values()))
      sys.stderr.write("\n  Unassigned reads:  %d" % (quality_reads - sum(barcode_to_count.values())))
      sys.stderr.write("\n  Avg reads/barcode: %d\n" % mean(barcode_to_count.values()))

if __name__ == "__main__":
  main()
