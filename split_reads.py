#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "6/22/2015"
__version__ = 1.0

import re
from optparse import OptionParser, OptionGroup
from fastq import fast_fastq, fastq_to_fasta
import sys
import os

def load_destinations(fname, formats):
  dests = {}

  # sample destination_file
  for l in open(fname, "r"):
    if l.startswith("#"):
      continue

    l = l.strip().split("\t")

    if l[0] in dests:
      raise ValueError("Sample %s appears twice in destination file" % l[0])

    if l[1] not in formats:
      raise ValueError("Sample %s has destination %s which is not in formats file" % (l[0], l[1]))

    dests[l[0]] = l[1]

  return dests

def load_formats(fname):
  formats = {}

  # format_name fastq_or_fasta min_read_length max_read_length regexp
  for l in open(fname, "r"):
    if l.startswith("#") or not l.strip():
      continue

    l = l.strip().split("\t")

    if l[0] in formats:
      raise ValueError("Format %s appears more than one in formats file" % l[0])

    if l[1] not in ("fastq", "fasta"):
      raise ValueError("Format %s: %s not one of fastq or fasta" % (l[0], l[1]))

    try:
      l[2] = int(l[2])
      l[3] = int(l[3])
    except ValueError:
      raise ValueError("Format %s: %s and %s must both be integers" % (l[0], l[2], l[3]))

    if l[4] == "*":
      regexp = False
    else: 
      regexp = re.compile(l[4])

    formats[l[0]] = (l[1], l[2], l[3], regexp)

  return formats

def open_files(formats):
  fps = {}

  for f in ["notfound"] + formats.keys():
    if options.merged_fname:
      out_path = os.path.join(options.output_dir, "%s.merged.%s" % (f, "fastq" if f[1] == "fastq" else "fa"))

      if os.path.isfile(out_path):
        raise ValueError("Destination %s already exists!" % out_path)

      fps[f] = (open(out_path, "w"), )
    else:
      fwd_path = os.path.join(options.output_dir, "%s.fwd.%s" % (f, "fastq" if f[1] == "fastq" else "fa"))
      rev_path = os.path.join(options.output_dir, "%s.rev.%s" % (f, "fastq" if f[1] == "fastq" else "fa"))

      if os.path.isfile(fwd_path):
        raise ValueError("Destination %s already exists!" % fwd_path)

      if os.path.isfile(rev_path):
        raise ValueError("Destination %s already exists!" % rev_path)

      fps[f] = (open(fwd_path, "w"), open(rev_path, "w"))

  return fps

def split_file(input_fp, fps, dests, formats, read_index, stats):
  for read in fast_fastq(input_fp):
    sample = read.id.split("_")[0]

    try:
      read_dest = dests[sample]
    except KeyError:
      #sys.stderr.write("read: %s notfound\n" % read.id)
      stats["notfound"]["matched"] += 1
      fps["notfound"][read_index].write(read.raw())
      continue

    # fastq_or_fasta min_read_length max_read_length regexp
    #sys.stderr.write("read: %s -> %s\n" % (read.id, read_dest))
    read_format = formats[read_dest]

    # read too short
    if read_format[1] != -1 and len(read.sequence) < read_format[1]:
      stats[read_dest]["short"] += 1
      #sys.stderr.write("  len(read) = %s < %s\n" % (len(read.sequence), read_format[1]))
      continue

    # read too long
    if read_format[2] != -1 and len(read.sequence) > read_format[2]:
      stats[read_dest]["long"] += 1
      #sys.stderr.write("  len(read) = %s > %s\n" % (len(read.sequence), read_format[2]))
      continue

    # read didn't match regular expression
    if read_format[3]:
      if not read_format[3].match(read.sequence):
        stats[read_dest]["regexp"] += 1
        #sys.stderr.write("  regexp match failed\n")
        continue

    # convert to fasta or don't
    if read_format[0] == "fasta":
      #sys.stderr.write("  convert to fasta\n")
      raw_read = ">%s\n%s\n" % (read.id, read.sequence)
    else:
      raw_read = read.raw()

    # write it out
    stats[read_dest]["matched"] += 1
    #sys.stderr.write("  write to %s\n\n" % fps[read_dest][read_index])
    fps[read_dest][read_index].write(raw_read) 

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <-f fwd_reads.fastq -r rev_reads.fastq|-m merged_reads.fastq> <destinations.txt> <formats.txt>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./split_reads]",
                    default="./split_reads",
                    help="write output files to this directory")

  parser.add_option("-f",
                    dest="fwd_fname",
                    type="str",
                    default=False,
                    help="path to FASTQ containing forward reads")

  parser.add_option("-r",
                    dest="rev_fname",
                    type="str",
                    default=False,
                    help="path to FASTQ containing reverse reads")

  parser.add_option("-m",
                    dest="merged_fname",
                    type="str",
                    default=False,
                    help="path to FASTQ containing merged reads")

  options, args = parser.parse_args(arguments)

  if len(args) <> 2:
    print "Error: Incorrect number of arguments"
    parser.print_help()
    sys.exit(1)

  if options.merged_fname and (options.fwd_fname or options.rev_fname):
    print "Error: When specifying -m, cannot specify -f or -r"
    parser.print_help()
    sys.exit(1)
  elif not (options.merged_fname or options.fwd_fname or options.rev_fname):
    print "Error: Must specify either -m or -f and -r"
    parser.print_help()
    sys.exit(1)
  elif bool(options.fwd_fname) ^ bool(options.rev_fname):
    print "Error: When using -f or -r, both must be specified"
    parser.print_help()
    sys.exit(1)

  for fname in (options.merged_fname, options.fwd_fname, options.rev_fname):
    if fname != False and \
       not os.path.isfile(fname):
      print "Error: Specified FASTQ file (%s) does not exist" % fname
      parser.print_help()
      sys.exit(1)

  if not os.path.isfile(args[0]):
    print "Error: Specified destinations file (%s) does not exist" % args[0]
    parser.print_help()
    sys.exit(1)

  if not os.path.isfile(args[1]):
    print "Error: Specified formats file (%s) does not exist" % args[1]
    parser.print_help()
    sys.exit(1)

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specified path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

def main():
  parse_options(sys.argv[1:])

  formats = load_formats(args[1])
  dests = load_destinations(args[0], formats)
  fps = open_files(formats)

  # make dict to keep track of stats
  stats = {}

  for f in ["notfound"] + formats.keys():
    stats[f] = {"matched": 0,
                "short": 0,
                "long":  0,
                "regexp": 0}

  if options.merged_fname:
    sys.stderr.write("Filtering merged reads...")
    split_file(open(options.merged_fname, "r"), fps, dests, formats, 0, stats)

  else:
    sys.stderr.write("Filtering forward reads...")
    split_file(open(options.fwd_fname, "r"), fps, dests, formats, 0, stats)
    
    sys.stderr.write("\n\nFiltering reverse reads...")
    split_file(open(options.rev_fname, "r"), fps, dests, formats, 1, stats)

  # write stats out to file
  with open(os.path.join(options.output_dir, "split_stats.log"), "w") as fp:
    fp.write("\t".join(["#format", "matched", "too_short", "too_long", "regexp_no_match"]) + "\n")

    for f in stats:
      fp.write("\t".join(map(str, [f, stats[f]["matched"], stats[f]["short"], stats[f]["long"], stats[f]["regexp"]])) + "\n")

  sys.stderr.write("\n")

if __name__ == "__main__":
  main()
