#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/26/2014"
__version__ = 1.0

import resource
from optparse import OptionParser, OptionGroup
from numpy import mean
from fastq import fast_fastq, fastq_to_fasta
from itertools import izip, izip_longest
import tempfile
import subprocess
import time
import sys
import os

PATH_TO_USEARCH = "usearch"

def usearch(input_fname, output_fname, params = []):
  usearch = subprocess.Popen([PATH_TO_USEARCH,
                              "-usearch_global", input_fname,
                              "-db", args[0],
                              "-blast6out", output_fname,
                              "-strand", "both",
                              "-id", str(options.identity),
                              "-threads", str(options.num_threads)] + \
                             params,
                             stderr=subprocess.STDOUT,
                             stdout=open(os.devnull))

  while usearch.poll() == None:
    time.sleep(0.1)

  if usearch.returncode <> 0:
    raise ValueError("USEARCH exited with an error (%s)" % usearch.returncode)

def extract_unaligned_lomem(b6_fname, in_fname, out_fname):
  # extracts unaligned reads using a small memory footprint--probably the
  # smallest possible that doesnt require reads to be sorted

  b6_cache = set([])
  in_cache = {}

  with open(b6_fname, "r") as b6_fp, \
       open(in_fname, "r") as in_fp, \
       open(out_fname, "w") as out_fp:
    for b6_l, in_r in izip_longest(b6_fp, fast_fastq(in_fp), fillvalue=False):
      if b6_l:
        b6_l = b6_l.split("\t", 1)

        try:
          del in_cache[b6_l[0]]
        except KeyError:
          b6_cache.add(b6_l[0])

      if in_r:
        try:
          b6_cache.remove(in_r.id)
        except KeyError:
          in_cache[in_r.id] = in_r.raw()

    # dump all the reads that remain
    out_fp.write("".join(in_cache.values()))

def parse_usearch(fwd_b6, rev_b6, out_fname, merged_b6 = False):
  count = 0

  if merged_b6:
    with open(merged_b6, "r") as in_fp, \
         open(out_fname, "w") as out_fp:
      for l in in_fp:
        l = l.strip().split("\t")

        # can be multiple 16S sequences per taxon, e.g. 10F2_1 and 10F2_2
        l[1] = l[1].rsplit("_", 1)[0]

        out_fp.write("%s\t%s\n" % tuple(l[:2]))
        count += 1

      return count
  else:
    with open(fwd_b6, "r") as fwd_fp, \
         open(rev_b6, "r") as rev_fp, \
         open(out_fname, "w") as out_fp:
      fwd_to_taxon = {}
      rev_to_taxon = {}

      for f_l, r_l in izip(fwd_fp, rev_fp):
        f_l = f_l.strip().split("\t")
        r_l = r_l.strip().split("\t")

        # can be multiple 16S sequences per taxon, e.g. 10F2_1 and 10F2_2
        f_l[1] = f_l[1].rsplit("_", 1)[0]
        r_l[1] = r_l[1].rsplit("_", 1)[0]

        try:
          # reverse read matches forward read
          if rev_to_taxon[f_l[0]] == f_l[1]:
            #print "match! %s, fwd = %s, rev = %s" % (f_l[0], f_l[1], rev_to_taxon[f_l[0]])
            out_fp.write("%s\t%s\n" % (f_l[0], f_l[1]))
            del rev_to_taxon[f_l[0]]
            count += 1
        except KeyError:
          # we haven't come across the matching reverse read yet
          #print "saving fwd %s -> %s" % (f_l[0], f_l[1])
          fwd_to_taxon[f_l[0]] = f_l[1]

        try:
          # forward read matches reverse read
          if fwd_to_taxon[r_l[0]] == r_l[1]:
            #print "match! %s, rev = %s, fwd = %s" % (r_l[0], r_l[1], fwd_to_taxon[r_l[0]])
            out_fp.write("%s\t%s\n" % (r_l[0], r_l[1]))
            del fwd_to_taxon[r_l[0]]
            count += 1
        except KeyError:
          # we haven't come across the matching foward read yet
          #print "saving rev %s -> %s" % (r_l[0], r_l[1])
          rev_to_taxon[r_l[0]] = r_l[1]

    return (fwd_to_taxon, rev_to_taxon, count)

def parse_options(arguments):
  global options, args

  parser = OptionParser(usage="%prog [options] <-f fwd_reads.fastq -r rev_reads.fastq|-m merged_reads.fastq> <db.fa>",
                        version="%prog " + str(__version__))

  parser.add_option("-o",
                    dest="output_dir",
                    metavar="[./intermediate]",
                    default="./intermediate",
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

  parser.add_option("-p",
                    dest="num_threads",
                    type="int",
                    metavar="[1]",
                    default=1,
                    help="number of threads used during analysis")

  parser.add_option("--id",
                    dest="identity",
                    type="int",
                    metavar="[95]",
                    default=95,
                    help="minimum percent identity for matches")

  parser.add_option("--save-unaligned",
                    dest="save_unaligned",
                    action="store_true",
                    default=False,
                    help="seperate out reads that didn't align")

  options, args = parser.parse_args(arguments)

  if len(args) <> 1:
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
    print "Error: Specified database file (%s) does not exist" % args[0]
    parser.print_help()
    sys.exit(1)

  options.identity /= 100.0

  if not os.path.exists(options.output_dir):
    os.mkdir(options.output_dir)
  elif not os.path.isdir(options.output_dir):
    print "Error: Specified path exists and is not a directory"
    parser.print_help()
    sys.exit(1)

def main():
  parse_options(sys.argv[1:])

  if options.merged_fname:
    sys.stderr.write("Converting merged reads to FASTA...")
    merged_fasta = fastq_to_fasta(options.merged_fname)
    merged_b6 = os.path.join(options.output_dir, "merged_reads.usearch.b6")
    sys.stderr.write("\nRunning USEARCH on merged reads...")
    usearch(merged_fasta.name, merged_b6)

    sys.stderr.write("\n\nParsing USEARCH results...")
    aligned_reads = parse_usearch(False, False, os.path.join(options.output_dir, "assigned_taxa.txt"), merged_b6)

    sys.stderr.write("\n\nSummary")
    sys.stderr.write("\n  Total reads:       %d" % merged_fasta.read_count)
    sys.stderr.write("\n  Aligned reads:     %d" % aligned_reads)
    sys.stderr.write("\n  Unaligned reads:   %d\n" % (merged_fasta.read_count - aligned_reads))

    if options.save_unaligned:
      sys.stderr.write("\nFinding and writing unaligned reads...")
      unaligned_fasta = os.path.join(options.output_dir, "merged_reads.unaligned.fa")
      extract_unaligned_lomem(merged_b6, options.merged_fname, unaligned_fasta)
      sys.stderr.write("\n")
  else:
    sys.stderr.write("Converting forward reads to FASTA...")
    fwd_fasta = fastq_to_fasta(options.fwd_fname)
    fwd_b6 = os.path.join(options.output_dir, "fwd_reads.usearch.b6")
    sys.stderr.write("\nRunning USEARCH on forward reads...")
    usearch(fwd_fasta.name, fwd_b6)

    sys.stderr.write("\n\nConverting reverse reads to FASTA...")
    rev_fasta = fastq_to_fasta(options.rev_fname)
    rev_b6 = os.path.join(options.output_dir, "rev_reads.usearch.b6")
    sys.stderr.write("\nRunning USEARCH on reverse reads...")
    usearch(rev_fasta.name, rev_b6)

    sys.stderr.write("\n\nParsing USEARCH results...")
    fwd_discord, rev_discord, aligned_pairs = parse_usearch(fwd_b6, rev_b6, os.path.join(options.output_dir, "assigned_taxa.txt"))

    total_reads = fwd_fasta.read_count + rev_fasta.read_count
    aligned_reads = aligned_pairs * 2 + len(fwd_discord) + len(rev_discord)
    discordant_pairs = set(fwd_discord.keys()).intersection(rev_discord.keys())
    orphaned_reads = set(fwd_discord.keys()).symmetric_difference(rev_discord.keys())

    sys.stderr.write("\n\nSummary")
    sys.stderr.write("\n  Total reads:       %d" % total_reads)
    sys.stderr.write("\n  Aligned reads:     %d" % aligned_reads)
    sys.stderr.write("\n  Unaligned reads:   %d\n" % (total_reads - aligned_reads))
    sys.stderr.write("\n  Concordant pairs:  %d" % aligned_pairs)
    sys.stderr.write("\n  Discordant pairs:  %d" % len(discordant_pairs))
    sys.stderr.write("\n  Orphaned reads:    %d\n" % len(orphaned_reads))

    if options.save_unaligned:
      sys.stderr.write("\nFinding and writing unaligned forward reads...")
      fwd_unaligned_fasta = os.path.join(options.output_dir, "fwd_reads.unaligned.fa")
      extract_unaligned_lomem(fwd_b6, options.fwd_fname, fwd_unaligned_fasta)

      sys.stderr.write("\nFinding and writing unaligned reverse reads...")
      rev_unaligned_fasta = os.path.join(options.output_dir, "rev_reads.unaligned.fa")
      extract_unaligned_lomem(rev_b6, options.rev_fname, rev_unaligned_fasta)
      sys.stderr.write("\n")

if __name__ == "__main__":
  main()
