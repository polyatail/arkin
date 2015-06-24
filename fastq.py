# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "2/26/2014"
__version__ = 1.0

import tempfile
import copy

def fast_fastq(fp):
    buf = []

    for line in fp:
        buf.append(line)

        if len(buf) == 4:
            yield fastq_record(buf)
            buf = []

class fastq_record(object):
    def __init__(self, lines):
        lines = [x.strip() for x in lines]

        self.id = lines[0][1:]
        self.sequence = lines[1]
        self.quals = lines[3]

    def __getitem__(self, val):
      new = copy.copy(self)
      new.sequence = self.sequence[val]
      new.quals = self.quals[val]

      return new

    def raw(self):
        return "\n".join(["@%s" % (self.id,), self.sequence, "+", self.quals, ""])

def fastq_to_fasta(input_fname):
  outfile = tempfile.NamedTemporaryFile(delete=True)

  count = 0

  for read in fast_fastq(open(input_fname, "r")):
    outfile.write(">%s\n%s\n" % (read.id, read.sequence))
    count += 1

  outfile.flush()
  outfile.read_count = count

  return outfile

def count_reads_fasta(input_fname):
  outfile = tempfile.NamedTemporaryFile(delete=True)

  count = 0

  for l in open(input_fname, "r"):
    if l.startswith(">"):
      count += 1

  outfile.read_count = count

  return outfile
