#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrewscz@gmail.com"
__date__ = "3/19/2014"
__version__ = 1.0

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
