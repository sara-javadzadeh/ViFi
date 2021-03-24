#!/usr/bin/python
# -*- coding: utf-8 -*-

import pysam, os
import argparse
from time import clock
from collections import Counter
from distutils.version import LooseVersion


def is_float(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def read_scores_file(hmm_file):
  input = open(hmm_file, 'r')
  scores = dict()
  for line in input:
    results = line.strip().split(',')
    scores[results[0]] = [float(x) if is_float(x) else x for x in results]
  input.close()
  return scores

def read_map(map_file, scores, delimiter = '\t', direction = 'reverse'):
  input = open(map_file, 'r')
  map = dict()
  for line in input:
    results = line.strip().split(delimiter)
    if results[1] in scores:
      results[0] = results[0].replace('/1','').replace('/2','').replace('@','')
      if results[0] in map and map[results[0]][2] < scores[results[1]][2]:
        map[results[0]] = scores[results[1]]
      elif results[0] not in map:
        map[results[0]] = scores[results[1]]
  input.close()
  return map


parser = \
  argparse.ArgumentParser(description='Merges viral reads identified via HMMs'
              )

parser.add_argument(
  '--trans',
  dest='transName',
  help='Output BAM file of trans reads between hg19 and viral',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )


parser.add_argument(
  '--unknown',
  dest='unknownName',
  help='BAM file of reads with only one read mapped to between hg19 and other to unknown',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--threshold',
  dest='threshold',
  help='E-value threshold to accept as virus read',
  metavar='float',
  action='store',
  type=float,
  default=0.01,
  )


parser.add_argument(
  '--reduced',
  dest='reducedName',
  help='CSV file containing HMM scores',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--map',
  dest='mapName',
  help='Map file containing mapping of scores and sequences',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--output',
  dest='outputName',
  help='Output BAM file of updated trans reads',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )
args = parser.parse_args()

transFile = pysam.Samfile(args.transName[0], 'rb')

# Add a reference to the samfile for each hmm used in scores.
hmm_references = []
num_non_hmm_refs = None
scores = read_scores_file(args.reducedName[0])
hmm_indices_with_mapped_reads = list(set([int(score[3]) for score in scores.values()]))
for hmm_index in hmm_indices_with_mapped_reads:
  ref_name = 'viral_hmm_{}'.format(str(hmm_index))
  ref_len = 8000
  hmm_references.append({'LN': ref_len, 'SN': ref_name})

#Hack to deal with pysam 0.14 or greater not being able to edit headers
if LooseVersion(pysam.__version__) <= LooseVersion("0.13.0"):
  references = transFile.header
  num_non_hmm_refs = len(references['SQ'])
  for ref in hmm_references:
      references['SQ'].append(ref)
else:
  references = transFile.header.to_dict()
  num_non_hmm_refs = len(references['SQ'])
  for ref in hmm_references:
      references['SQ'].append(ref)

  outputFile = pysam.Samfile(args.outputName[0], 'wb', header=references)
  outputFile.close()
  os.system('samtools reheader %s %s > %s.fixed' % (args.outputName[0], args.unknownName[0], args.unknownName[0]))
  os.system('mv %s.fixed %s' % (args.unknownName[0], args.unknownName[0]))


unknownFile = pysam.Samfile(args.unknownName[0], 'rb')
outputFile = pysam.Samfile(args.outputName[0], 'wb', header=references)
mapping = read_map(args.mapName[0], scores)
for read in unknownFile.fetch(until_eof=True):
  # Note: threshold is set and checked in run_hmms.py through nhmmer internal tools.
  if read.qname in mapping:
    read_or_mate_tid = num_non_hmm_refs + hmm_indices_with_mapped_reads.index(mapping[read.qname][3])
    if read.is_unmapped:
      read.tid = read_or_mate_tid
      read.is_unmapped = False
      read.pos = mapping[read.qname][4]
      read.cigartuples = [(0,read.qlen)]
    else:
      read.mate_is_unmapped = False
      read.mpos = mapping[read.qname][4]
      read.mrnm = read_or_mate_tid
    outputFile.write(read)

for read in transFile.fetch(until_eof=True):
  outputFile.write(read)

outputFile.close()

