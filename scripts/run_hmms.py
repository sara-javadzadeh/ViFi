import sys, os, time,pysam
import argparse
import tempfile,shutil

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--threads', type=int, action='store',
                      help='''Threads to use''', required=False, default = 1)
  parser.add_argument('--threshold', type=float, action='store',
                      help='''Threshold for E-values of the reads to be included in nhmmer. See http://eddylab.org/software/hmmer/Userguide.pdf''', required=False, default = 0.02)
  parser.add_argument('-b', '--bamfile', type=str, action='store',
                      help='''Bamfile to analyze''', required=True)
  parser.add_argument('-d', '--directory', default=".",
                      help='''output directory for all files created during run (default: current directory)''')
  parser.add_argument('-H', '--hmm_list', default=None, type=str, required=True,
                      help='''File containing HMMs to use''')
  options = parser.parse_args()
  return options

def read_hmm_file(hmm_list):
  hmms = {}
  input = open(hmm_list, "r");
  idx = 0;
  for line in input:
    #Todo, instead of index, note which viral family the hmm came from
    hmms[idx] = line.strip()
    idx+=1
  input.close()
  return hmms

def run_pipeline(options):
  hmms = read_hmm_file(options.hmm_list)
  if not os.path.exists('%s/logs' % options.directory):
    os.mkdir('%s/logs' % options.directory)
  if not os.path.exists('%s/temp' % options.directory):
    os.mkdir('%s/temp' % options.directory)
  #if (not os.path.exists('%s/temp/unmapped.fas' % options.directory)) or (os.path.getsize('%s/temp/unmapped.fas' % options.directory) == 0):
  prepare_unmapped_sequences(options)
  #Now search HMMs against reads
  print "Running HMMs"
  start_time = time.time()
  for i in hmms.keys():
    start_time_for_each_hmm = time.time()
    print "\tRunning HMM %s" % hmms[i]
    # --noali is to indicate that alignments are not specified in output to reduce output volume
    # --incE Use an E-value of <= <x> as the inclusion threshold.
    # The default is 0.01, meaning that on average, about 1 false positive
    # would be expected in every 100 searches with different query sequences
    command = 'nhmmer -o %s/temp/hmmsearch.%s --noali --incE %.2f --cpu %d %s %s/temp/unmapped.fas' % (options.directory, str(i), options.threshold, options.threads, hmms[i], options.directory)
    print(command)
    os.system(command)
    duration_time = time.time() - start_time_for_each_hmm;
    print "Finished running against HMMs %s: %fs" % (hmms[i], duration_time)
  end_time = time.time() - start_time;
  print "Finished running against HMMs: %fs" % end_time
  print "Processing results\n";
  scores = {}
  # This is to store hmm ids that improved scores for some read
  for i in hmms.keys():
    scores = read_nhmmer_result("%s/temp/hmmsearch.%s" % (options.directory, i), i, scores)
  output = open('%s/temp/reduced.csv' % options.directory, 'w')
  for score in scores.keys():
    output.write("%s\n" % (','.join([str(s) for s in scores[score]])))
  output.close()
def read_nhmmer_result(file, hmm_index, scores):
  input = open(file, 'r')
  start_line = 'Query:'
  start = False
  for line in input:
    if start == False and line.find(start_line) == 0:
      start = True
      foo = input.next()
      foo = input.next()
      foo = input.next()
    if start == True and line.find('read_') != -1:
      res = line.split()
      # line header: E-value    score   bias    Sequence    start   end     description
      # Appending to scores: (sequence_id, E-value, float(score), hmm_index, start, end)
      candidate_entry = (res[3], res[0], float(res[1]), hmm_index, res[4], res[5])
      if scores.setdefault(res[3], list(candidate_entry))[2] < float(res[1]):
        scores[res[3]] = list(candidate_entry)
    elif line.strip() == "" and start == True:
      break
    # Reads that pass the inclusion threshold will come before the line:
    #  ------ inclusion threshold ------
    if line.find("------ inclusion threshold ------") != -1:
      break
  return scores

def prepare_unmapped_sequences(options):
  start_time = time.time()
  counter = 0;
  bam = pysam.Samfile(options.bamfile, 'rb')
  fas = open('%s/temp/unmapped.fas' % options.directory, 'wb')
  map = open('%s/temp/unmapped.map' % options.directory, 'wb')

  for read in bam:
      fas.write('>read_%d\n%s\n' % (counter, read.seq))
      map.write('%s\tread_%d\n' % (read.qname, counter))
      counter+=1
  fas.close()
  map.close()
  bam.close()
  end_time = time.time() - start_time;
  print "Prepared sequences for searching against HMMs: %fs" % end_time

if __name__ == '__main__':
  start_time = time.time()
  options = parse_args()
  run_pipeline(options)
