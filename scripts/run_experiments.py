#!/usr/bin/python

import sys, subprocess

def run_doe(doe_in, output_file):
	with open(doe_in) as fin, open(output_file, 'w') as fout:
		next(fin)
		fout.write("size threads sequential_time parallel_time builtin_time\n")
		for row in fin:
			size, tlevel = map(float, row.strip().split())
			batcmd = "./src/parallelQuicksort %d %d" % (size, tlevel)
			result = subprocess.check_output(batcmd, shell=True).strip()
			fout.write("%d %d %s\n" % (size, 2**(tlevel-1), result))

if len(sys.argv) != 3:
	print "usage: python run.py <experiments.csv> <output.csv>"
	exit(1)

print "Start running the experiments.."
run_doe(sys.argv[1], sys.argv[2])
print "Execution is done."
