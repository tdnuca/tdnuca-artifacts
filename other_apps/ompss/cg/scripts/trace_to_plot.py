#!/usr/bin/env python3

from __future__	import print_function
from math		import ldexp
from os.path	import splitext, exists, basename
from os			import remove, SEEK_SET
from sys		import argv, exit, stderr

import re

conv_data = re.compile(":9300001:(\d+):9300002:(\d+):9300003:(\d+)")
test_data = re.compile("930000[0-3]")

for (i,name) in enumerate(argv[1:], 1):
	base, ext = splitext(name)
	symfilename=base+'.pcf'
	plotname=base + '.plot'
	print("reading {}: checking {} for symbols and writing plot data to {}".format(name, symfilename, plotname), file=stderr)

	assert ext == '.prv', "File has wrong extension to be a paraver trace"

	try:
		with open(symfilename, 'r') as symbolfile:
			symbols=sum(1 for l in symbolfile if test_data.search(l))

			if symbols == 0:
				print("WARNING none of the symbols used in this script are defined in the symbol file {} ".format(symfilename))
			elif symbols != 4:
				print("WARNING some symbols used in this script are not defined in the symbol file {}".format(symfilename))
	except IOError:
		print("WARNING could not open the symbol file {} ".format(symfilename))
		continue

	try:
		with open(plotname, 'w') as plotfile, open(name, 'r') as tracefile:
			print("iteration	norm_gradient_squared	time (ns)", file=plotfile);

			errors=0
			start_solving=0
			n=0
			pos=0

			# see when first data appears to know first ns-timestamp
			while True :
				(pos,line) = (tracefile.tell(),tracefile.readline())
				if conv_data.search(line):
					start_solving=int(line.split(":", 6)[5])
					break

			# rewind to before that line so we can plot it as well
			tracefile.seek(pos, SEEK_SET)

			# now iterate over data in file and add to our output files
			for line in tracefile:
				findings = conv_data.search(line)
				if findings:
					numbers=list(map(int, findings.groups()))
					time=int(line.split(":", 6)[5]) - start_solving
					print("{}	{}	{}".format(numbers[0], ldexp(numbers[2],numbers[1]-1023), time), file=plotfile)
	except IOError:
		print("WARNING could not open either input or output file {}.{{prv,plot}}".format(base))
		continue

	title=basename(base).replace('_', ' ')

	time="set title 'Convergence plot'; set ylabel 'log(error)'; set xlabel 'time (s)'; plot '{}' using (\\$3/1e6):(log(\\$2)) with lines lt {} title '{}'".format(plotname, i, title)
	it  ="set title 'Convergence plot'; set ylabel 'log(error)'; set xlabel 'iteration'; plot '{}' using 1:(log(\\$2)) with lines lt {} title '{}'".format(plotname, i, title)

	print("Success ! To plot convergence over time (resp. per iteration), use the following :\ngnuplot -p -e \"{}\"\ngnuplot -p -e \"{}\"".format(time, it))

