import os
from collections import OrderedDict
from .RunCmdsMP import  run_cmd
def bin_data(data, fout, bin_size=10000):
	'''
three format:
1. chrom, pos, 
2. chrom, pos, key
3. chrom, pos, key, value
'''
	d_bin = OrderedDict()
	d_keys = OrderedDict()
	for line in data:
		chrom, pos, *key = line
		if len(key) == 0:
			key = 'None'
			count = 1
		elif len(key) == 2:
			key, count = key
			count = float(count)
		else:
			key = key[0]
			count = 1
		bin = bin_line(chrom, pos, bin_size=bin_size)
		if bin not in d_bin:
			d_bin[bin] = {}
		# count
		try: d_bin[bin][key] += count
		except KeyError: d_bin[bin][key] = count
		# record
		try: d_keys[key] += count
		except KeyError: d_keys[key] = count

	keys = list(d_keys.keys())
	line = ['#chrom', 'start', 'end'] + keys
	fout.write('\t'.join(line)+'\n')

	for bin, d_count in d_bin.items():
		chrom, start = bin
		start = start*bin_size
		end = start + bin_size
		values = [d_count.get(key, 0) for key in keys]
		line = [chrom, start, end] + values
		line = map(str, line)
		fout.write('\t'.join(line)+'\n')

def bin_line(chrom, pos, key=None, bin_size=1000):
	bin = chrom, pos // bin_size
	return bin

def bin_bam(bamfile, outfile, bin_size=10000, opts='--ignoreDuplicates', ncpu=8):
	import pysam
	bam = pysam.Samfile(bamfile, 'rb')
	try: bam.check_index()
	except ValueError:
		pysam.index(bamfile, '-@', str(ncpu))
	bam.close()

	opts += ' --numberOfProcessors {}'.format(ncpu)	
	cmd = 'bamCoverage --bam {} --outFileName {} --outFileFormat bedgraph --binSize {} {}'.format(
		bamfile, outfile, bin_size, opts)
	run_cmd(cmd, log=True, fail_exit=True)
	
