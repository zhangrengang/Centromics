
import os
import sys
import itertools
import collections
from .small_tools import mkdirs
from .RunCmdsMP import run_job, logger

bindir = os.path.dirname(os.path.realpath(__file__))
juicebox_bin = f'{bindir}/bin/juicebox_tools.jar'

def run_juicerbox(inHic, chr_combs, outdir='matrix', res=10000, norm="NONE", cmd_opts={}):
	mkdirs(outdir)
	d_files = {}
	cmds = []
	for chr1, chr2 in chr_combs:
		outfile = '%s/%s-%s.%s.mat' % (outdir, chr1, chr2, res)
		cmd = f'java -jar {juicebox_bin} dump observed {norm} {inHic} {chr1} {chr2} BP {res} {outfile}' 
		#print(cmd)
		d_files[(chr1, chr2)] = outfile
		cmds += [cmd]
	cmd_file = outdir.rstrip('/') + f'{res}.sh'
	run_job(cmd_file, cmd_list=cmds, **cmd_opts)
	return d_files
	
def juicer2centrome(chr_combs, outbed, outmatix, ChrLst, d_id, outdir='matrix', res=10000):
	f1 = open(outbed, 'w')
	f2 = open(outmatix, 'w')
	i = 0
	d_id0 = {}
	points = []
	lines = []
	for chr1, chr2 in chr_combs:
		outfile = '%s/%s-%s.%s.mat' % (outdir, chr1, chr2, res)
		for line in open(outfile):
			pos1, pos2, cov = line.rstrip().split()
			pos1 = int(pos1)
			pos2 = int(pos2)
			lines.append([(chr1, pos1), (chr2, pos2), cov])
			if not (chr1, pos1) in d_id0:
				points.append((chr1, pos1))
			if not (chr2, pos2) in d_id0:
				points.append((chr2, pos2))
			d_id0[(chr1, pos1)] = i
			d_id0[(chr2, pos2)] = i
			i = i+2
	
	for (chr, pos), i in list(d_id.items()):
		line = [chr, pos, pos+res, i]
		line = list(map(str, line))
		print('\t'.join(line), file=f1)
	for (chr1, pos1), (chr2, pos2), cov in sorted(lines, key=lambda x: (ChrLst.index(x[0][0]), x[0][1], ChrLst.index(x[1][0]), x[1][1])):
		id1 = d_id[(chr1, pos1)]
		id2 = d_id[(chr2, pos2)]
		line3 = [id1, id2, cov]
		line3 = list(map(str, line3))
		print('\t'.join(line3), file=f2)
	f1.close()
	f2.close()
def parse_chr(inChrLst, res=10000):
	d = collections.OrderedDict()
	i = 0
	for line in open(inChrLst):
		temp = line.strip().split()
		chr, length = temp[:2]
		length = int(length)
		for pos in range(0, length+1, res):
			key = (chr, pos)
			d[key] = i
			i += 1
	return d
def run_centrion(outbed, outmatix, outfig, res=10000):
	cmd = 'plot_finding_centromeres.py %s %s %s %s 2> %s.err' % (outbed, outmatix, outfig, res, outfig)
	print(cmd)
	os.system(cmd)
def hic2matrix(inHic=None, inChrLst=None, prefix=None, res=100000, figfmt='pdf', self=True, cent=False, norm="NONE", cmd_opts={}):
	if self:
		combinations = itertools.combinations_with_replacement
	else:
		combinations = itertools.combinations
	outbed = '%s.%s.bed' % (prefix, res)
	outmatix = '%s.%s.matrix' % (prefix, res)
	outfig = '%s.%s.%s' % (prefix, res, figfmt)
	if isinstance(inChrLst, str):
		ChrLst = [line.strip().split()[0] for line in open(inChrLst)]
	else:
		ChrLst = inChrLst
	chr_combs = list(combinations(ChrLst, 2))
	outdir = '%s.matrix' % (prefix, )
	d_files = run_juicerbox(inHic, chr_combs, outdir=outdir, res=res, norm=norm, cmd_opts=cmd_opts)
	if cent:
		d_id = parse_chr(inChrLst, res=res)
		juicer2centrome(chr_combs, outbed, outmatix, ChrLst, d_id, outdir=outdir, res=res)
		run_centrion(outbed, outmatix, outfig, res=res) # res not > 40k # reviced: res=res0, not applied
	return d_files
def hic2signals(fout1, fout2, res=1000, **kargs):
	d_files = hic2matrix(res=res, **kargs)
	logger.info('Count inter and intra-chromosomal signals')
	distance = 20*res
	for d_count, fout in [(count_diff_chrom(d_files), fout1), (count_same_chrom(d_files, distance), fout2)]:
		for (chr, bin), count in sorted(d_count.items()):
			signal = sum(count) / len(count) if isinstance(count, list) else count # normalize
			line = [chr, bin, bin+res, signal]
			line = map(str, line)
			print('\t'.join(line), file=fout)

def count_diff_chrom(d_files):
	d_count = {}
	for (chr1, chr2), matfile in d_files.items():
		if chr1 == chr2: 
			continue
		for (bin1, bin2, observed) in JuicerMatrix(matfile):
			for chr, bin in ([chr1, bin1], [chr2, bin2]):
				key = (chr, bin)
				try: d_count[key] += observed
				except KeyError: d_count[key] = observed
	return d_count
def count_same_chrom(d_files, distance):
	d_count = {}
	for (chr1, chr2), matfile in d_files.items():
		if chr1 != chr2:
			continue
		lines = list(JuicerMatrix(matfile))
		bins = sorted(set([b1 for b1, *_ in lines]))
		for bin in bins:
			key = (chr1, bin)
			left, right, cross = 0, 0, 0
			x, y, z = 0,0,0
			for bin1, bin2, observed in lines:
				assert bin2 >= bin1, [chr1, bin1, bin2]
				if bin1+distance < bin < bin2-distance:
					#try: d_count[key] += [observed]
					#except KeyError: d_count[key] = [observed]
					cross += observed
					z += 1
				elif bin2-distance < bin and bin2-bin1 > distance:	# left
					left += observed
					x += 1
				elif bin1+distance > bin and bin2-bin1 > distance:
					right += observed
					y += 1
				else:
					continue
			if z == 0:
				continue
			signal = cross /( left + right) 
			d_count[key] = signal
	return d_count

class JuicerMatrix:
	def __init__(self, matfile):
		self.matfile = matfile
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.matfile):
		#	yield JuicerMatrixLine(line)
			temp = line.strip().split()
			yield int(temp[0]), int(temp[1]), float(temp[2])

class JuicerMatrixLine:
	def __init__(self, line):
		temp = line.strip().split()
		self.bin1, self.bin2, self.observed = int(temp[0]), int(temp[1]), float(temp[2])

if __name__ == '__main__':
	inHic=sys.argv[1]
	inChrLst=sys.argv[2]
	try: args = {'res': int(sys.argv[4]), 'self':True}
	except: args = {}
	try: args['cent'] = sys.argv[5]
	except: pass
	hic2matrix(**args)
