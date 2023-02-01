import numpy as np

class CustomBdgParser:
	def __init__(self, inBdg):
		self.bdg = inBdg
	def parse(self):	# return regions by chrom
		lastchr = ''
		regions = []
		keys = [None]
		d_regions = {}
		for line in open(self.bdg):
			temp = line.strip().split('\t')
			if line.startswith('#'):
				keys = temp[3:]
				continue
			if not d_regions:
				d_regions = {key: {} for key in keys}
			chrom, start, end = temp[0], int(temp[1]), int(temp[2])
			if not chrom == lastchr:
				i = 0
			i += 1
			values = map(float, temp[3:])
			for key, value in zip(keys, values):
				region = Region(chrom, start, end, value, index=i)
				try: d_regions[key][chrom] += [region]
				except KeyError: d_regions[key][chrom] = [region]
			lastchr = chrom
		return d_regions

class Regions():
	def __init__(self, regions):
		self.regions = regions
		self.data = [region.value for region in regions]
	def __iter__(self):
		return iter(self.regions)
	def percentile(self, tile, Abs=True):
		return np.percentile(self.data, tile)
	def max(self):
		return np.max(self.data)
	def min(self):
		return np.min(self.data)
class Region():
	def __init__(self, *args, **kargs):
		self.chr = args[0]
		self.start = args[1]
		self.end = args[2]
		self.value = args[3]
		self.index = kargs['index']
	def __hash__(self):
		return hash((self.chr, self.start, self.end))
	def __eq__(self, other):
		return (self.chr, self.start, self.end) == (other.chr, other.start, other.end)
def find_peaks(d_bdgfiles, outbed, **kargs):
	#print(d_bdgfiles)
	lines = []
	for label, bdgfile in d_bdgfiles.items():
		if bdgfile is None:	
			continue
		d_regions = CustomBdgParser(bdgfile).parse()
		for key, d_chrs in d_regions.items():
			_from = '{}-{}'.format(label, key) if key else label
#			print(_from)
			for chrom, regions in d_chrs.items():
				peaks = findPeak(Regions(regions), **kargs)
				for chrom, start, end, pvalue,svalue in peaks:
					line = [chrom, start, end, _from, pvalue, svalue]
					lines += [line]
	lines = sorted(lines, )
	line = ['#chrom', 'start', 'end', 'data_from', 'peak_value', 'sum_value']
	outbed.write('\t'.join(line)+'\n')
	for line in lines:
		line = map(str, line)
		outbed.write('\t'.join(line)+'\n')

def findPeak(regions, upper=0.6, lower=0.2, lower_tile=10, max_dist=1, min_nsites=3, bin_size=100e3, **kargs):
	'''same chrom'''
	_max, _min = regions.max(), regions.min()
	__min = _min
#	_min = min(max(_min, _max*lower), regions.percentile(lower_tile))
	_min = max(_min, min(regions.percentile(lower_tile), _max*0.5))
	diff = _max - _min
#	upthd = _min + diff*upper
	upthd = _max*upper
	lowthd = _min + diff*lower
#	print(_max, _min, __min, regions.percentile(lower_tile), upthd, lowthd)
	up_regions = set([]) # [site1, site2]
	for rc in regions:
		if rc.value >= upthd:
			up_regions.add(rc)
	sig_sites = [] # [[site1, site2],[site5, site6]]
	for rc in regions:
		if rc.value < lowthd:
			continue
		try:
			last = sig_sites[-1][-1]
			if 0 <= rc.start - last.end <= max_dist*bin_size:
				sig_sites[-1] += [rc]
			else:
				sig_sites += [[rc]]
		except IndexError: # first value
			sig_sites = [[rc]]
	peaks = []
	for sites in sig_sites:
		if len(sites) < min_nsites:
			continue
		if not set(sites) & up_regions:
			continue
		chr = sites[0].chr
		start = sites[0].start
		end = sites[-1].end
		values = [site.value for site in sites]
		line = [chr, start, end, np.max(values), np.sum(values)]
		peaks += [line]
	return peaks
