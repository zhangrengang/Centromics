from .FileIO import Mnd
from .Bin import bin_data

def count_links(mnd, fout, bin_size=10000, ncpu=8):
	'''count links between different chromosomes'''
	d_count = Mnd(mnd, ncpu=ncpu, method='imap_unordered', chunksize=2000).count_links(
				  bin_size=bin_size, diff_chr=True)
	data = [(chr, pos, 'diff_chrom', count) for (chr, pos), count in d_count.items()]
	data = sorted(data, key=lambda x: (x[0], x[1]))
	bin_data(data, fout, bin_size=bin_size)
