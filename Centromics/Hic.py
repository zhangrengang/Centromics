from .FileIO import Mnd
from .Bin import bin_data
from .Juicer import hic2signals

def count_links(mnd, fout, bin_size=10000, ncpu=8):
	'''count links between different chromosomes'''
	d_count = Mnd(mnd, ncpu=ncpu, method='imap_unordered', chunksize=2000).count_links(
				  bin_size=bin_size, diff_chr=True)
	data = [(chr, pos, 'diff_chrom', count) for (chr, pos), count in d_count.items()]
	data = sorted(data, key=lambda x: (x[0], x[1]))
	bin_data(data, fout, bin_size=bin_size)

def count_obs(inHic, chrLst, prefix='pre', tmpdir='tmp', norm='NONE', bin_size=10000, ncpu=8, overwrite=1):
	out1, out2 = prefix + '.inter_chr', prefix + '.intra_chr'
	cmd_opts = {'tc_tasks': ncpu, 'mode': 'local', 'retry': 2, 'cont': overwrite, 'fail_exit':False}
	with open(out1, 'w') as fout1, open(out2, 'w') as fout2:
		hic2signals(fout1, fout2, overwrite=overwrite,
			inHic=inHic, inChrLst=chrLst, prefix=tmpdir, res=bin_size, norm=norm, cmd_opts=cmd_opts)
	return out1, out2
