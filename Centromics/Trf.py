import sys,os,glob
from collections import OrderedDict
from Bio import SeqIO
from .Repeat import TRF
from .split_records import split_fastx_by_chunk_num, bin_split_fastx_by_chunk_num
from .small_tools import mkdirs
from .RunCmdsMP import run_job, logger, run_cmd
from .Blast import blast,BlastOut
from .Bin import bin_data

#OPT='2 7 7 80 10 50 1000 -d'
OPT='1 1 2 80 5 200 2000 -d -h'

def run_trf(inseq, nbins=20, tmpdir='tmp', tr_opts=OPT, window_size=1e6, overwrite=False):
	mkdirs(tmpdir)
	prefix = '{}/chunks'.format(tmpdir, )
	_,_,_, chunk_files = bin_split_fastx_by_chunk_num(
			inseq, prefix=prefix, chunk_num=nbins, window_size=1e6, window_ovl=0,
			seqfmt='fasta', suffix='')
	cmds = []
	for chunk_file in chunk_files:
		cmd = 'cd {tmpdir} && trf {genome} {opts} > /dev/null; ls {genome}.*.dat'.format(
				tmpdir=tmpdir, genome=os.path.realpath(chunk_file), opts=tr_opts)
		cmds +=[cmd]
	cmd_file = '{}/cmds.list'.format(tmpdir)
	run_job(cmd_file, cmd_list=cmds, tc_tasks=nbins/2, fail_exit=False, mode='local', cont=(not overwrite))
	datfiles = glob.glob('{}.*.dat'.format(prefix))
	return datfiles
class Trf:
	def __init__(self, inseq, tr_opts=OPT, nbins=20, tmpdir='tmp', window_size=1e6, overwrite=False):
		self.inseq = inseq
		self.tr_opts = tr_opts
		self.nbins = nbins
		self.tmpdir = tmpdir
		self.overwrite = overwrite
		self.window_size = window_size
	def __iter__(self):
		return self._run_trf()
	def _run_trf(self):
		datfiles = run_trf(self.inseq, nbins=self.nbins, tmpdir=self.tmpdir, tr_opts=self.tr_opts,
					window_size=self.window_size, overwrite=self.overwrite)
		for datfile in datfiles:
			for rc in TRF(datfile):
				yield rc
	def reads_trf(self, fout, min_cov=0.8):
		d_cons = {}
		d_len = {}
		i = 0
		for rc in self.get_best_trf(min_cov=min_cov):
			if rc.cons_seq in d_cons:
				id = d_cons[rc.cons_seq]
				d_len[id] += len(rc)
				continue
			i += 1
			id = 'TR{}'.format(i)
			d_cons[rc.cons_seq] = id
			d_len[id] = len(rc)
			desc = 'source={};cov={:.2f};copy={};size={}'.format(rc.chrom, rc.cov, rc.copy_number, rc.consensus_size)
			print('>{} {}\n{}'.format(id, desc, rc.cons_seq), file=fout)
		return d_len
	def get_best_trf(self, min_cov=0.8):
		d_seqlen = self.seqlen
		rcs = []
		last_chrom = ''
		for rc in self:
			cov = len(rc) / d_seqlen[rc.chrom]
			if cov < min_cov:
				continue
			rc.cov = cov
			if last_chrom and rc.chrom != last_chrom:
				yield self.best_trf(rcs)
				rcs = []
			rcs += [rc]
			last_chrom = rc.chrom
		yield self.best_trf(rcs)
	def best_trf(self, rcs):
		return max(rcs, key=lambda x:(x.cov, -x.consensus_size))
	@property
	def seqlen(self):
		return {rc.id: len(rc.seq) for rc in SeqIO.parse(self.inseq, 'fasta')}

def filter_trf_family(trfseq, trfmcl, d_trf_len, fout, total_len=None, min_ratio=0.05):
	from REPcluster.Mcl import MclGroup
	f = open(trfseq+'.info', 'w')
	rcs = []
	lens, ratios = [], []
	for rc, grp in zip(SeqIO.parse(trfseq, 'fasta'), MclGroup(trfmcl)):
		grp_len = sum([d_trf_len[trf] for trf in grp])
		ratio = grp_len/total_len
		rc.ratio = ratio
		lens += [len(rc)]
		ratios += [ratio*100]
		rcs += [rc]
		line = [rc.id, len(rc), len(grp), ratio]
		print('\t'.join(map(str, line)), file=f)
	f.close()
	plot_trf(lens, ratios, outfig=trfseq+'.png')

	# get top
	for i, rc in enumerate(sorted(rcs, key=lambda x:-x.ratio)):
		if i == 0:
			max_ratio= rc.ratio
		if rc.ratio / max_ratio < min_ratio:
			continue
		rc.description += ';ratio={:.3%}'.format(rc.ratio)
		SeqIO.write(rc, fout, 'fasta')

def trf_map(trfseq, genome, fout, min_cov=0.9, ncpu=4, window_size=10000, tmpdir='tmp/',
		blast_opts='-task blastn-short -word_size 9 -dust no -soft_masking false'):
	blast_outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand'"
	blast_out = trfseq + '.blastout'
	# db=trfseq, query=genome
	blast_out = blast(trfseq, genome, seqtype='nucl', blast_out=blast_out, blast_outfmt=blast_outfmt, 
				blast_opts=blast_opts, ncpu=ncpu)

	logger.info('Parse blast out')	
	data = []
	for rc in BlastOut(blast_out, blast_outfmt):
#		print(dir(rc))
		if rc.scov < min_cov:
			continue
		line = (rc.qseqid, rc.qstart, rc.sseqid)
		data += [line]

	# sort
	seqids = {rc.id:i for i, rc in enumerate(SeqIO.parse(trfseq, 'fasta'))}
	data = sorted(data, key=lambda x: (x[0], seqids[x[2]], x[1]))

	# bin
	bin_data(data, fout, bin_size=window_size)
	
	
def plot_trf(lens, ratios, outfig, xlab='Repeat length', ylab='Genomic ratio (%)'):
	import matplotlib.pyplot as plt
	plt.figure()	
	for x,y in zip(lens, ratios):
		plt.vlines(x, 0, y)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.savefig(outfig)

def trf_plot(inseq, inbed='', outbed=sys.stdout, nbins=64, window_size=20000, min_percent=5, min_copy=10, min_unit_size=5):
	seqs = []
	for rc in SeqIO.parse(inseq, 'fasta'):
		seqs += [(rc.id, len(rc.seq))]
	if os.path.exists(inbed):
		d_trf = load_from_bed(inbed)
		outfig = inbed + '.pdf'
	else:
		d_trf = load_from_seq(inseq)
		outfig = inseq + '.trf.pdf'
	#print >> sys.stderr, seqs, d_trf.keys()
	d_coord = {}
	f = open(inseq+'.coord', 'w')
	for chrom, _ in seqs:
		if chrom not in d_trf:
			continue
		d_bins = OrderedDict()
		lines = d_trf[chrom]
		lines = resolve_overlaps(lines)
		last_end = 0
		for rc in lines:
			line = [chrom, start, end, period_size, copy_number, cons_seq]  = rc
			if end <= last_end:
				continue
			if start < last_end:
				start = last_end
			last_end = end

			if not os.path.exists(inbed):
				line = map(str, line)
				print >>outbed, '\t'.join(line)
			bin1 = start/window_size
			bin2 = end/window_size
			last = start-1
			for bin in range(bin1, bin2+1):
				_end = (bin+1)*window_size
				span = min(_end, end) - last
#				span = range(last, min(end, _end))
				try: d_bins[bin] += span
				except KeyError: d_bins[bin] = span
				last = _end
		for bin, span in d_bins.items():
#			span = len(set(span))
			percent = 1e2*span/window_size
			if percent < min_percent:
				continue
			start = bin *window_size
			try: d_coord[chrom] += [(start, percent)]
			except KeyError: d_coord[chrom] = [(start, percent)]
			line = [chrom, start, start+window_size, percent]
			line = map(str, line)
			print >>f, '\t'.join(line)
	f.close()

	#print >> sys.stderr, d_coord
	# seqs: [(id, len), ...]
	# d_coord: {id: (pos, value), ...}
	plot_telomere(seqs, d_coord, outfig)

def main():
	inseq=sys.argv[1]
	try: inbed = sys.argv[2]
	except IndexError: inbed = ''
	trf_plot(inseq, inbed)

def resolve_overlaps(lines, max_ovl=10, ):
	'''assume only one chromsome'''
	last_line = None
	discards = []
	ie, io = 0, 0
	for line in sorted(lines, key=lambda x:x[1]):
		discard = None
	#	print(last_line, line)
		if last_line:
			if line == last_line:	# equal
				ie += 1
				line_pair = [last_line, line]	# retain, discard
			else:
				if overlap(line, last_line) > max_ovl:
					io += 1
					if (line[2]-line[1]) > (last_line[2]-last_line[1]):	# length is prior
						line_pair = [line, last_line]
					else: 
						line_pair = [last_line, line]
				else:	# no overlap or too short overlap
					last_line = line
					continue
			
			retain, discard = line_pair

			discards += [discard]

		if not last_line or discard != line:
			last_line = line
	logger.info('discard {} equal and {} overlapped hits; {} in total'.format(ie, io, ie+io))
	return sorted(set(lines) - set(discards), key=lambda x:x[1])
def overlap(self, other):
	ovl = max(0, min(self[2], other[2]) - max(self[1], other[1]))
	return 100*ovl/(min((self[2]-self[1]+1), (other[2]-other[1]+1)))

def load_from_seq(inseq, nbins=32, min_percent=5, min_copy=10, min_unit_size=5):
	datfiles = run_trf(inseq, nbins)

	d_trf = {}
	for datfile in datfiles:
		for rc in TRF(datfile):
			if rc.copy_number < min_copy:
				continue
			if rc.period_size < min_unit_size:
				continue
			if rc.chrom not in d_trf:
				d_trf[rc.chrom] = []
			line = (rc.chrom, rc.start, rc.end, rc.period_size, rc.copy_number, rc.cons_seq)
			d_trf[rc.chrom] += [line]
	return d_trf
def load_from_bed(inded):
	d_trf = {}
	for line in open(inded):
		line = line.strip().split('\t')
		line[1:3] = map(int, line[1:3])
		chrom = line[0]
		if chrom not in d_trf:
			d_trf[chrom] = []
		d_trf[chrom] += [line]
	return d_trf
if __name__ == '__main__':
	main()
