from xopen import xopen as open
from .RunCmdsMP import pool_func

class Fastx:
	def __init__(self, seqfile, seqfmt=None):
		self.seqfile = seqfile
		self.seqfmt = seqfmt
		if seqfmt is not None:
			assert seqfmt in {'fasta', 'fastq'}
	def __iter__(self):
		return self._parse()

	def _parse(self):
		if isinstance(self.seqfile, str):
			handle = open(self.seqfile)
		else:	# file
			handle = self.seqfile

		i = 0
		lines = []
		for line in handle:
			i += 1
			if i == 1 and self.seqfmt is None:
				if line[0] == '>':
					self.seqfmt = 'fasta'
				elif line[0] == '@':
					self.seqfmt = 'fastq'
				else:
					raise ValueError('Unknown sequence format, neither fasta nor fastq')
			if lines and (\
					(self.seqfmt == 'fasta' and line[0] == '>') or (self.seqfmt == 'fastq' and i % 4 == 1)):
				yield FastxRecord(lines, self.seqfmt)
				lines = []
			lines.append(line)
		yield FastxRecord(lines, self.seqfmt)

class FastxRecord:
	def __init__(self, lines, seqfmt='fasta'):
		self.description = lines[0][1:].rstrip()
		self.seqfmt = seqfmt
		if seqfmt == 'fasta':
			self.seq = ''.join(lines[1:]).rstrip()
		elif seqfmt == 'fastq':
			self.seq = lines[1].rstrip()
			self.qual = lines[3].rstrip()
	@property
	def id(self):
		return self.description.split()[0]
	def __len__(self):
		return len(self.seq)
	def write(self, fout, seqfmt=None):
		if seqfmt is None:
			seqfmt = self.seqfmt
		if seqfmt == 'fasta':
			print('>{}\n{}'.format(self.description, self.seq), file=fout)
		elif seqfmt == 'fastq':
			print('>{}\n{}\n+\n{}'.format(self.description, self.seq, self.qual), file=fout)

class Mnd:
	def __init__(self, mnd, ncpu=4, method='imap', chunksize=1000):
		self.mnd = mnd
		self.ncpu = ncpu
		self.method = method
		self.chunksize  = chunksize
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.mnd):
			yield MndLongLine(line)
#		iterable = open(self.mnd)
#		jobs = pool_func(_parse_mnd_line, iterable, 
#				processors=self.ncpu, method=self.method, chunksize=self.chunksize)
#		for job in jobs:
#			yield job
	def count_links(self, diff_chr=True, bin_size=1000):
		d_count = {}
		for rc in self:
			if diff_chr and rc.chr1 == rc.chr2:
				continue
			for chr, pos in [(rc.chr1, rc.pos1), (rc.chr2, rc.pos2)]:
				bin = chr, pos // bin_size * bin_size
				try: d_count[bin] += 1
				except KeyError: d_count[bin] = 1
		return d_count
def _parse_mnd_line(line):
	return MndLongLine(line)
class MndLongLine:
	keys = ['str1', 'chr1', 'pos1', 'frag1', 'str2', 'chr2', 'pos2', 'frag2',
		'mapq1', 'cigar1', 'sequence1', 'mapq2', 'cigar2', 'sequence2', 'readname1', 'readname2']
	types = [str, str, int, int, str, str, int, int, int, str, str, int, str, str, str, str]
	def __init__(self, line):
		self._line = line
		self.line = line.strip().split()
		self.set_attr()
	def set_attr(self):
		for key, type, value in zip(self.keys, self.types, self.line):
			setattr(self, key, type(value))
	def write(self, fout):
		fout.write(self._line)

