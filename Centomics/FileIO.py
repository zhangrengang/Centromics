from xopen import xopen as open

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
