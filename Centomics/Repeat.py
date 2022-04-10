import sys, os
import re
import glob
import argparse
try: import networkx as nx
except: pass
from itertools import combinations
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
try: from Gff import GffLine, GffExons
except: pass
from .RunCmdsMP import run_cmd, logger

def makeArgparse():
	parser = argparse.ArgumentParser( \
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("genome", action="store",type=str,
					help="input genome sequence in fastaformat [required]")
	addRepeatArgs(parser)
	args = parser.parse_args()
	return args

def addRepeatArgs(parser):
	parser.add_argument('-dr_opts', type=str, default='-p -d -l 50 -e 3',
					help="options for dispersed repeats by vmatch [default='%(default)s']")
	parser.add_argument('-tr_opts', type=str, default='2 7 7 80 10 50 800 -d',
					help="options for tandem repeats by TRF [default='%(default)s']")
	parser.add_argument('-sr_opts', type=str, default='1-10 2-5 3-4 4-3 5-3 6-3',
					help="options for simple repeats by regular expression [default='%(default)s']")


class RepeatPipeline():
	def __init__(self, genome, tmpdir='/dev/shm/tmp', 
				prefix=None,
				dr_opts='-p -d -l 50 -e 3',
				tr_opts='2 7 7 80 10 50 800 -d',
				sr_opts='1-10 2-5 3-4 4-3 5-3 6-3',
				**kargs
				):
		self.genome = os.path.realpath(genome)
		if prefix is None:
			prefix = os.path.basename(self.genome)
		#self.prefix = prefix
		self.tmpdir = tmpdir
		self.prefix = '{}/{}'.format(self.tmpdir, prefix)
		self.dr_opts = dr_opts
		if tr_opts.find('-d') < 0 :
			tr_opts += ' -d'
		self.tr_opts = tr_opts
		self.sr_opts = sr_opts
		
	def run(self):
		# dispersed repeats: DR
		vmatch_out = self.run_vmatch()
		d_no = self.get_seq_no()
		graph1 = Vmatch(vmatch_out, d_no=d_no).cluster()

		# tandem repeats: TR
		trf_datfile = self.run_trf()
		graph2 = TRF(trf_datfile).cluster()
		
		# simple repeats: SR
		graph3 = SSR(self.genome).cluster()

		nodes = graph1.nodes() + graph2.nodes() + graph3.nodes()
		nodes = sorted(nodes, key=lambda x:(x.chrom, x.start, -x.end))
		dedup_nodes = [nodes[0]]
		for node in nodes[1:]:
			if dedup_nodes[-1].contains(node):
				continue
			dedup_nodes += [node]
		dedup_nodes = self.remove_singleton(dedup_nodes)
		records = []	
		for node in dedup_nodes:
			record = node.to_exons()
			records += [record]
		for record in sorted(records, key=lambda x:x.start):
			record.write(sys.stdout)
		return records	
	def remove_singleton(self, nodes, source='vmatch'):
		d_names = {}
		for node in nodes:
			if not node.source == source:
				continue
			try: d_names[node.name] += [node]
			except: d_names[node.name] = [node]
		remove_nodes = []
		for xnodes in d_names.values():
			if len(xnodes) == 1:
				remove_nodes += xnodes
		remove_nodes = set(remove_nodes)
		left_nodes = []
		for node in nodes:
			if node in remove_nodes:
				continue
			left_nodes += [node]
		return left_nodes
	def run_vmatch(self):
		cmd = 'mkvtree -db {} -indexname {} -dna -allout -pl'.format(
				self.genome, self.prefix)
		run_cmd(cmd, log=True)
		vmatch_out = '{}.vmatch'.format(self.prefix)
		cmd = 'vmatch {} {} > {}'.format(self.dr_opts, self.prefix, vmatch_out)
		run_cmd(cmd, log=True)
		return vmatch_out
	def run_trf(self):
		cmd = 'cd {tmpdir} && trf {genome} {opts}'.format(
				tmpdir=self.tmpdir, genome=self.genome, opts=self.tr_opts)
		run_cmd(cmd, log=True)
		datfiles = glob.glob('{}/{}.*.dat'.format(self.tmpdir, os.path.basename(self.genome)))
		assert len(datfiles) == 1
		return datfiles[0]
	def get_seq_no(self):
		return {i:rc.id for i, rc in enumerate(SeqIO.parse(self.genome, 'fasta'))}

class RepeatGraph(nx.Graph):
	def __init__(self):
		super(RepeatGraph, self).__init__()
	def number_repeat(self, prefix='REP'):
		d_names = {}
		for i, cmpt in enumerate(nx.connected_components(self)):
			name = '{}{}'.format(prefix, i+1)
			for j, node in enumerate(sorted(cmpt, key=lambda x:x.start)):
				rid = '{}-{}'.format(name, j+ 1)
				node.id = rid
				node.name = name
				d_names[node] = (rid, name)
		for node in self.nodes():
			node.id, node.name = d_names[node]

	def set_node_attributes(self, **kargs):
		for node in self.nodes():
			for key, value in kargs.items():
				setattr(node, key, value)

class RepeatSegemnt:
	def __init__(self, chrom, start, end):	# 1-based
		self.chrom, self.start, self.end = chrom, start, end
		self.strand = '+'
		self.type = 'repeat_region'
	def __hash__(self):
		return hash(self.key)
	def __eq__(self, other):
		if self.key == other.key:
			return True
		return False
	def __len__(self):
		return self.end - self.start +1
	def __str__(self):
		return '{chrom}:{start}-{end}'.format(**self.__dict__)
	@property
	def key(self):
		return (self.chrom, self.start, self.end)
	def to_gff(self):
		source = self.source
		attributes = OrderedDict()
		attributes['ID'] = self.id #getattr(self, 'id', None)
		attributes['Name'] = self.name # getattr(self, 'name', None)
		attributes['rpt_type'] = self.rpt_type
		attributes['rpt'] = getattr(self, 'rpt', None)
		line = [self.chrom, source, self.type, self.start, self.end, 
					'.', self.strand, '.', attributes]
		return GffLine(line)
	def to_exons(self):
		record = GffExons([self.to_gff()])
		record.chrom = self.chrom
		record.source = self.source
		record.start = self.start
		record.end = self.end
		record.strand = self.strand
		record.gene = record.name = self.name
		record.id = self.id
		record.product = None
		record.gene_id = self.name
		record.rna_id = self.id
		record.rna_type = 'repeat'
		record.trans_splicing = None
		return record

	def has_overlap(self, other):
		if not self.chrom == other.chrom:
			return False
		return max(0, min(self.end, other.end) - max(self.start, other.start))
	def contains(self, other):
		if not self.chrom == other.chrom:
			return False
		if self.start <= other.start <= self.end and self.start <= other.end <= self.end:
			return True
		return False

class Vmatch:
	def __init__(self, matchfile, d_no=None, prefix='DR', rpt_type='dispersed', source='vmatch'):
		self.matchfile = matchfile
		self.d_no = d_no
		self.prefix = prefix
		self.rpt_type = rpt_type
		self.source = source
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.matchfile):
			if line.startswith('#'):
				continue
			line = line.strip().split()
			if not len(line) == 11:
				continue
			record = VmatchRecord(line, self.d_no)
			yield record
	def remove_duplicates(self):
		matches = sorted(self._parse(), key=lambda x: (x.left_id, x.left_start, -x.left_len))
		
	def cluster(self, G=None):
		if G is None:
			G = RepeatGraph()
		for match in self:
			G.add_edge(match.left_segment, match.right_segment, match=match)
		G.number_repeat(prefix=self.prefix)

		G.set_node_attributes(rpt_type=self.rpt_type, source=self.source)
		return G
	

class VmatchRecord(object):
	def __init__(self, line, d_no=None):
		self.title = ['left_len', 'left_no', 'left_pos', 'type', 
					'right_len', 'right_no', 'right_pos',
					'distance', 'evalue', 'score', 'identity']
		self.types = [int, int, int, str,
					int, int, int,
					int, float, float, float]
		self.line = line
		self._set_attr(self.title, self.line, self.types)
#		print >>sys.stderr, vars(self)
		if d_no is not None:
			self.left_id, self.right_id = d_no[self.left_no], d_no[self.right_no]
		self.left_start, self.left_end = self.left_pos+1, self.left_pos+self.left_len
		self.right_start, self.right_end = self.right_pos+1, self.right_pos+self.right_len
		self.left_segment = RepeatSegemnt(self.left_id, self.left_start, self.left_end)
		self.right_segment = RepeatSegemnt(self.right_id, self.right_start, self.right_end)
		#self.left_segment.type = self.type
	def _set_attr(self, keys, values, types):
		for key, value, type in zip(keys, values, types):
			setattr(self, key, type(value))

class TRF(object):
	def __init__(self, datfile, prefix='TR', rpt_type='tandem', source='trf'):
		self.datfile = datfile
		self.prefix = prefix
		self.rpt_type = rpt_type
		self.source = source
	def __iter__(self):
		return self._parse()
#	def __len__(self):
#		return len(list(self))
	def _parse(self):
		for line in open(self.datfile):
			line = line.strip().split()
			if not line:
				continue
			if line[0] == 'Sequence:':
				chrom = line[1]
			elif len(line) == 15:
				yield TRFRecord(chrom, line)
	def cluster(self, G=None):
		if G is None:
			G = RepeatGraph()
		for match1, match2 in combinations(self, 2):
			G.add_node(match1.segment)
			G.add_node(match2.segment)
			if match1.unit_key == match2.unit_key:
				G.add_edge(match1.segment, match2.segment, match=(match1,match2))
		G.number_repeat(prefix=self.prefix)
		G.set_node_attributes(rpt_type=self.rpt_type, source=self.source)
		return G
				
class TRFRecord(VmatchRecord):
	def __init__(self, chrom, line, base=1):
		self.chrom = chrom
		self.line = line
		self.title = ['start', 'end', 'period_size', 'copy_number', 'consensus_size',
					'percent_matches', 'percent_indels', 'score', 'A', 'C', 'G', 'T', 'entropy',
					'cons_seq', 'seq']
		self.types = [int, int, int, float, int,
					int, int, int, int, int, int, int, float,
					str, str]
		self._set_attr(self.title, self.line, self.types)
		self.segment = RepeatSegemnt(self.chrom, self.start, self.end)
		self.unit = self.cons_seq
		self.segment.unit = self.unit
		self.base = base
		self.add_bin()
	@property
	def unit_key(self):
		key = self.unit, str(Seq(self.unit).reverse_complement())
		return tuple(sorted(key))	
	def __len__(self):
		return self.end-self.start+1
	def add_bin(self):
		try:
			bin_chrom, bin_start, bin_end = re.compile(r'(\S+):(\d+)-(\d+)').match(self.chrom).groups()
			bin_start = int(bin_start)
			self.chrom = bin_chrom
			self.start += bin_start - self.base
			self.end += bin_start - self.base
		except AttributeError: pass
class SSR(TRF):
	def __init__(self, genome, prefix='SR', rpt_type='tandem', 
				source='regex',
				definition='1-10 2-5 3-4 4-3 5-3 6-3'):
		self.genome = genome
		self.definition = [map(int, pair.split('-')) for pair in definition.split()]
		self.prefix = prefix
		self.rpt_type = rpt_type
		self.source = source
	def _parse(self):
		for rc in SeqIO.parse(self.genome, 'fasta'):
			seq = str(rc.seq.upper())
#			seq = seq[45954:45978]
			for unit_size,min_repeats in self.definition:
				pattern = r'([ATCG]{%s})(\1{%s,})' % (
							unit_size, min_repeats-1)
				for match in re.compile(pattern).finditer(seq):
					start = match.start() + 1
					end = match.end()
					unit = match.groups()[0]
					mseq = match.group()
					redundant = False
#					if unit_size == 5:
#						print >>sys.stderr, rc.id, start, end, unit, mseq
					for i in range(1, unit_size):
						redmotif = r'^([ATCG]{%s})(\1{%s})$' % (i, unit_size/i-1)
						if re.compile(redmotif).match(unit):
							redundant = True
							break
					if redundant:
						continue
					yield SSRRecord(rc.id, start, end, unit, mseq)
					
SSR_TYPE = {1:'mono', 2:'di', 3:'tri', 4:'tetra', 5:'penta', 6:'hexa'}
class SSRRecord(TRFRecord):
	def __init__(self, chrom, start, end, unit, seq):
		self.chrom, self.start, self.end = chrom, start, end
		self.unit, self.seq = unit, seq
		self.repeats = len(seq) / len(unit)
		self.ssr = '({}){}'.format(unit, self.repeats)
		self.type = SSR_TYPE[len(unit)]
		self.segment = RepeatSegemnt(self.chrom, self.start, self.end)
		self.segment.unit = self.unit
		self.segment.rpt = self.ssr

def main():
	args = makeArgparse()
	print >>sys.stderr, args.__dict__
	pipeline = RepeatPipeline(**args.__dict__)
	pipeline.run()

if __name__ == '__main__':
	main()
