import sys
import math
from Bio import SeqIO
from Bio.Seq import Seq
from xopen import xopen as open

def main():
	seqfile = sys.argv[1]
	try: fold = int(sys.argv[2])
	except KeyError: fold = 100
	outfile = sys.stdout

def multi_seqs(seqfiles=None, seqlist=None, outfile=sys.stdout, fold=2, min_length=100):
	if not seqfiles:
		seqfiles = []
	if seqlist is not None:
		seqfiles += list([line.strip().split()[0] for line in open(seqlist)])
	d_seqs = {}
	for seqfile in seqfiles:
		for rc in SeqIO.parse(seqfile, 'fasta'):
			d_seqs[rc.id] = rc.seq
			if len(rc.seq) < min_length/fold:
				_fold = int(math.ceil(1.0*min_length/len(rc.seq)))
				fold = max(fold, _fold)
			rc.seq = Seq(''.join([str(rc.seq)] * fold))
			rc.description = rc.id
			SeqIO.write(rc, outfile, 'fasta')
	return d_seqs
