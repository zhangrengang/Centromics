from .FileIO import Fastx

def subsample_seqs(seqfiles, fout, L=1e6, N=100, outfmt='fasta'):
	l_per_file = L / len(seqfiles) if L is not None else 0
	n_per_file = N // len(seqfiles) if N is not None else 0

	total_num, total_len = 0, 0
	for seqfile in seqfiles:
		_len = 0
		for i, rc in enumerate(Fastx(seqfile)):
			_len += len(rc)
			if (L is not None and _len < l_per_file) or (N is not None and i < n_per_file):
				rc.write(fout, outfmt)
				total_num += 1
				total_len += len(rc)
			else:
				break
	return total_num, total_len
