import sys, os, re
import argparse
import shutil
import multiprocessing
from xopen import xopen as open
from Bio import SeqIO
from REPcluster.Mcl import MclGroup
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp, test_s, get_suffix
from .RunCmdsMP import logger, run_cmd
from .sample_seqs import subsample_seqs
from .multi_seqs import multi_seqs
from .Trf import Trf, filter_trf_family, trf_map
from .Bin import bin_bam
from .Hic import count_obs
from .Peaks import find_peaks
from . import Circos
from .__version__ import version

NCPU = multiprocessing.cpu_count()
bindir = os.path.dirname(os.path.realpath(__file__))

def makeArgparse():
	parser = argparse.ArgumentParser( 
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Cluster Repeat Sequences.',
		)
	# input
	group_in = parser.add_argument_group('Input', )
	group_in.add_argument('-g', '-genome', default=None, metavar='FILE', 
					dest='genome', 
					help="Genome FASTA file")
	group_in.add_argument('-l', '-long', default=None, metavar='FILE', nargs='+', 
					required=True, dest='long',
					help="Long whole-genome-shotgun reads such as PacBio CCS/CLR or ONT reads \
in fastq or fasta format [required]")
	group_in.add_argument('-hic', default=None, metavar='FILE', 
					help="Hi-C data alignments by juicer")
	group_in.add_argument('-chip', default=None, metavar='FILE', 
					help="ChIP data alignments in bam format (sorted) or in BedGraph format")
	group_in.add_argument('-chip_input', default=None, metavar='FILE',
					help="Input data alignments in bam format (sorted)")	
	# output
	group_out = parser.add_argument_group('Output')
	group_out.add_argument('-pre', '-prefix', default='centomics', dest='prefix', metavar='STR',
					help="Prefix for output [default=%(default)s]")
	group_out.add_argument('-o', '-outdir', default='cent-output', dest='outdir', metavar='DIR',
					help="Output directory [default=%(default)s]")
	group_out.add_argument('-tmpdir', default='tmp', type=str, metavar='DIR',
					help="Temporary directory [default=%(default)s]")
	# tandem repeat
	group_tr = parser.add_argument_group('Kmer matrix',)
	group_tr.add_argument('-subsample_x', type=int, default=5, metavar='INT',
					 help="Subsample long reads up to X depth (prior to `-subsample_n`) [default=%(default)s]")
	group_tr.add_argument('-subsample_n', type=int, default=100000, metavar='INT',
					 help="Subsample long reads up to N reads [default=%(default)s]")
	group_tr.add_argument('-trf_opts', default='1 1 2 80 5 200 2000 -d -h', type=str,metavar='STR',
					help="TRF options to identify tandem repeats on a read [default='%(default)s']")
	group_tr.add_argument('-min_cov', type=float, default=0.9, metavar='FLOAT',
					 help="Minimum coverage of tandem repeats for a read [default=%(default)s]")
	group_tr.add_argument('-min_len', type=int, default=100, metavar='INT',
					 help="Minimum length of tandem repeats for a read [default=%(default)s]")
	group_tr.add_argument('-min_monomer_len', type=int, default=1, metavar='INT',
					 help="Minimum monomer length of a tandem repeat [default=%(default)s]")
	group_tr.add_argument('-clust_opts', default='-m jaccard -k 15 -c 0.2 -x 2 -I 2', type=str, metavar='STR',
					help="REPclust options to cluster tandem repeat units [default='%(default)s']")
	group_tr.add_argument('-min_ratio', type=float, default=0.1, metavar='FLOAT',
					 help="Minimum relative mass ratio to filter tandem repeats [default=%(default)s]")

	# circos
	group_circ = parser.add_argument_group('Circos', 'Options for circos plot')
	group_circ.add_argument('-window_size', type=int, default=200000, metavar='INT',
					help="Window size (bp) for circos plot [default=%(default)s]")
	group_circ.add_argument("-chr_prefix", metavar='STR',
					default='chr[\dXYZW]+',
					help='match chromosome to only plot chromosomes [default="%(default)s"]')
					
	# others
	group_other = parser.add_argument_group('Other options')
	group_other.add_argument('-p', '-ncpu', type=int, default=NCPU, metavar='INT', dest='ncpu',
					 help="Maximum number of processors to use [default=%(default)s]")
	group_other.add_argument('-cleanup', action="store_true", default=False,
					help="Remove the temporary directory [default=%(default)s]")	
	group_other.add_argument('-overwrite', action="store_true", default=False,
					help="Overwrite even if check point files existed [default=%(default)s]")
	group_other.add_argument('-v', '-version', action='version', version=version)
	
	args = parser.parse_args()
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	if args.prefix is not None:
		args.prefix = args.prefix.replace('/', '_')
	return args

		
class Pipeline:
	def __init__(self, **kargs):
		self.__dict__.update(**kargs)
		
	def run(self):
		# mkdir
		self._outdir = self.outdir = os.path.realpath(self.outdir)
		self._tmpdir = self.tmpdir = os.path.realpath(self.tmpdir)
		mkdirs(self.outdir)
		mkdirs(self.tmpdir)
		self.outdir += '/'
		self.tmpdir += '/'
		if not re.search(r'[\.\-_]$', self.prefix):
			self.prefix += '.'
		if self.prefix is not None:
			self.outdir = self.outdir + self.prefix
			self.tmpdir = self.tmpdir + self.prefix

		d_beds = {}
		# identify centromere tandem repeats
		tr_bed, tr_labels = self.run_long()
		d_beds['TR'] = tr_bed
		if not self.genome:
			return
		# hic
		hic_bed = self.run_hic()
		try: d_beds['HiC'] = hic_bed[0]
		except TypeError: pass
		# chip
		chip_bed = self.run_chip()
		d_beds['ChIP'] = chip_bed

		# out peak
		outbed = self.outdir + 'candidate_peaks.bed'
		with open(outbed, 'w') as fout:
			find_peaks(d_beds, fout, bin_size=self.window_size)

		# circos
		self.run_circos(tr_bed=tr_bed, tr_labels=tr_labels,
			hic_bed=hic_bed, # inter & intra
			chip_bed=chip_bed, 
			prefix=self.outdir + 'circos',
			chr_prefix=self.chr_prefix, window_size=self.window_size)
		
		self.step_final()
		logger.info('Pipeline completed.')
		return
		
	def run_circos(self, *args, **kargs):
		# circos
		circos_dir = bindir+'/circos'
		wkdir = self.outdir + 'circos' #self._outdir+'/circos'
		rmdirs(wkdir)
		try:
			logger.info('Copy `{}` to `{}`'.format(circos_dir, self._outdir))
			shutil.copytree(circos_dir, wkdir)
		except FileExistsError:
			pass
		Circos.centomics_plot(self.genome, wkdir, *args, **kargs)

	def get_trf_ids(self, trf_count):
		for line in open(trf_count):
			return line.strip().split()[3:]
	def run_long(self):
		logger.info('##Step: Processing long reads data')
		if self.genome is not None:
			logger.info('Loading {}'.format(self.genome))
			self.d_seqL = {rc.id: len(rc.seq) for rc in SeqIO.parse(open(self.genome), 'fasta')}
		
		trf_count = self.outdir + 'trf.count'
		ckp_file = self.mk_ckpfile(trf_count)
		if check_ckp(ckp_file, overwrite=self.overwrite):
			if not test_s(trf_count):
				tmpdir = '{}blast'.format(self.tmpdir)
				trf_famx = '{}/{}.blastqry'.format(tmpdir, self.prefix)
				with open(trf_count, 'w') as fout:
					trfids = trf_map(trf_famx, self.genome, fout, min_cov=0.8, ncpu=self.ncpu, window_size=10000, overwrite=0)
			return trf_count, self.get_trf_ids(trf_count)

#		if self.genome is not None:
#			logger.info('Loading {}'.format(self.genome))
#			self.d_seqL = {rc.id: len(rc.seq) for rc in SeqIO.parse(open(self.genome), 'fasta')}	
		# sample reads
		if self.genome is not None:
			genome_size = sum(self.d_seqL.values())
			logger.info('Genome size: {:,} bp'.format(genome_size))
			logger.info('Subsample {}x reads from {}'.format(self.subsample_x, self.long))
			L, N = self.subsample_x*genome_size, None
		else:
			L, N = None, self.subsample_n
			logger.info('Subsample {} reads from {}'.format(self.subsample_n, self.long))
			
		reads_fa = '{}sample.fa'.format(self.tmpdir)
		with open(reads_fa, 'w') as fout:
			total_num, total_len = subsample_seqs(self.long, fout, L=L, N=N)
			logger.info('Subsampled {:,} reads({:,} bases)'.format(total_num, total_len))
			
		# run trf
		logger.info('Run TRF to identify tandem repeats in reads')
		tmpdir = '{}trf'.format(self.tmpdir)
		trf = Trf(reads_fa, tmpdir=tmpdir, tr_opts=self.trf_opts, overwrite=self.overwrite, 
				nbins=self.ncpu*2, window_size=1e9, )
		trf_fa = '{}trf.fa'.format(self.tmpdir)
		with open(trf_fa, 'w') as fout:
			d_trf_len = trf.reads_trf(fout, min_cov=self.min_cov)
			
		# cluster
		logger.info('Cluster tandem repeats to identify TR families')
		tmpdir = '{}clust'.format(self.tmpdir)
		opts = '-pre {2} -outdir {0} -tmpdir {0} -p {1}'.format(tmpdir, self.ncpu, self.prefix)
		if self.overwrite:
			opts += ' -overwrite'
		cmd = 'REPclust {} {} {}'.format(trf_fa, self.clust_opts, opts)
		run_cmd(cmd, log=True, fail_exit=True)
		trfmcl = '{}/{}.mcl'.format(tmpdir, self.prefix)
		trfseq = '{}/{}.clust'.format(tmpdir, self.prefix)
		
		# filter trf family
		logger.info('Filter tandem repeats as putive centromeric')
		trf_fam = self.outdir + 'trf.fa'
		with open(trf_fam, 'w') as fout:
			filter_trf_family(trfseq, trfmcl, d_trf_len, fout, total_len=total_len, min_ratio=self.min_ratio)
			
		if self.genome is None:
			return trf_fam, None
		# align with genome and count density
		logger.info('Align with genome and count')
		tmpdir = '{}blast'.format(self.tmpdir)
		mkdirs(tmpdir)
		trf_famx = '{}/{}.blastqry'.format(tmpdir, self.prefix)
		with open(trf_famx, 'w') as fout:
			multi_seqs(seqfiles=[trf_fam], outfile=fout, 
						fold=1, min_length=50)
		with open(trf_count, 'w') as fout:
			trfids = trf_map(trf_famx, self.genome, fout, min_cov=0.8, ncpu=self.ncpu, window_size=10000)
		mk_ckp(ckp_file)
		return trf_count, trfids

	def run_chip(self):
		if not self.chip:
			return
		logger.info('##Step: Processing ChIP-seq data')
		if get_suffix(self.chip).lower() not in {'bam'}: #in {'.bdg', '.bedgraph'}:
			return self.chip

		chip_count = self.outdir + 'chip.bdg'
		ckp_file = self.mk_ckpfile(chip_count)
		if check_ckp(ckp_file, overwrite=self.overwrite):
			return chip_count
			
		bin_bam(self.chip, chip_count, input_bam=self.chip_input, ncpu=self.ncpu, bin_size=10000)
		mk_ckp(ckp_file)
		return chip_count
	def run_hic(self):
		if not self.hic:
			return
		logger.info('##Step: Processing Hi-C data')
		for res in (2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000):
			if res <= self.window_size:
				break
		hic_count = self.outdir + 'hic.count.{}'.format(res)
		
		ckp_file = self.mk_ckpfile(hic_count)
		ckp = check_ckp(ckp_file, overwrite=self.overwrite)
		if ckp:
			return ckp
			
		#with open(hic_count, 'w') as fout:
		#	count_links(self.hic, fout, bin_size=10000, ncpu=self.ncpu)
		
		tmpdir = '{}hic'.format(self.tmpdir)
		out1, out2 = count_obs(self.hic, chrLst=self.d_seqL.keys(), prefix=hic_count, tmpdir=tmpdir,
			bin_size=res, ncpu=self.ncpu, overwrite=self.overwrite)
		mk_ckp(ckp_file, out1, out2)
		return out1, out2
	def mk_ckpfile(self, file):
		return '{}{}.ok'.format(self.tmpdir, os.path.basename(file))
	def step_final(self):
		# cleanup
		if self.cleanup:
			logger.info('Cleaning {}'.format(self._tmpdir))
			rmdirs(self._tmpdir)
			
def main():
	args = makeArgparse()
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(version))
	logger.info('Arguments: {}'.format(args.__dict__))
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()

if __name__ == '__main__':
	main()
