import sys,os
import copy
import argparse
import shutil
import json
import math
import multiprocessing
from collections import OrderedDict, Counter
from xopen import xopen as open
from Bio import SeqIO
from REPcluster.Mcl import MclGroup
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp, test_s
from .RunCmdsMP import logger, run_cmd
from .sample_seqs import subsample_seqs
from .multi_seqs import multi_seqs
from .Trf import Trf, filter_trf_family, trf_map
from .Bin import bin_bam
from .Hic import count_links
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
					help="ChIP data alignments in bam format (sorted)")
					
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
	group_tr.add_argument('-subsample_x', type=int, default=10, metavar='INT',
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
					default='chr\d+',
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
	#	if not self.prefix.endswith('.'):
	#		self.prefix += '.'
		if self.prefix is not None:
			self.outdir = self.outdir + self.prefix
			self.tmpdir = self.tmpdir + self.prefix

		# identify centromere tandem repeats
		
		tr_bed = self.run_long()
		if not self.genome:
			return
		# hic
		hic_bed = self.run_hic()
		# chip
		chip_bed = self.run_chip()
		
		self.run_circos(tr_bed=tr_bed, hic_bed=hic_bed, chip_bed=chip_bed, 
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

	def run_long(self):
		logger.info('##Step: Processing long reads data')
		trf_count = self.outdir + 'trf.count'
		ckp_file = self.mk_ckpfile(trf_count)
		if check_ckp(ckp_file) and not self.overwrite:
			return trf_count
		
		# sample reads
		if self.genome is not None:
			genome_size = sum([len(rc.seq) for rc in SeqIO.parse(self.genome, 'fasta')])
			logger.info('Subsample {}x reads from {}'.format(self.subsample_x, self.long))
			L, N = self.subsample_x*genome_size, None
		else:
			L, N = None, self.subsample_n
			logger.info('Subsample {} reads from {}'.format(self.subsample_n, self.long))
			
		reads_fa = '{}sample.fa'.format(self.tmpdir)
		with open(reads_fa, 'w') as fout:
			total_num, total_len = subsample_seqs(self.long, fout, L=L, N=N)
			
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
		run_cmd(cmd, log=True)
		trfmcl = '{}/{}.mcl'.format(tmpdir, self.prefix)
		trfseq = '{}/{}.clust'.format(tmpdir, self.prefix)
		
		# filter trf family
		logger.info('Filter tandem repeats as putive centromeric')
		trf_fam = self.outdir + 'trf.fa'
		with open(trf_fam, 'w') as fout:
			filter_trf_family(trfseq, trfmcl, d_trf_len, fout, total_len=total_len, min_ratio=self.min_ratio)
			
		if self.genome is None:
			return trf_fam
		# align with genome and count density
		logger.info('Align with genome and count')
		tmpdir = '{}blast'.format(self.tmpdir)
		mkdirs(tmpdir)
		trf_famx = '{}/{}.blastqry'.format(tmpdir, self.prefix)
		with open(trf_famx, 'w') as fout:
			multi_seqs(seqfiles=[trf_fam], outfile=fout, 
						fold=1, min_length=50)
		with open(trf_count, 'w') as fout:
			trf_map(trf_famx, self.genome, fout, min_cov=0.8, ncpu=self.ncpu, window_size=10000)
		mk_ckp(ckp_file)
		return trf_count
	def run_chip(self):
		if not self.chip:
			return
		logger.info('##Step: Processing ChIP-seq data')
		chip_count = self.outdir + 'chip.count'
		ckp_file = self.mk_ckpfile(chip_count)
		if check_ckp(ckp_file) and not self.overwrite:
			return chip_count
			
		bin_bam(self.chip, chip_count, ncpu=self.ncpu, bin_size=10000)
		mk_ckp(ckp_file)
		return chip_count
	def run_hic(self):
		if not self.hic:
			return
		logger.info('##Step: Processing Hi-C data')
		hic_count = self.outdir + 'hic.count'
		ckp_file = self.mk_ckpfile(hic_count)
		if check_ckp(ckp_file) and not self.overwrite:
			return hic_count
			
		with open(hic_count, 'w') as fout:
			count_links(self.hic, fout, bin_size=10000, ncpu=self.ncpu)
		mk_ckp(ckp_file)
		return hic_count
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