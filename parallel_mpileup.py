#!/usr/bin/env python
import subprocess, sys, os, argparse, shutil

def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def main(cmdargs):
	
	if which('vcfutils.pl'):	pass
	else:	
		print 'Error: bcftools not found'
		sys.exit()

	if which('samtools'):	pass
	else:	
		print 'Error: samtools not found'
		sys.exit()

	if which('vcf-concat'):	pass
	else:	
		print 'Error: vcftools/PERL5LIB not found'
		sys.exit()

	if which('VarScan.v2.3.6.jar'): pass
	else:
		print 'Error: VarScan.v2.3.6.jar not found'
		sys.exit()

	args = ['vcfutils.pl', 'splitchr', '-l', cmdargs.chunk, cmdargs.fasta + '.fai']
	pro = subprocess.Popen(args, stdout=subprocess.PIPE)
	out, err = pro.communicate()
	parts = out.split('\n')
	parts = parts[:-1]

	makefile = open('Makefile', 'w')
	print >> makefile, 'all:',
	for each in range(len(parts)):
		print >> makefile, 'region-%s.vcf' %(each),
	print >> makefile, '\n'
	for i, each in enumerate(parts):
		if cmdargs.method == 'bcftools':
			print >> makefile, '''region-%s.vcf:
	samtools mpileup %s -DSIuf -r '%s' -b %s | bcftools view -cegv - > region-%s.vcf
	''' %(i, cmdargs.fasta, each, cmdargs.bamfiles, i)
		if cmdargs.method == 'varscan':
			print >> makefile, '''region-%s.vcf:
	samtools mpileup -f %s -r '%s' -b %s | java -jar bin/VarScan.v2.3.6.jar mpileup2snp --min-coverage 2 --min-var-freq 0.9 --p-value 0.005 --output-vcf --vcf-sample-list %s > region-%s.vcf
	''' %(i, cmdargs.fasta, each, cmdargs.bamfiles, cmdargs.bamfiles, i)
	makefile.close()

	args = ['make', '-j', cmdargs.threads]
	subprocess.call(args, stdout=subprocess.PIPE)
	
	list = []
	for each in range(len(parts)):
		list.append('region-%s.vcf' %(each))
	args = ['vcf-concat'] + list

	with open(cmdargs.output + '.vcf', 'w') as outfile:
		subprocess.call(args, stdout=outfile)

	# Clean up temporary files
	for each in list:
		os.remove(each)
	os.remove('Makefile')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script is designed to speed up samtools mpileup by running a subregion per thread then merging the output.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('fasta', help="Reference fasta")
	parser.add_argument('bamfiles', help="File, list of input BAM files")
	parser.add_argument('method', choices=['bcftools', 'varscan'], help="Default settings are: 'bcftools view -cegv -' and 'VarScan.v2.3.6.jar mpileup2snp --min-coverage 2 --min-var-freq 0.9 --p-value 0.005'")
	parser.add_argument('-t', '--threads', metavar='', help="Number of threads", default='32')
	parser.add_argument('-c', '--chunk', metavar='', help="Chunk size in bases", default='250000')
	parser.add_argument('-o', '--output', metavar='', help="Output suffix", default='merged')
	cmdargs = parser.parse_args()
	main(cmdargs)
